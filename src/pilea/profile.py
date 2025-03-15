import sys
import re
import os
import glob
import pickle
import numpy as np

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from tqdm.contrib.logging import logging_redirect_tqdm
from collections import defaultdict
from functools import partial
from itertools import chain
from sklearn.linear_model import RANSACRegressor
from scipy.stats import median_abs_deviation

from .log import log
from .kmc import KMC
from .ztp import ZTP


class GrowthProfiler:
    def __init__(self, files, outdir, database, force=False, single=False, threads=os.cpu_count()):
        '''
        Profile bacterial growth dynamics.
        '''
        self.force = force
        self.outdir = outdir
        self.threads = threads
        self.database = database

        ## sanity check
        if os.path.isdir(outdir):
            log.warning(f'Output directory <{outdir}> exists.')
        else:
            os.makedirs(outdir, exist_ok=True)

        if not os.path.isdir(database):
            log.critical(f'Database directory <{database}> does not exist.')
            sys.exit(2)
        else:
            for file in ['parameter.tab', 'taxonomy.tab', 'sketch.uni', 'sketch.dup', 'sketch.kmc_pre', 'sketch.kmc_suf']:
                if not os.path.isfile(f'{database}/{file}'):
                    log.critical(f'File <{file}> is missing from database directory <{database}>.')
                    sys.exit(2)

            for file in files:
                if not os.path.isfile(file):
                    log.critical(f'File <{file}> does not exist.')
                    sys.exit(2)

        ## determine file format
        samples = [re.sub('.gz$', '', os.path.basename(file)) for file in sorted(files)]
        extensions = {sample.split('.')[-1] for sample in samples}
        if len(extensions) != 1 or not {'fasta', 'fa', 'fastq', 'fq'} & extensions:
            log.critical('Input files need to end with <fa|fq|fasta|fastq>.')
            sys.exit(2)
        extension = extensions.pop()

        if not single:
            samples = [re.sub(rf'(_(1|2|R1|R2|fwd|rev))?.{extension}$', '', sample) for sample in samples]
        else:
            samples = [re.sub(rf'.{extension}$', '', sample) for sample in samples]

        self.files = defaultdict(list)
        for sample, file in zip(samples, files):
            self.files[sample].append(file)

        if len({len(file) for file in self.files.values()}) != 1:
            log.warning('Files are mixed with single/paired-end. Check whether <--single> is needed.')

        ## auto detect format for KMC
        self.f = 'fm' if extension in {'fa', 'fasta'} else 'fq'
        with open(f'{database}/parameter.tab') as f:
            for line in f:
                if line[0] == 'k':
                    self.k = int(line.split()[-1])

    def count(self):
        '''
        Run KMC for k-mer counting.
        '''
        log.info('Counting k-mers ...')
        queue = tqdm(self.files.items(), leave=False)
        with logging_redirect_tqdm():
            for sample, file in queue:
                queue.set_description(f'==> Processing <{sample}>')
                if os.path.isfile(f'{self.outdir}/{sample}.sketch.txt') and not self.force:
                    log.info(f'File <{self.outdir}/{sample}.sketch.txt> exists, skip. Use <--force> for overwritting.')
                else:
                    kmc = KMC(file=file, prefix=f'{self.outdir}/{sample}', threads=self.threads)
                    kmc.count(k=self.k, f=self.f, ci=1)
                    kmc.intersect(in_prefix=f'{self.outdir}/{sample}', ref_prefix=f'{self.database}/sketch', out_prefix=f'{self.outdir}/{sample}.sketch')
                    kmc.dump(f'{self.outdir}/{sample}.sketch')
                for file in glob.glob(f'{self.outdir}/{sample}.*kmc_*'):
                    os.remove(file)

    @staticmethod
    def _parse(kval, min_cont, accession2info, kmer2accession):
        sample, kcnt = kval
        kset = defaultdict(set)
        for key, val in kcnt.items():
            for accession in kmer2accession.get(key):
                kset[accession.split('|', 1)[0]].add(key)
        kset = {key: val for key, val in kset.items() if len(val) / accession2info.get(key)[-1] > min_cont}

        data = []
        while kset:
            ## pop an accession with the highest containment
            accession = max(kset, key=lambda x: len(kset[x]) / accession2info.get(x)[-1])

            ## record kmers' information
            vset = kset.get(accession)
            data.append([sample, *accession2info.get(accession), [(z[1], kcnt.get(x)) for x in vset for y in kmer2accession.get(x) if (z := y.split('|', 1))[0] == accession]])

            ## update by set subtraction
            kset = {key: newval for key, val in kset.items() if (len(newval := val - vset) / accession2info.get(key)[-1]) > min_cont}
        return data

    @staticmethod
    def _filter(row, max_disp, min_dept, min_cont):
        def _trim(x, return_limits=False):
            if x.size == 0:
                return x
            x = np.log2(x)
            q1, q3 = np.percentile(x, [25, 75])
            lower, upper = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)

            x = np.rint(2 ** x[(x >= lower) & (x <= upper)])
            return (2 ** lower, 2 ** upper) if return_limits else x

        kcnt = defaultdict(list)
        for x in row[-1]:
            if x[0][-1] == '+':
                kcnt[x[0].rsplit('|', 1)[0]].append(x[1])

        lower, upper = _trim(np.asarray(list(chain(*kcnt.values()))), return_limits=True)
        depths, dispersions, observations = [], [], []
        for key, val in kcnt.items():
            val = np.asarray(val)
            val = _trim(val[(val >= lower) & (val <= upper)])
            if len(val) > max(1, int(key.rsplit('|', 1)[-1]) / 2):
                mean, var = np.mean(val), np.var(val, ddof=1)
                depths.append(mean)
                dispersions.append(var / mean)
                observations.append(val)

        if observations:
            if (
                (depth := np.median(depths)) > min_dept and
                (dispersion := np.median(dispersions)) < max_disp and
                (containment := sum(len(x) for x in observations) / row[3]) > min_cont
            ):
                return row[:3] + [depth, dispersion, containment, observations]

    def parse(self, max_disp=25, min_dept=5, min_cont=0.5):
        '''
        Parse and filter outputs of KMC.
        '''
        ## get taxonomy/accession mapping for observed kmers
        log.info('Loading mapping files ...')
        kval = defaultdict(dict)
        for sample in self.files.keys():
            file = f'{self.outdir}/{sample}.sketch.txt'
            with open(file) as f:
                for line in f:
                    kval[sample][line[:self.k]] = int(line[self.k + 1:-1])

        accession2info = dict()
        with open(f'{self.database}/taxonomy.tab') as f:
            for line in f:
                ls = line.rstrip().split('\t')
                accession2info[ls[0]] = ls[1], ls[2], int(ls[3])

        kmer2accession = dict()
        aset = {x for key, val in kval.items() for x in val.keys()}
        with open(f'{self.database}/sketch.uni') as f:
            for line in f:
                if line[:self.k] in aset:
                    kmer2accession[line[:self.k]] = [line[self.k + 1:-1]]

        with open(f'{self.database}/sketch.dup') as f:
            for line in f:
                if line[:self.k] in aset:
                    kmer2accession[line[:self.k]] = line[self.k + 1:-1].split(',')

        log.info('Parsing and filtering outputs ...')
        ## allocate kmers to species based on containment ANI
        fun = partial(self._parse, min_cont=min_cont, accession2info=accession2info, kmer2accession=kmer2accession)
        self.data = [row for sub in process_map(fun, kval.items(), max_workers=self.threads, chunksize=1, leave=False) for row in sub]

        ## prune local/global outliers
        fun = partial(self._filter, max_disp=max_disp, min_dept=min_dept, min_cont=min_cont)
        self.data = [row for row in process_map(fun, self.data, max_workers=self.threads, chunksize=1, leave=False) if row]

    @staticmethod
    def _fit(observation, components, max_iter, tol):
        Y = []
        for x in observation:
            r = []
            for n in range(1, components + 1):
                r.append(ZTP(x=x).fit(components=n, max_iter=max_iter, tol=tol))
            lmds, weights, _ = sorted(r, key=lambda x: x[-1])[0]
            Y.append(lmds[np.argmax(weights)])

        Y = np.log2(sorted(Y))
        X = np.asarray(range(len(Y))).reshape(-1, 1)
        return [2 ** (RANSACRegressor(residual_threshold = median_abs_deviation(Y) / 5, random_state=i).fit(X, Y).estimator_.coef_[0] * (len(Y) - 1)) for i in range(100)]

    def fit(self, components=5, max_iter=np.inf, tol=1e-5):
        log.info('Fitting counts ...')
        fun = partial(self._fit, components=components, max_iter=max_iter, tol=tol)
        for row, ptr in zip(self.data, process_map(fun, [row[-1] for row in self.data], max_workers=self.threads, chunksize=1, leave=False)):
            row[-1] = np.median(ptr) if median_abs_deviation(ptr) < 0.5 else None

    def write(self):
        with open(f'{self.outdir}/output.tsv', 'w') as f:
            f.write('\t'.join(['sample', 'genome', 'taxonomy', 'depth', 'dispersion', 'containment', 'ptr']) + '\n')
            for row in sorted(self.data):
                if row[-1] is not None:
                    f.write('\t'.join(row[:3] + [f'{x:.4f}' for x in row[3:]]) + '\n')
        log.info('Done.')

def profile(files, outdir, database, force=False, single=False, max_disp=25, min_dept=5, min_cont=0.5, components=5, max_iter=np.inf, tol=1e-5, threads=os.cpu_count()):
    '''
    Profile bacterial growth dynamics.
    '''
    gp = GrowthProfiler(files=files, outdir=outdir, database=database, force=force, single=single, threads=threads)
    gp.count()
    gp.parse(max_disp=max_disp, min_dept=min_dept, min_cont=min_cont)
    gp.fit(components=components, max_iter=max_iter, tol=tol)
    gp.write()