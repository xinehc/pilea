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
from scipy.stats import median_abs_deviation, iqr
from statsmodels.nonparametric.smoothers_lowess import lowess

from .log import log
from .kmc import KMC
from .ztp import ZTP


class GrowthProfiler:
    def __init__(self, files, outdir, database, force=False, single=False, threads=os.cpu_count()):
        '''
        Profile bacterial growth dynamics.
        '''
        self.data = []
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
    def _filter(row, min_dept, max_disp, min_frac):
        def _trim(x, return_limits=False):
            if x.size == 0:
                return x
            x = np.log2(x)
            q1, q3 = np.percentile(x, [25, 75])
            lower, upper = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)

            x = np.rint(2 ** x[(x >= lower) & (x <= upper)]).astype(np.int32)
            return (2 ** lower, 2 ** upper) if return_limits else x

        def _debias(x, y, frac=0.25):
            t = lowess(exog=x, endog=y, frac=frac)[:, 1]
            y = np.rint(2 ** (np.asarray(y) - t + np.mean(t))).astype(np.int32)
            return y

        A = []
        for val in row[-1]:
            s = val[0].rsplit('|', 2)
            A.append((int(s[1]), np.log2(val[1]), s[0]))

        ## sort sketches by their GC content
        A = sorted(A, key=lambda x: x[0])
        Y = _debias(x=[a[0] for a in A], y=[a[1] for a in A])

        kcnt = defaultdict(list)
        for a, y in zip(A, Y):
            if y != 0:
                kcnt[a[-1]].append(y)

        ## trim global/local outliers
        lower, upper = _trim(np.asarray(list(chain(*kcnt.values()))), return_limits=True)
        depths, dispersions, observations = [], [], []
        for key, val in kcnt.items():
            val = np.asarray(val)
            val = _trim(val[(val >= lower) & (val <= upper)])
            if len(val) > max(1, int(key.rsplit('|', 1)[-1]) / 5):
                mean, var = np.mean(val), np.var(val, ddof=1)
                depths.append(mean)
                dispersions.append(var / mean)
                observations.append(val)

        if observations:
            if (
                (depth := np.median(depths)) > min_dept and
                (dispersion := np.median(dispersions)) < max_disp and
                (fraction := len(observations) / row[3]) > min_frac
            ):
                return row[:3] + [depth, dispersion, fraction] + [row[4], observations]

    def parse(self, min_dept=5, max_disp=np.inf, min_frac=0.75, min_cont=0.75):
        '''
        Parse and filter outputs of KMC.
        '''
        log.info('Loading mapping files ...')
        info = dict()
        with open(f'{self.database}/taxonomy.tab') as f:
            for i, line in enumerate(f.readlines()):
                ls = line.rstrip().split('\t')
                info[str(i)] = ls[0], ls[4], int(ls[2]), int(ls[3])

        ## load observed counts
        obs = defaultdict(dict)
        for sample in self.files.keys():
            file = f'{self.outdir}/{sample}.sketch.txt'
            with open(file) as f:
                for line in f:
                    obs[sample][line[:self.k]] = int(line[self.k + 1:-1])

        ## load sketches' mapping files
        kuni, kdup = dict(), dict()
        kset = {x for key, val in obs.items() for x in val.keys()}
        with open(f'{self.database}/sketch.uni') as f, open(f'{self.database}/sketch.dup') as g:
            for line in f:
                if line[:self.k] in kset:
                    kuni[line[:self.k]] = line[self.k + 1:-1]

            for line in g:
                if line[:self.k] in kset:
                    kdup[line[:self.k]] = line[self.k + 1:-1].split(',')

        log.info('Parsing and filtering outputs ...')
        queue = tqdm(obs.items(), leave=False)
        for sample, kcnt in queue:
            queue.set_description(f'==> Processing <{sample}>')

            ## assign unique sketches directly
            ku, kd = defaultdict(list), defaultdict(set)
            for key, val in kcnt.items():
                if (i := kuni.get(key)):
                    i = i.split('|', 1)
                    ku[i[0]].append((i[1], val))
                else:
                    for i in kdup[key]:
                        kd[i.split('|', 1)[0]].add(key)

            ## assign shared sketches based on containment
            ka = dict()
            kc = {key: len(val) / info[key][-1] for key, val in ku.items()}
            while kd:
                accession = max(kd, key = lambda key: kc.get(key, 0) + len(kd[key]) / info[key][-1])
                ks = kd.pop(accession)
                kd = {key: nval for key, val in kd.items() if kc.get(key, 0) + len(nval := val - ks) / info[key][-1] > min_cont / 2}

                kc[accession] = kc.get(accession, 0) + len(ks) / info[accession][-1]
                if kc[accession] > min_cont:
                    ka[accession] = ks

            for accession, ks in ka.items():
                for idx, val in [(kdup[x], kcnt[x]) for x in sorted(ks)]:
                    if len({s for i in idx if kc.get(s := i.split('|', 1)[0], 0) > min_cont / 2}) == 1:
                        for i in idx:
                            s = i.split('|', 1)
                            if s[0] == accession and s[1][-1] == '+':
                                ku[s[0]].append((s[1], val))

            self.data.extend([[sample, *info[key][:-1], cont, val] for key, val in ku.items() if (cont := kc[key]) > min_cont])

        ## prune local/global outliers
        fun = partial(self._filter, max_disp=max_disp, min_dept=min_dept, min_frac=min_frac)
        self.data = [row for row in process_map(fun, self.data, max_workers=self.threads, chunksize=1, leave=False) if row]

    @staticmethod
    def _fit(observation, components, max_iter, tol):
        def _sample(u, seed):
            rng = np.random.RandomState(seed)
            return np.log2(sorted(rng.choice(y, p=p) for y, p in u))

        u, v = [], []
        for x in observation:
            fits = [ZTP(x=x).fit(components=n, max_iter=max_iter, tol=tol) for n in range(1, components + 1)]
            lmds, weights, _ = sorted(fits, key=lambda x: x[-1])[0]
            u.append((lmds, weights))

        N = len(u)
        X = np.asarray(range(N)).reshape(-1, 1)
        for i in range(100):
            Y = _sample(u, seed=i)
            R = median_abs_deviation(Y) / 5
            v.append(2 ** (RANSACRegressor(residual_threshold=R, random_state=i).fit(X, Y).estimator_.coef_[0] * (N - 1)))

        if iqr(v, rng=(10, 90)) < 1:
            Y = np.log2(sorted(y[np.argmax(p)] for y, p in u))
            R = median_abs_deviation(Y) / 5
            return np.median([2 ** (RANSACRegressor(residual_threshold=R, random_state=i).fit(X, Y).estimator_.coef_[0] * (N - 1)) for i in range(100)])

    def infer(self, components=5, max_iter=np.inf, tol=1e-5):
        log.info('Fitting counts ...')
        fun = partial(self._fit, components=components, max_iter=max_iter, tol=tol)
        for row, ptr in zip(self.data, process_map(fun, [row[-1] for row in self.data], max_workers=self.threads, chunksize=1, leave=False)):
            row[-1] = ptr

    def write(self):
        with open(f'{self.outdir}/output.tsv', 'w') as f:
            f.write('\t'.join(['sample', 'genome', 'taxonomy', 'depth', 'dispersion', 'fraction', 'containment', 'ptr']) + '\n')
            for row in sorted(self.data, key=lambda row: (row[0], row[2], row[1])):
                if row[-1] is not None:
                    f.write('\t'.join(row[:3] + [f'{x:.4f}' for x in row[3:]]) + '\n')
        log.info('Done.')

def profile(files, outdir, database, force=False, single=False, min_dept=5, max_disp=np.inf, min_frac=0.75, min_cont=0.75, components=5, max_iter=np.inf, tol=1e-5, threads=os.cpu_count()):
    '''
    Profile bacterial growth dynamics.
    '''
    gp = GrowthProfiler(files=files, outdir=outdir, database=database, force=force, single=single, threads=threads)
    gp.count()
    gp.parse(max_disp=max_disp, min_dept=min_dept, min_frac=min_frac, min_cont=min_cont)
    gp.infer(components=components, max_iter=max_iter, tol=tol)
    gp.write()