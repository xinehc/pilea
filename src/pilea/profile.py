import sys
import re
import os
import mmap
import glob
import numpy as np
import multiprocessing as mp

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from tqdm.contrib.logging import logging_redirect_tqdm
from collections import defaultdict
from functools import partial
from itertools import chain
from sklearn.linear_model import RANSACRegressor
from scipy.stats import median_abs_deviation, iqr
from statsmodels.nonparametric.smoothers_lowess import lowess
from needletail import parse_fastx_file
import packaging

from .log import log
from .ztp import ZTP
from .utils import u8, u3
from .kmc import count64

from . import __version__


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
            for file in ['parameters.tab', 'genomes.tab', 'sketches.pdb']:
                if not os.path.isfile(f'{database}/{file}'):
                    log.critical(f'File <{file}> is missing from database directory <{database}>.')
                    sys.exit(2)

            for file in files:
                if not os.path.isfile(file):
                    log.critical(f'File <{file}> does not exist.')
                    sys.exit(2)

        ## determine file format
        samples = [re.sub('.gz$', '', os.path.basename(file)) for file in files]
        extensions = {sample.split('.')[-1] for sample in samples}
        if len(extensions) != 1 or not {'fasta', 'fa', 'fastq', 'fq'} & extensions:
            log.critical('Input files need to end with <fa|fq|fasta|fastq>.')
            sys.exit(2)
        extension = extensions.pop()

        if not single:
            samples = [re.sub(rf'(_(1|2|R1|R2|fwd|rev))?.{extension}$', '', sample) for sample in samples]
        else:
            samples = [re.sub(rf'.{extension}$', '', sample) for sample in samples]

        self.items = defaultdict(list)
        for sample, file in sorted(zip(samples, files)):
            self.items[sample].append(file)

        for sample, files in self.items.items():
            if len(files) > 2:
                log.critical(f'Sample <{sample}> has more than two files. Check file format.')
                sys.exit(2)

        if len({len(file) for file in self.items.values()}) != 1:
            log.warning('Files are mixed with single-/paired-end. Check whether <--single> is needed.')

        ## load database parameters
        with open(f'{database}/parameters.tab') as f:
            for line in f:
                if line[0] == 'k':
                    self.k = int(line.split()[-1])
                elif line[0] == 's':
                    self.s = int(line.split()[-1])
                    self.m = ((1 << 64) - 1) // self.s
                elif line[0] == 'w':
                    self.w = int(line.split()[-1])
                elif line[0] == 'v':
                    self.v = packaging.version.parse(line.split()[-1])

        if (v := packaging.version.parse(__version__)) < self.v:
            log.critical(f'Database <{database}> requires <v{self.v}> or above.')
            sys.exit(2)

        self.meta = dict()
        with open(f'{self.database}/genomes.tab') as f:
            for i, line in enumerate(f.readlines()):
                sp = line.rstrip().split('\t')
                self.meta[int(i)] = sp[0], sp[4], int(sp[3]), int(sp[2])

    @staticmethod
    def _collect(items, k, m, outdir, force, REC=11):
        sample, files = items
        kmc = f'{outdir}/{sample}.kmc'
        if os.path.isfile(kmc) and not force:
            with open(kmc, 'rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                n = mm.size() // REC
                keys = np.ndarray(shape=(n,), dtype='>u8', buffer=mm, offset=0, strides=(REC,))
                cnts = np.ndarray(shape=(n, 3), dtype='>u1',  buffer=mm, offset=8, strides=(REC, 1)).astype(np.uint32, copy=False)
                counts = {int(k): int(c) for k, c in zip(keys, (cnts[:, 0] << 16) | (cnts[:, 1] << 8) | cnts[:, 2])}
        else:
            counts = {}
            if len(files) == 1:
                it = parse_fastx_file(files[0])
                for record in it:
                    count64(record.seq, '', k, m, counts)
            else:
                it1 = parse_fastx_file(files[0])
                it2 = parse_fastx_file(files[1])
                for record1, record2 in zip(it1, it2):
                    count64(record1.seq, record2.seq, k, m, counts)
        return sample, counts

    @staticmethod
    def _load(file, arr, REC=18, CHUNK_RECS=1000000):
        pdb = defaultdict(list)
        if arr.size == 0:
            return pdb

        with open(file, 'rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            nrec = mm.size() // REC
            for i in tqdm(range(0, nrec, CHUNK_RECS), leave=False):
                n = min(CHUNK_RECS, nrec - i)
                keys = np.ndarray(shape=(n,), dtype='>u8', buffer=mm, offset=i * REC, strides=(REC,))

                pos  = arr.searchsorted(keys)
                ins  = np.minimum(pos, arr.size - 1)
                mask = (pos < arr.size) & (arr[ins] == keys)
                for j in np.flatnonzero(mask).tolist():
                    o = (i + j) * REC + 8

                    idx = (mm[o + 0] << 16) | (mm[o + 1] << 8) | mm[o + 2]
                    widx = (mm[o + 3] << 8) | mm[o + 4]
                    sins = (mm[o + 5] << 8) | mm[o + 6]
                    gc = (mm[o + 7] << 8) | mm[o + 8]
                    uq = (mm[o + 9] != 0)

                    pdb[int(keys[j])].append((idx, widx, sins, gc, uq))
        return pdb

    @staticmethod
    def _filter(row, min_cove, max_disp, min_frac, min_cont):
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
            A.append((val[0][-1], np.log2(val[1]), val[0][:2]))

        ## sort kmers by their GC content
        A = sorted(A, key=lambda x: x[0])
        Y = _debias(x=[a[0] for a in A], y=[a[1] for a in A])

        kcnt = defaultdict(list)
        for a, y in zip(A, Y):
            if y != 0:
                kcnt[a[-1]].append(y)

        ## trim global/local outliers
        lower, upper = _trim(np.asarray(list(chain(*kcnt.values()))), return_limits=True)
        coverages, dispersions, observations = [], [], []
        for key, val in kcnt.items():
            val = np.asarray(val)
            val = _trim(val[(val >= lower) & (val <= upper)])
            if len(val) > max(1, key[-1] * min_cont):
                mean, var = np.mean(val), np.var(val, ddof=1)
                coverages.append(mean)
                dispersions.append(var / mean)
                observations.append(val)

        if observations:
            if (
                (coverage := np.median(coverages)) > min_cove and
                (dispersion := np.median(dispersions)) < max_disp and
                (fraction := len(observations) / row[3]) > min_frac
            ):
                return row[:3] + [coverage, dispersion, fraction] + [row[4], observations]

    @staticmethod
    def _fit(observation, components, max_iter, tol):
        def _sample(u, seed):
            rng = np.random.RandomState(seed)
            return np.log2(sorted(rng.choice(y, p=p) for y, p in u))

        def _ransac(x, y, r, seed):
            fit = RANSACRegressor(residual_threshold=r, random_state=seed).fit(x, y)
            return 2 ** (fit.estimator_.coef_[0] * (len(y) - 1))

        u, v, w = [], [], []
        for x in observation:
            fits = [ZTP(x=x).fit(components=n, max_iter=max_iter, tol=tol) for n in range(1, components + 1)]
            lmds, weights, _ = sorted(fits, key=lambda x: x[-1])[0]
            u.append((lmds, weights))

        X = np.arange(len(u)).reshape(-1, 1)
        for i in range(100):
            Y = _sample(u, seed=i)
            R = median_abs_deviation(Y) / 5
            v.append(_ransac(X, Y, R, seed=i))

        if iqr(v, rng=(10, 90)) < 1:
            Y = np.log2(sorted(y[np.argmax(p)] for y, p in u))
            R = median_abs_deviation(Y) / 5
            for i in range(100):
                w.append(_ransac(X, Y, R, seed=i))
            return np.median(w)

    def count(self, min_cove=5, max_disp=np.inf, min_frac=0.75, min_cont=0.25):
        '''
        Count, parse and filter.
        '''
        log.info('Counting k-mers ...')
        kmc = [file for file in glob.glob(f'{self.outdir}/*.kmc') if os.path.basename(file).split('.kmc')[0] in self.items]
        if kmc and not self.force:
            log.info(f"File {' '.join(f'<{os.path.basename(file)}>' for file in kmc)} exists, skip. Use <--force> for overwritting.")

        arr = set()
        obs = dict()
        fun = partial(self._collect, k=self.k, m=self.m, outdir=self.outdir, force=self.force)
        with logging_redirect_tqdm():
            with mp.Pool(self.threads) as pool:
                for sample, counts in tqdm(pool.imap_unordered(fun, self.items.items(), chunksize=1), total=len(self.items), leave=False):
                    obs[sample] = counts
                    arr.update(counts.keys())

        arr = np.fromiter(arr, dtype='>u8', count=len(arr))
        arr.sort()

        log.info('Loading database ...')
        pdb = self._load(f'{self.database}/sketches.pdb', arr)

        log.info('Parsing and filtering outputs ...')
        for sample, kmc in tqdm(obs.items(), leave=False):
            kmc = {k: v for k, v in kmc.items() if k in pdb}  # recast
            if not os.path.isfile(f'{self.outdir}/{sample}.kmc') or self.force:
                with open(f'{self.outdir}/{sample}.kmc', 'wb') as f:
                    for key, cnt in kmc.items():
                        f.write(u8(key))
                        f.write(u3(cnt))

            ku, kd = defaultdict(list), defaultdict(set)
            for key, cnt in kmc.items():
                if len(pos := pdb.get(key)) == 1:  # unique
                    ku[pos[0][0]].append((pos[0][1:4], cnt))
                else:
                    for idx, *_ in pos:
                        kd[idx].add(key)

            ka = {}
            kc = {idx: len(val) / self.meta[idx][-1] for idx, val in ku.items()}
            for idx in list(kd):
                if kc.get(idx, 0) + len(kd[idx]) / self.meta[idx][-1] <= min_cont:
                    del kd[idx]

            while kd:
                ba = max(kd, key=lambda idx: kc.get(idx, 0) + len(kd[idx]) / self.meta[idx][-1])
                bs = kd.pop(ba)

                sa = set()
                for s in bs:
                    for idx, *_ in pdb[s]:
                        sa.add(idx)  # relevant idx

                for idx in sa & kd.keys():
                    kd[idx] -= bs
                    if kc.get(idx, 0) + len(kd[idx]) / self.meta[idx][-1] <= min_cont:
                        del kd[idx]

                kc[ba] = kc.get(ba, 0) + len(bs) / self.meta[ba][-1]
                ka[ba] = sorted(bs)

            for ba, bs in ka.items():
                for s in bs:
                    cand = [pos for pos in pdb[s] if pos[0] in ka]
                    if len(cand) == 1 and cand[0][4]:  # exactly one candidate and flag is +
                        ku[ba].append((cand[0][1:4], kmc[s]))

            self.data.extend([[sample, *self.meta[idx][:-1], cont, val] for idx, val in ku.items() if (cont := len(val) / self.meta[idx][-1]) > min_cont])

        ## prune local/global outliers
        fun = partial(self._filter, min_cove=min_cove, max_disp=max_disp, min_frac=min_frac, min_cont=min_cont)
        self.data = [row for row in process_map(fun, self.data, max_workers=self.threads, chunksize=1, leave=False) if row]

    def infer(self, components=5, max_iter=np.inf, tol=1e-5):
        '''
        Estiamte PTRs by fitting counts with ZTP.
        '''
        log.info('Fitting counts ...')
        fun = partial(self._fit, components=components, max_iter=max_iter, tol=tol)
        for row, ptr in zip(self.data, process_map(fun, [row[-1] for row in self.data], max_workers=self.threads, chunksize=1, leave=False)):
            row[-1] = ptr

    def write(self):
        '''
        Write valid PTRs to output.
        '''
        with open(f'{self.outdir}/output.tsv', 'w') as f:
            f.write('\t'.join(['sample', 'genome', 'taxonomy', 'coverage', 'dispersion', 'fraction', 'containment', 'ptr']) + '\n')
            for row in sorted(self.data, key=lambda row: (row[0], row[2], row[1])):
                if row[-1] is not None:
                    f.write('\t'.join(row[:3] + [f'{x:.4f}' for x in row[3:]]) + '\n')
        log.info('Done.')

def profile(files, outdir, database, force=False, single=False, min_cove=5, max_disp=np.inf, min_frac=0.75, min_cont=0.25, components=5, max_iter=np.inf, tol=1e-5, threads=os.cpu_count()):
    '''
    Profile bacterial growth dynamics.
    '''
    gp = GrowthProfiler(files=files, outdir=outdir, database=database, force=force, single=single, threads=threads)
    gp.count(max_disp=max_disp, min_cove=min_cove, min_frac=min_frac, min_cont=min_cont)
    gp.infer(components=components, max_iter=max_iter, tol=tol)
    gp.write()
