import sys
import re
import os
import mmap
import glob
import struct
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

import packaging

from .log import log
from .ztp import ZTP
from .utils import u4, u8
from .kmc import count64
from .parse import parse_fastx_file

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
        if not extensions == extensions & {'fasta', 'fa', 'fastq', 'fq'}:
            log.critical('Input files need to have extensions in <fa|fq|fasta|fastq>.')
            sys.exit(2)

        if len(extensions) != 1:
            log.critical('Input files have mixed extensions.')
            sys.exit(2)

        if len(samples) > len(set(samples)):
            log.critical('Input files have duplicated names.')
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
                self.meta[i] = sp[0], sp[4], int(sp[3]), int(sp[2])

        if os.path.isdir(outdir):
            log.warning(f'Output directory <{outdir}> exists.')
        else:
            os.makedirs(outdir, exist_ok=True)

    @staticmethod
    def _collect(items, k, m, outdir, force, REC=struct.Struct('<QI')):
        sample, files = items
        msg = ''
        kmc = f'{outdir}/{sample}.kmc'
        if os.path.isfile(kmc) and not force:
            msg += f"File <{kmc}> exists, skip. Use <--force> for overwriting."
            with open(kmc, 'rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                counts = {key: cnt for (key, cnt) in REC.iter_unpack(mm)}
        else:
            counts = {}
            if len(files) == 1:
                it = parse_fastx_file(files[0])
                for record in it:
                    count64(record, None, k, m, counts)
            else:
                it1 = parse_fastx_file(files[0])
                it2 = parse_fastx_file(files[1])
                for record1, record2 in zip(it1, it2):
                    count64(record1, record2, k, m, counts)

                if next(it1, None) is not None or next(it2, None) is not None:
                    msg += f'Sample <{sample}> is not properly paired, skip. Make sure <{files[0]}> and <{files[1]}> have equal numbers of reads.'
                    counts.clear()

            if counts:
                with open(kmc, 'wb') as f:
                    for key, cnt in counts.items():
                        f.write(u8(key))
                        f.write(u4(cnt))
        return sample, counts, msg

    @staticmethod
    def _load(file, arr, REC=19, CHUNK_RECS=2500000):
        pdb = defaultdict(bytearray)
        if arr.size == 0:
            return pdb

        buf = bytearray(CHUNK_RECS * REC)
        total_recs = os.path.getsize(file) // REC
        total_chunks = (total_recs + CHUNK_RECS - 1) // CHUNK_RECS
        with open(file, 'rb') as f, tqdm(total=total_chunks, leave=False) as pbar:
            while True:
                nbytes = f.readinto(buf)
                if nbytes <= 0:
                    break

                n = nbytes // REC
                rows = np.frombuffer(buf, dtype=np.uint8, count=n * REC).reshape(n, REC)
                keys = rows[:, :8].view('<u8').ravel()

                order = np.argsort(keys)
                mask = np.isin(keys[order], arr, assume_unique=True)

                for i in order[np.flatnonzero(mask)]:
                    pdb[int(keys[i])].extend(rows[i, 8:19])

                pbar.update(1)
        return pdb

    @staticmethod
    def _assign(kmc, pdb, meta, min_cont):
        def _decode(b, gid, REC=struct.Struct('<IHHHB'), GID=struct.Struct('<I7x')):
            return GID.iter_unpack(b) if gid else REC.iter_unpack(b)

        ku, kd = defaultdict(list), defaultdict(set)
        for key, cnt in kmc.items():
            m = len(b := pdb[key]) // 11
            if m == 1:
                p = next(_decode(b, gid=False))
                ku[p[0]].append((p[1:4], cnt))
            else:
                for (gid,) in _decode(b, gid=True):
                    kd[gid].add(key)

        ka = {}
        kc = {gid: len(val) / meta[gid][-1] for gid, val in ku.items()}
        for gid in list(kd):
            if kc.get(gid, 0) + len(kd[gid]) / meta[gid][-1] <= min_cont:
                del kd[gid]

        while kd:
            ba = max(kd, key=lambda gid: (kc.get(gid, 0) + len(kd[gid]) / meta[gid][-1], gid))
            bs = kd.pop(ba)

            sa = set()
            for s in bs:
                sa.update(gid for (gid,) in _decode(pdb[s], gid=True))  # relevant gid

            for gid in sa & kd.keys():
                kd[gid] -= bs
                if kc.get(gid, 0) + len(kd[gid]) / meta[gid][-1] <= min_cont:
                    del kd[gid]

            kc[ba] = kc.get(ba, 0) + len(bs) / meta[ba][-1]
            ka[ba] = sorted(bs)

        for ba, bs in ka.items():
            for s in bs:
                i = (p for p in _decode(pdb[s], gid=False) if p[0] in ka)
                p = next(i)
                if next(i, None) is None and p[-1]:  # unique
                    ku[ba].append((p[1:4], kmc[s]))
        return ku

    @staticmethod
    def _filter(row, min_cove, max_disp, min_frac, min_cont):
        def _trim(x, return_limits=False):
            if x.size == 0:
                return x
            x = np.log2(x)
            q1, q3 = np.percentile(x, [25, 75])
            lower, upper = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)

            x = np.rint(2 ** x[(x >= lower) & (x <= upper)]).astype(np.uint32)
            return (2 ** lower, 2 ** upper) if return_limits else x

        def _debias(x, y, frac=0.25):
            t = lowess(exog=x, endog=y, frac=frac)[:, 1]
            y = np.rint(2 ** (np.asarray(y) - t + np.mean(t))).astype(np.uint32)
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
        arr = set()
        obs = dict()
        fun = partial(self._collect, k=self.k, m=self.m, outdir=self.outdir, force=self.force)
        with logging_redirect_tqdm():
            with mp.Pool(self.threads) as pool:
                for sample, counts, msg in tqdm(pool.imap_unordered(fun, self.items.items(), chunksize=1), total=len(self.items), leave=False):
                    obs[sample] = counts
                    arr.update(counts.keys())
                    if msg:
                        log.warning(msg)

        arr = np.fromiter(arr, dtype='<u8', count=len(arr))
        arr.sort()

        log.info('Loading database ...')
        pdb = self._load(f'{self.database}/sketches.pdb', arr)

        log.info('Parsing and filtering outputs ...')
        for sample in tqdm(list(obs.keys()), leave=False):
            kmc = {k: v for k, v in obs.pop(sample).items() if k in pdb}  # recast
            self.data.extend([[sample, *self.meta[gid][:-1], cont, p] for gid, p in self._assign(kmc, pdb, self.meta, min_cont).items() if (cont := len(p) / self.meta[gid][-1]) > min_cont])

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
