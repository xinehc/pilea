import os
import glob
import shutil
import sys
import re

from collections import defaultdict
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from functools import partial
import packaging

from .utils import ctar
from .sketch import sketch
from .kmc import KMC
from .log import log
from . import __version__


def index(files, outdir, taxonomy=None, compress=False, database=None, k=31, s=250, w=25000, threads=os.cpu_count()):
    '''
    Update an existing database or construct a new one.
    '''
    pattern = '\\.(fa|fna|fasta)(\\.gz)?$'

    if len(files) == 1 and files[0].endswith('.txt'):
        if not os.path.isfile(files[0]):
            log.critical(f'File list <{files[0]}> does not exist.')
            sys.exit(2)
        else:
            with open(files[0]) as f:
                files = []
                for line in f:
                    files.append(line.rstrip())

    files = sorted([file for file in files if re.search(pattern, file) and '*' not in file], key=os.path.basename)
    names = [re.sub(pattern, '', os.path.basename(file)) for file in files]
    if not files:
        log.critical(f'No files with extension <*.(fa|fna|fasta).(gz)?> are given.')
        sys.exit(2)

    if database is not None:
        if not os.path.isdir(database):
            log.critical(f'Folder <{database}> does not exist.')
            sys.exit(2)

        for file in ['parameter.tab', 'taxonomy.tab', 'sketch.uni', 'sketch.dup', 'sketch.kmc_pre', 'sketch.kmc_suf']:
            if not os.path.isfile(f'{database}/{file}'):
                log.critical(f'File <{file}> does not exist under <{database}>.')
                sys.exit(2)

        with open(f'{database}/parameter.tab') as f:
            for line in f:
                if line[0] == 'v':
                    version = packaging.version.parse(line.split()[-1])
                    if packaging.version.parse(__version__) < version:
                        log.critical(f'Database <{database}> requires <v{version}> or above.')
                        sys.exit(2)
                else:
                    if int(line.split()[-1]) != locals()[line[0]]:
                        log.critical(f'Parameter <{line[0]}> does not match.')
                        sys.exit(2)

        with open(f'{database}/taxonomy.tab') as f:
            n = sum(1 for _ in f)
    else:
        n = 0

    tmpdir = f'{outdir}/tmp'
    if os.path.isdir(outdir):
        shutil.rmtree(tmpdir, ignore_errors=True)
    os.makedirs(tmpdir)

    if taxonomy is None:
        metadata = {name: 'n/a' for name in names}
    elif not os.path.isfile(taxonomy):
        log.critical(f'Taxonomy mapping file <{taxonomy}> does not exist.')
        sys.exit(2)
    else:
        with open(taxonomy) as f:
            metadata = dict(line.rstrip().split('\t')[:2] for line in f)

        ## special format for GTDB representatives
        if os.path.basename(taxonomy) == 'bac120_taxonomy.tsv':
            names = [name[:15] for name in names]
            metadata = {key[3:]: val for key, val in metadata.items()}

    log.info('Sketching ...')
    files = [(*file, n + idx) for idx, file in enumerate([(file, name) for file, name in zip(files, names) if name in metadata])]
    files = process_map(partial(sketch, folder=tmpdir, k=k, s=s, w=w), files, max_workers=threads, chunksize=1, leave=False)

    log.info('Preparing mapping files ...')
    with open(f'{outdir}/tmp.fa', 'wb') as f:
        for file in files:
            with open(f'{outdir}/tmp/{file[0]}.fa', 'rb') as g:
                shutil.copyfileobj(g, f, length=1024 * 1024)

    kmc = KMC(file=f'{outdir}/tmp.fa', prefix=f'{outdir}/tmp', threads=threads)
    kmc.count(k=k, f='fm', ci=1)
    kmc.compact(in_prefix=f'{outdir}/tmp', out_prefix=f'{outdir}/sketch')

    kmc = KMC(file=f'{outdir}/tmp.fa', prefix=f'{outdir}/tmp.sketch', threads=threads)
    kmc.count(k=k, f='fm', ci=2)
    kmc.dump(f'{outdir}/tmp.sketch')

    log.info('Writing sequences ...')
    with open(f'{outdir}/tmp.sketch') as f:
        kdup = {line[:k] for line in f}

    klst = defaultdict(list)
    with open(f'{outdir}/sketch.uni', 'w') as f, open(f'{outdir}/tmp.fa') as g:
        for line in g:
            if line[0] == '>':
                acc = line[1:].rstrip()
            else:
                if (seq := line.rstrip()) in kdup:
                    klst[seq].append(acc)
                else:
                    f.write(f'{seq}\t{acc}\n')

    with open(f'{outdir}/sketch.dup', 'w') as f:
        for key, val in klst.items():
            f.write(f'{key}\t{",".join(list(dict.fromkeys(val)))}\n')

    with open(f'{outdir}/taxonomy.tab', 'w') as f:
        for file in [file[1] for file in files]:
            f.write(f'{file[0]}\t{file[1]}\t{file[2]}\t{file[3]}\t{metadata[file[0]]}\n')

    with open(f'{outdir}/parameter.tab', 'w') as f:
        f.write('\n'.join((f'k\t{k}', f's\t{s}', f'w\t{w}', f'v\t{__version__}')) + '\n')

    ## merge databases if necessary
    if database is not None:
        log.info(f'Merging ...')
        for file in ['sketch.uni', 'sketch.dup', 'taxonomy.tab']:
            os.rename(f'{outdir}/{file}', f'{outdir}/tmp.{file}')

        kmc.intersect(in_prefix=f'{outdir}/sketch', ref_prefix=f'{database}/sketch', out_prefix=f'{outdir}/tmp.sketch')
        kmc.dump(prefix=f'{outdir}/tmp.sketch')
        with open(f'{outdir}/tmp.sketch') as f:
            kdup = {line[:k] for line in f}

        klst = defaultdict(list)
        for suffix in ['uni', 'dup']:
            with open(f'{outdir}/sketch.{suffix}', 'w') as f:
                for file in [f'{database}/sketch.{suffix}', f'{outdir}/tmp.sketch.{suffix}']:
                    with open(file) as g:
                        for line in g:
                            if line[:k] not in kdup:
                                f.write(line)
                            else:
                                klst[line[:k]].extend(line[k + 1:-1].split(','))

        with open(f'{outdir}/sketch.dup', 'a') as f:
            for key, val in klst.items():
                f.write(f'{key}\t{",".join(val)}\n')

        kmc.union(in_prefix=f'{outdir}/sketch', ref_prefix=f'{database}/sketch', out_prefix=f'{outdir}/sketch.union')
        for suffix in ['kmc_pre', 'kmc_suf']:
            shutil.move(f'{outdir}/sketch.union.{suffix}', f'{outdir}/sketch.{suffix}')

        with open(f'{outdir}/taxonomy.tab', 'wb') as f:
            for file in [f'{database}/taxonomy.tab', f'{outdir}/tmp.taxonomy.tab']:
                with open(file, 'rb') as g:
                    shutil.copyfileobj(g, f, length=1024 * 1024)

    if compress:
        log.info('Compressing ...')
        os.makedirs(f'{outdir}/database', exist_ok=True)
        for file in glob.glob(f'{outdir}/*.tab') + glob.glob(f'{outdir}/sketch.*'):
            shutil.copyfile(file, f'{outdir}/database/{os.path.basename(file)}')
        ctar(file=f'database', folder=outdir)
        shutil.rmtree(f'{outdir}/database')

    ## clean up
    for file in glob.glob(f'{outdir}/tmp*'):
        try:
            shutil.rmtree(file)
        except NotADirectoryError:
            os.remove(file)

    log.info('Done.')