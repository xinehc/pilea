import os
import glob
import shutil
import sys
import re

from tqdm.contrib.concurrent import process_map
from functools import partial
import packaging

from .utils import ctar
from .sketch import sketch
from .log import log
from . import __version__


def index(files, outdir, taxonomy=None, compress=False, database=None, k=31, s=500, w=25000, threads=os.cpu_count()):
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

        for file in ['parameters.tab', 'genomes.tab', 'sketches.pdb']:
            if not os.path.isfile(f'{database}/{file}'):
                log.critical(f'File <{file}> does not exist under <{database}>.')
                sys.exit(2)

        with open(f'{database}/parameters.tab') as f:
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

        with open(f'{database}/genomes.tab') as f:
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
    files = [(*file, n + gid) for gid, file in enumerate([(file, name) for file, name in zip(files, names) if name in metadata])]
    files = process_map(partial(sketch, folder=tmpdir, k=k, s=s, w=w), files, max_workers=threads, chunksize=1, leave=False)

    log.info('Writing ...')
    with open(f'{outdir}/parameters.tab', 'w') as f:
        f.write(f'k\t{k}\ns\t{s}\nw\t{w}\nv\t{__version__}\n')

    with open(f'{outdir}/genomes.tab', 'w') as f:
        for file in [file[1] for file in files]:
            f.write(f'{file[0]}\t{file[1]}\t{file[2]}\t{file[3]}\t{metadata[file[0]]}\n')

    with open(f'{outdir}/sketches.pdb', 'wb') as f:
        for file in files:
            with open(f'{outdir}/tmp/{file[0]}.pdb', 'rb') as g:
                shutil.copyfileobj(g, f, length=1024 * 1024)

    ## merge databases if necessary
    if database is not None:
        log.info(f'Merging ...')
        for file in ['genomes.tab', 'sketches.pdb']:
            os.rename(f'{outdir}/{file}', f'{outdir}/tmp.{file}')
            with open(f'{outdir}/{file}', 'wb') as f:
                for file in [f'{database}/{file}', f'{outdir}/tmp.{file}']:
                    with open(file, 'rb') as g:
                        shutil.copyfileobj(g, f, length=1024 * 1024)

    if compress:
        log.info('Compressing ...')
        os.makedirs(f'{outdir}/database', exist_ok=True)
        for file in ['parameters.tab', 'genomes.tab', 'sketches.pdb']:
            shutil.copyfile(f'{outdir}/{file}', f'{outdir}/database/{file}')
        ctar(file=f'database', folder=outdir, threads=threads)
        shutil.rmtree(f'{outdir}/database')

    ## clean up
    shutil.rmtree(f'{outdir}/tmp', ignore_errors=True)
    for file in glob.glob(f'{outdir}/tmp*'):
        os.remove(file)

    log.info('Done.')
