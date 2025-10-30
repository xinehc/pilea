import os
import glob

from .utils import wget, xtar
from .index import index
from .log import log


def rebuild(outdir, k=31, s=500, w=25000, threads=os.cpu_count()):
    '''
    Build a reference database by sketching GTDB assemblies.
    '''
    ## download and extract necessary files
    os.makedirs(outdir, exist_ok=True)
    log.info('Downloading assemblies ...')
    for file in ['bac120_taxonomy.tsv', 'genomic_files_reps/gtdb_genomes_reps.tar.gz']:
        wget(file=f'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/{file}', folder=outdir)

    log.info('Extracting ...')
    xtar(file=f'{outdir}/gtdb_genomes_reps.tar.gz', folder=outdir)

    ## sketch
    index(files=glob.glob(f'{outdir}/gtdb_genomes_reps*/**/*.fna.gz', recursive=True), outdir=outdir, taxonomy=f'{outdir}/bac120_taxonomy.tsv', compress=True, database=None, k=k, s=s, w=w, threads=threads)
