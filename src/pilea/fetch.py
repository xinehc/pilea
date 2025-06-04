import os
import sys
import subprocess

from .utils import wget, xtar
from .log import log


def fetch(outdir):
    '''
    Donwload pre-built database from Zenodo.
    '''
    record = 15469467
    stdout = subprocess.run([
        'wget', '-qO-', f'https://doi.org/10.5281/zenodo.{record}',
    ], check=True, capture_output=True, text=True).stdout
    file = [x for x in stdout.splitlines() if 'database.tar.gz' in x][0].split('"')[-2]

    log.info('Downloading ...')
    wget(file=file, folder=outdir)

    log.info('Extracting files ...')
    xtar(file=f'{outdir}/database.tar.gz', folder=outdir)
    os.remove(f'{outdir}/database.tar.gz')

    log.info('Done.')
