import subprocess
import os
import sys
import shutil
import numpy as np

from .log import log


def wget(file, folder):
    stderr = subprocess.run([
        'wget', '-S', '--spider', file,
    ], check=True, capture_output=True, text=True).stderr

    for line in stderr.splitlines():
        if 'content-length' in line.lower():
            file_zize = int(line.split(':', 1)[1].strip())
    
    file_path = f'{folder}/{file.split("/")[-1]}'
    if os.path.isfile(file_path) and file_zize == os.path.getsize(file_path):
        log.info(f'File <{file_path}> exist, skip.')
    else:
        subprocess.run([
            'wget', '-q', '--show-progress', file, '-O', file_path
        ], check=True)
        sys.stdout.write("\033[F\033[K")
        sys.stdout.flush()

def xtar(file, folder):
    subprocess.run([
        'tar', '-C', folder, '-zxf', file,
    ], check=True)

def ctar(file, folder):
    subprocess.run([
        'tar', '-C', folder, '--sort=name', '-zcf', f'{folder}/{file}.tar.gz', file
    ], check=True)
