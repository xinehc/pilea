import subprocess
import os
import sys

from .log import log


def u1(x): return x.to_bytes(1, 'big', signed=False)
def u2(x): return x.to_bytes(2, 'big', signed=False)
def u3(x): return x.to_bytes(3, 'big', signed=False)
def u8(x): return x.to_bytes(8, 'big', signed=False)

def wget(file, folder):
    os.makedirs(folder, exist_ok=True)
    stderr = subprocess.run([
        'wget', '-S', '--spider', file,
    ], check=True, capture_output=True, text=True).stderr

    file_size = 0
    for line in stderr.splitlines():
        if 'content-length' in line.lower():
            file_size = int(line.split(':', 1)[1].strip())
            break

    file_path = f'{folder}/{file.split("/")[-1]}'
    if os.path.isfile(file_path) and file_size == os.path.getsize(file_path):
        log.info(f'File <{file_path}> exist, skip.')
    else:
        subprocess.run([
            'wget', '-q', '--show-progress', file, '-O', file_path
        ], check=True)
        sys.stdout.write("\033[F\033[K")
        sys.stdout.flush()

def xtar(file, folder, threads=os.cpu_count()):
    subprocess.run([
        'tar', '-I', f'pigz -p {threads} -dc', '-C', folder, '-xf', file,
    ], check=True)

def ctar(file, folder, threads=os.cpu_count()):
    subprocess.run([
        'tar', '-I', f'pigz -p {threads} --best', '-C', folder, '--sort=name', '-cf', f'{folder}/{file}.tar.gz', file
    ], check=True)
