import subprocess
import os
import shutil


class KMC:
    def __init__(self, file, prefix, threads=os.cpu_count()):
        """
        https://github.com/refresh-bio/KMC
        """
        self.file = file
        self.prefix = prefix
        self.threads = min(threads, 128)

        os.makedirs(f'{self.prefix}.kmc', exist_ok=True)
        if isinstance(self.file, list):
            with open(f'{self.prefix}.kmc/id', 'w') as f:
                f.write('\n'.join(self.file) + '\n')
            self.file = f'@{self.prefix}.kmc/id'

    def count(self, k=31, f='fm', ci=1):
        subprocess.run([
            'kmc', f'-k{k}', f'-{f}', f'-ci{ci}', '-cs2147483647', f'-t{self.threads}', self.file, self.prefix, f'{self.prefix}.kmc'
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.rmtree(f'{self.prefix}.kmc', ignore_errors=True)

    def compact(self, in_prefix, out_prefix):
        subprocess.run([
            'kmc_tools', f'-t{self.threads}', 'transform', in_prefix, 'compact', out_prefix
        ], check=True, stderr=subprocess.DEVNULL)

    def dump(self, prefix):
        subprocess.run([
            'kmc_tools', f'-t{self.threads}', 'transform', prefix, 'dump', f'{prefix}'
        ], check=True, stderr=subprocess.DEVNULL)

    def union(self, in_prefix, ref_prefix, out_prefix):
        subprocess.run([
            'kmc_tools', f'-t{self.threads}', 'simple', in_prefix, ref_prefix, 'union', out_prefix, '-ocmin'
        ], check=True, stderr=subprocess.DEVNULL)

    def intersect(self, in_prefix, ref_prefix, out_prefix):
        subprocess.run([
            'kmc_tools', f'-t{self.threads}', 'simple', in_prefix, ref_prefix, 'intersect', out_prefix, '-ocmax'
        ], check=True, stderr=subprocess.DEVNULL)
