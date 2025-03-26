import os
import gzip

from mmh3 import hash64
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser

revcomp = bytes.maketrans(b"ACTG", b"TGAC")


def canonicalize(seq, trans=revcomp):
    return min(seq, seq[::-1].translate(trans))

def sketch(file, folder, k, s, w):
    file, name, idx = file
    maxhash = int((2 ** 64 - 1) / s)

    ctgs = dict()
    kcnt = Counter()
    with gzip.open(file, 'rt') if '.gz' in file else open(file, 'r') as f:
        for n, (_, seq) in enumerate(SimpleFastaParser(f)):
            klst = [canonicalize(seq[i:i + k]) for i in range(len(seq) - k + 1)]
            kcnt.update(klst)
            ctgs[n] = klst

    wids = 0
    seqs = set()
    kdup = {key for key, val in kcnt.items() if val > 1}
    with open(f'{folder}/{idx}.fa', 'w') as f:
        for n, klst in ctgs.items():
            for i, j in enumerate(range(0, len(klst), w)):
                if ksub := [x for x in klst[j:j + w] if hash64(x, signed=False)[0] < maxhash]:
                    wids += 1
                    seqs.update(ksub)
                    sign = ['-' if x in kdup else '+' for x in ksub]
                    f.write('\n'.join(f'>{idx}|{n}|{i}|{sign.count('+')}|{y}\n' + x for x, y in zip(ksub, sign)) + '\n')
    return idx, (name, len(ctgs), wids, len(seqs))
