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
        for n, (_, fseq) in enumerate(SimpleFastaParser(f)):
            fseq = fseq.upper()
            kpos = [(canonicalize(fseq[i:i + k]), i + k // 2) for i in range(len(fseq) - k + 1)]
            kcnt.update([kmer for kmer, _ in kpos])
            ctgs[n] = (fseq, kpos)

    wins = 0
    skes = set()
    kdup = {key for key, val in kcnt.items() if val > 1}
    with open(f'{folder}/{idx}.fa', 'w') as f:
        for _, (fseq, kpos) in ctgs.items():
            for i in range(0, len(kpos), w):
                if ksub := [(kmer, pos, kmer in kdup) for kmer, pos in kpos[i:i + w] if hash64(kmer, signed=False)[0] < maxhash]:
                    sins = sum(not dup for _, _, dup in ksub)
                    for kmer, pos, dup in ksub:
                        skes.add(kmer)
                        gseq = fseq[max(0, pos - 500):pos + 500]
                        gcnt = int((gseq.count('G') + gseq.count('C')) / len(gseq) * 1000)
                        f.write(f'>{idx}|{wins}:{sins}|{gcnt}|{"-" if dup else "+"}\n{kmer}\n')
                    wins += 1
    return idx, (name, len(ctgs), wins, len(skes))
