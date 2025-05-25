import os
import gzip

from mmh3 import hash64
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

revcomp = bytes.maketrans(b'ACGT', b'TGCA')
alphabet = set('ACGT')


def canon(seq, trans=revcomp):
    return min(seq, seq[::-1].translate(trans))

def count(seq, pos, flank=500):
    seq = seq[max(0, pos - flank):pos + flank]
    return int((seq.count('G') + seq.count('C')) / len(seq) * flank * 2)

def sketch(file, folder, k, s, w):
    file, name, idx = file
    maxhash = ((1 << 64) - 1) // s

    ctgs = defaultdict(list)
    kcnt = Counter()
    with gzip.open(file, 'rt') if '.gz' in file else open(file, 'r') as f:
        for n, (_, fseq) in enumerate(SimpleFastaParser(f)):
            fseq = fseq.upper()
            kpos = [(i // w, (kmer, count(fseq, i + k // 2))) for i in range(len(fseq) - k + 1) if hash64((kmer := canon(fseq[i:i + k])), signed=False)[0] < maxhash and set(kmer) <= alphabet]
            kcnt.update(kmer for _, (kmer, _) in kpos)
            for i, j in kpos:
                ctgs[f'{n}|{i}'].append(j)

    wins = 0
    kdup = {key for key, val in kcnt.items() if val > 1}
    with open(f'{folder}/{idx}.fa', 'w') as f:
        for widx, (_, kpos) in enumerate(ctgs.items()):
            if (sins := sum(not kmer in kdup for kmer, _ in kpos)) > 1:
                wins += 1
            for kmer, gc in kpos:
                f.write(f'>{idx}|{widx}:{sins}|{gc}|{"-" if kmer in kdup else "+"}\n{kmer}\n')
    return idx, (name, len(set(i.split('|')[0] for i in ctgs.keys())), wins, len(kcnt))
