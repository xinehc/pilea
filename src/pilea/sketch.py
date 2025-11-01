import os
import numpy as np

from .kmc import hash64
from .utils import u1, u2, u4, u8

from collections import Counter, defaultdict
from needletail import parse_fastx_file


def scan(center, gc_prefix, flank=500):
    n = len(gc_prefix) - 1
    l = 0 if center - flank < 0 else center - flank
    r = n if center + flank > n else center + flank
    return round((gc_prefix[r] - gc_prefix[l]) / (r - l) * 1e4)

def sketch(file, folder, k, s, w):
    file, name, gid = file
    maxhash = ((1 << 64) - 1) // s

    kctg = defaultdict(list)
    kcnt = Counter()
    for cid, record in enumerate(parse_fastx_file(file)):
        seq = record.seq
        s = np.frombuffer(seq.encode('ascii'), np.uint8)
        gc_prefix = np.r_[np.int32(0), np.cumsum((s == 71) | (s == 103) | (s == 67) | (s == 99), dtype=np.int32)]
        for i in range(len(seq) - k + 1):
            if (key := hash64(seq[i:i + k])) < maxhash:
                kcnt[key] += 1
                kctg[(cid, i // w)].append((key, scan(i + k // 2, gc_prefix)))

    win = 0
    kdup = {key for key, val in kcnt.items() if val > 1}
    with open(f'{folder}/{gid}.pdb', 'wb') as f:
        for wid, (_, kpos) in enumerate(kctg.items()):
            if (sin := sum(not key in kdup for key, _ in kpos)) > 1:
                win += 1

            for key, gc in kpos:
                f.write(u8(key))
                f.write(u4(gid))
                f.write(u2(wid))
                f.write(u2(sin))
                f.write(u2(gc))
                f.write(u1(0 if key in kdup else 1))

    return gid, (name, cid + 1, len(kcnt), win)
