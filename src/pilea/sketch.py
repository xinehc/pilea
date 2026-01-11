import os
import numpy as np
from collections import Counter, defaultdict

from .kmc import hash64
from .parse import parse_fastx_file
from .utils import u1, u2, u4, u8

GC_LUT = np.zeros(256, dtype=np.uint32)
GC_LUT[ord(b'G')] = 1; GC_LUT[ord(b'g')] = 1
GC_LUT[ord(b'C')] = 1; GC_LUT[ord(b'c')] = 1


def scan(center, gc_prefix, flank=250):
    n = len(gc_prefix) - 1
    l = 0 if center - flank < 0 else center - flank
    r = n if center + flank > n else center + flank
    return round((gc_prefix[r] - gc_prefix[l]) / (r - l) * 1e4)

def sketch(file, folder, k, s, w):
    file, name, gid = file
    maxhash = ((1 << 64) - 1) // s

    kctg = defaultdict(list)
    kcnt = Counter()
    for cid, seq in enumerate(parse_fastx_file(file)):
        s = np.frombuffer(seq, dtype=np.uint8)
        gc_prefix = np.r_[np.uint32(0), np.cumsum(GC_LUT[s], dtype=np.uint32)]
        for i, key in hash64(seq, k, maxhash):
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
