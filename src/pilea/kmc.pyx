# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False

from libc.stdint cimport uint64_t, uint8_t
from cpython.bytes cimport PyBytes_AsStringAndSize


cdef uint8_t BASE_MAP[256]
for i in range(256): BASE_MAP[i] = 0xFF
BASE_MAP[ord('A')] = 0; BASE_MAP[ord('a')] = 0
BASE_MAP[ord('C')] = 1; BASE_MAP[ord('c')] = 1
BASE_MAP[ord('G')] = 2; BASE_MAP[ord('g')] = 2
BASE_MAP[ord('T')] = 3; BASE_MAP[ord('t')] = 3

cdef inline uint64_t _rc(uint64_t r, uint64_t code, int top_shift) nogil:
    return (r >> 2) | ((code ^ 0x3) << top_shift)

cdef inline uint64_t _hash(uint64_t key) nogil:
    key = (~key) + (key << 21)
    key ^= key >> 24
    key += (key << 3) + (key << 8)
    key ^= key >> 14
    key += (key << 2) + (key << 4)
    key ^= key >> 28
    key += (key << 31)
    return key

ctypedef fused OutputContainer:
    set
    list

cdef void _scan(const uint8_t* s, Py_ssize_t n, int k, uint64_t maxhash, OutputContainer out) except *:
    cdef:
        uint64_t f = 0, r = 0, c, h, mask = (1ULL << (2 * k)) - 1
        int valid = 0, top_shift = 2 * (k - 1)
        Py_ssize_t i
        uint8_t code

    for i in range(n):
        code = BASE_MAP[s[i]]
        if code == 0xFF:
            valid = 0; f = 0; r = 0
            continue

        valid += 1
        f = ((f << 2) | code) & mask
        r = _rc(r, code, top_shift)
        if valid >= k:
            c = f if f <= r else r
            h = _hash(c)
            if h < maxhash:
                if OutputContainer is list:
                    out.append((i - k + 1, h))
                else:
                    out.add(h)

def count64(bytes seq_a, bytes seq_b, int k, uint64_t maxhash, dict counts):
    '''
    Count unique FracMinHash k-mers from single or paired reads.
    '''
    cdef:
        Py_ssize_t n
        char* b
        set tmp = set()
        uint64_t key

    PyBytes_AsStringAndSize(seq_a, &b, &n)
    _scan(<const uint8_t*>b, n, k, maxhash, tmp)

    if seq_b is not None:
        PyBytes_AsStringAndSize(seq_b, &b, &n)
        _scan(<const uint8_t*>b, n, k, maxhash, tmp)

    for key in tmp:
        if key in counts:
            counts[key] += 1
        else:
            counts[key] = 1

def hash64(bytes seq, int k, uint64_t maxhash):
    '''
    Store FracMinHash k-mers and their positions.
    '''
    cdef:
        Py_ssize_t n
        char* b
        list tmp = list()

    PyBytes_AsStringAndSize(seq, &b, &n)
    _scan(<const uint8_t*>b, n, k, maxhash, tmp)

    return tmp
