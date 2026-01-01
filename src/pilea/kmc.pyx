# cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3

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
    cdef uint64_t f = 0, r = 0, canon, h, mask
    cdef int valid = 0, top_shift
    cdef Py_ssize_t i
    cdef uint8_t code

    top_shift = 2 * (k - 1)
    mask = ((<uint64_t>1) << (2 * k)) - 1

    for i in range(n):
        code = BASE_MAP[s[i]]
        if code == 0xFF:
            valid = 0; f = 0; r = 0
            continue

        valid += 1
        f = ((f << 2) | code) & mask
        r = _rc(r, code, top_shift)
        if valid >= k:
            canon = f if f <= r else r
            h = _hash(canon)
            if h < maxhash:
                if OutputContainer is list:
                    out.append((i - k + 1, h))
                else:
                    out.add(h)

cpdef void count64(bytes a, bytes b, int k, uint64_t maxhash, dict counts) except *:
    '''
    Count unique FracMinHash k-mers from single or paired reads.
    '''
    cdef Py_ssize_t n
    cdef char* buffer
    cdef set tmp = set()

    PyBytes_AsStringAndSize(a, &buffer, &n)
    _scan(<const uint8_t*>buffer, n, k, maxhash, tmp)

    if b is not None:
        PyBytes_AsStringAndSize(b, &buffer, &n)
        _scan(<const uint8_t*>buffer, n, k, maxhash, tmp)

    cdef uint64_t key
    for key in tmp:
        counts[key] = counts.get(key, 0) + 1

cpdef list hash64(bytes seq, int k, uint64_t maxhash):
    '''
    Store FracMinHash k-mers and their positions.
    '''
    cdef Py_ssize_t n
    cdef char* buffer
    cdef list tmp = []

    PyBytes_AsStringAndSize(seq, &buffer, &n)
    _scan(<const uint8_t*>buffer, n, k, maxhash, tmp)

    return tmp