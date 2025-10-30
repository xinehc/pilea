# cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3

from libc.stdint cimport uint64_t, uint8_t
from cpython.unicode cimport PyUnicode_AsUTF8AndSize
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING


cdef uint8_t BASE_MAP[256]
for i in range(256): BASE_MAP[i] = 0xFF
BASE_MAP[ord('A')] = 0; BASE_MAP[ord('a')] = 0
BASE_MAP[ord('C')] = 1; BASE_MAP[ord('c')] = 1
BASE_MAP[ord('G')] = 2; BASE_MAP[ord('g')] = 2
BASE_MAP[ord('T')] = 3; BASE_MAP[ord('t')] = 3

cdef inline uint64_t _rc(uint64_t r, uint64_t code, int top_shift) nogil:
    return (r >> 2) | ((code ^ 0x3) << top_shift)

## https://github.com/lh3/minimap2/blob/master/sketch.c
cdef inline uint64_t _hash(uint64_t key) nogil:
    key = (~key) + (key << 21)
    key ^= key >> 24
    key += (key << 3) + (key << 8)
    key ^= key >> 14
    key += (key << 2) + (key << 4)
    key ^= key >> 28
    key += (key << 31)
    return key

cdef void _scan(const uint8_t* s, Py_ssize_t n, int k, uint64_t maxhash, set out) except *:
    cdef uint64_t f = 0
    cdef uint64_t r = 0
    cdef Py_ssize_t i
    cdef uint8_t code
    cdef int valid = 0
    cdef int top_shift = 2*(k-1)

    cdef uint64_t mask = ((<uint64_t>1) << (2*k)) - 1
    cdef uint64_t canon, h

    ## prime
    for i in range(k):
        code = BASE_MAP[s[i]]
        if code == 0xFF:
            valid = 0; f = 0; r = 0
        else:
            valid += 1
            f = ((f << 2) | code) & mask
            r = _rc(r, code, top_shift)

    if valid == k:
        canon = f if f <= r else r
        h = _hash(canon)
        if h < maxhash:
            out.add(h)

    ## slide
    for i in range(k, n):
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
                out.add(h)

cpdef void count64(str a, str b, int k, uint64_t maxhash, dict counts) except *:
    '''
    Count unique k-mers from single or paired reads.
    '''
    cdef Py_ssize_t n
    cdef const uint8_t* p
    cdef set tmp = set()

    p = <const uint8_t*> PyUnicode_AsUTF8AndSize(a, &n)
    _scan(p, n, k, maxhash, tmp)

    if b != '':
        p = <const uint8_t*> PyUnicode_AsUTF8AndSize(b, &n)
        _scan(p, n, k, maxhash, tmp)

    cdef uint64_t key
    for key in tmp:
        counts[key] = counts.get(key, 0) + 1

cpdef uint64_t hash64(str s):
    '''
    Return 64-bit hash of an A/C/G/T string; 0xFFFFFFFFFFFFFFFF if invalid.
    '''
    cdef Py_ssize_t n
    cdef const uint8_t* p = <const uint8_t*> PyUnicode_AsUTF8AndSize(s, &n)

    cdef int k = <int>n
    cdef int top_shift = 2*(k-1)
    cdef uint64_t mask = ((<uint64_t>1) << (2*k)) - 1

    cdef uint64_t f = 0
    cdef uint64_t r = 0
    cdef Py_ssize_t i
    cdef uint8_t code

    for i in range(n):
        code = BASE_MAP[p[i]]
        if code == 0xFF:
            return <uint64_t>0xFFFFFFFFFFFFFFFF
        f = ((f << 2) | code) & mask
        r = _rc(r, code, top_shift)

    cdef uint64_t canon = f if f <= r else r
    return _hash(canon)
