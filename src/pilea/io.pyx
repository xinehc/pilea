# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False

from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t
from libc.stdio cimport fopen, fread, fclose, FILE
from libc.stdlib cimport malloc, free


cdef extern from *:
    '''
    #include <zlib.h>
    #include <htslib/kseq.h>

    #ifndef KSEQ_INIT_GZ
    #define KSEQ_INIT_GZ
    KSEQ_INIT(gzFile, gzread)
    #endif
    '''
    ctypedef void *gzFile
    gzFile gzopen(const char *path, const char *mode)
    int gzclose(gzFile file)

    ctypedef struct kstring_t:
        size_t l
        char *s

    ctypedef struct kseq_t:
        kstring_t seq

    kseq_t *kseq_init(gzFile fd)
    void kseq_destroy(kseq_t *ks)
    int kseq_read(kseq_t *ks)

cdef extern from *:
    """
    #include <htslib/khash.h>

    KHASH_SET_INIT_INT64(kobs)
    static inline khash_t(kobs)* kobs_init(void) {
        return kh_init(kobs);
    }
    static inline void kobs_destroy(khash_t(kobs)* h) {
        kh_destroy(kobs, h);
    }
    static inline void kobs_add(khash_t(kobs)* h, uint64_t key) {
        int ret; kh_put(kobs, h, key, &ret);
    }
    static inline int kobs_has(khash_t(kobs)* h, uint64_t key) {
        return kh_get(kobs, h, key) != kh_end(h);
    }
    """
    ctypedef struct kobs_t "khash_t(kobs)":
        pass
    kobs_t* kobs_init()
    void kobs_destroy(kobs_t*)
    void kobs_add(kobs_t*, uint64_t)
    int  kobs_has(kobs_t*, uint64_t)

def parse_fastx_file(str filename):
    '''
    Yield sequence from fasta/fastq (gzip optional) in bytes.
    '''
    cdef:
        gzFile f
        kseq_t* k

    f = gzopen(filename.encode('utf-8'), 'r')
    k = kseq_init(f)

    try:
        while kseq_read(k) >= 0:
            yield k.seq.s[:k.seq.l]
    finally:
        if k != NULL: kseq_destroy(k)
        if f != NULL: gzclose(f)

def load_kmc(str filename):
    '''
    Load kmer counts into dict of int.
    '''
    cdef:
        FILE* f
        uint8_t* b
        size_t i, n, l, offset
        uint64_t key
        uint32_t val
        dict counts = {}
        int chunk = 10000000 * 12

    b = <uint8_t*> malloc(chunk)
    f = fopen(filename.encode('utf-8'), 'rb')
    try:
        while True:
            l = fread(b, 1, chunk, f)
            if l == 0: break

            n = l // 12
            for i in range(n):
                offset = i * 12
                key = (<uint64_t*> &b[offset])[0]
                val = (<uint32_t*> &b[offset+8])[0]
                counts[key] = val
        return counts
    finally:
        if f != NULL: fclose(f)
        if b != NULL: free(b)

def load_pdb(str filename, dict obs):
    cdef:
        FILE* f
        uint8_t* b
        size_t i, n, l, offset
        uint64_t key, observed_key
        kobs_t* k
        dict counts, pdb = {}
        int chunk = 10000000 * 19

    k = kobs_init()
    b = <uint8_t*> malloc(chunk)
    f = fopen(filename.encode('utf-8'), 'rb')
    try:
        for counts in obs.values():
            for observed_key in counts:
                kobs_add(k, <uint64_t>observed_key)

        while True:
            l = fread(b, 1, chunk, f)
            if l == 0: break

            n = l // 19
            for i in range(n):
                offset = i * 19
                key = (<uint64_t*> &b[offset])[0]
                if kobs_has(k, key):
                    if key in pdb:
                        pdb[key].extend(b[offset + 8 : offset + 19])
                    else:
                        pdb[key] = bytearray(b[offset + 8 : offset + 19])
        return pdb
    finally:
        if k != NULL: kobs_destroy(k)
        if f != NULL: fclose(f)
        if b != NULL: free(b)
