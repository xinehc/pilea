# cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False

from libc.stdint cimport uint8_t


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

def parse_fastx_file(str filename):
    '''
    Yield sequence in bytes.
    '''
    cdef gzFile fp
    cdef kseq_t *ks
    cdef bytes b_filename = filename.encode('utf-8')
    
    fp = gzopen(b_filename, 'r')
    ks = kseq_init(fp)
    
    try:
        while kseq_read(ks) >= 0:
            yield ks.seq.s[:ks.seq.l]
    finally:
        if ks != NULL: kseq_destroy(ks)
        if fp != NULL: gzclose(fp)
