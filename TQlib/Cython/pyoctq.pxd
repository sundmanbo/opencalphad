# declarations for liboctqisoc interface
cdef extern from "octqd.h":
    cdef struct gtp_equilibrium_data:
        pass
       # int status
       # int multiuse
       # int eqno
       # int cnext
       # char *eqname
       # double tpval, rtn
       # double *svfunres
       # void *lastcondition 
       # void *lastexperiment
       # void *complist
       # void **compstoi
       # void **invcompstoi
       # void *phase_varres
       # void *eq_tpres
       # double *cmuval
       # double xconv
       # double gmindif
       # int maxiter
       # int sysmatdim 
       # int nfixmu
       # int nfixph
       # int *fixmu
       # int *fixph
       # double **savesysmat
    
    cdef void  c_tqini(int n, void **ceq)
    cdef void  c_tqrfil(char *, void **ceq)
    cdef void  c_tqgcom(int *, char[][24], void **ceq)
    cdef void  c_tqgnp(int *, void **ceq)
    cdef void  c_tqgpn(int, char *, void **ceq)
    cdef void  c_tqgetv(char *, int, int, int *, double *, void **ceq)
    cdef void  c_tqsetc(char *, int, int, double, int *, void **ceq)
    cdef void  c_tqce(char *, int, int, double *, void **ceq)
    #cdef int c_ntup
    #cdef int c_nel
    #cdef int c_maxc
    #cdef char *c_cnam[]



