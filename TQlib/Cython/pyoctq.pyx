cimport pyoctq
from numpy cimport ndarray as arr
from numpy import array

cdef public gtp_equilibrium_data *ceq
cdef extern int c_ntup
cdef extern char *c_cnam[]
cdef extern int c_nel

def tqini(int n):
    #cdef gtp_equilibrium_data *ceq
    c_tqini(n, <void**> &ceq)

def tqrfil(filename):
    c_tqrfil(filename,<void**> &ceq)

def tqgpn(int i):
    cdef char *cphname=""
    c_tqgpn(i, cphname,<void**> &ceq)
    cdef bytes pyname = cphname
    return pyname

def tqsetc(cond, int n1, int n2, double tp):
    cdef int cnum
    c_tqsetc(<char *> cond, n1, n2,  tp, &cnum, <void **> &ceq)
    return cnum

def tqce(target, int n1, int n2, double value):
    c_tqce(<char*>target, n1, n2, &value, <void **> &ceq)
    return value

def tqgetv(statvar, int n1, int n2, int n3, arr[double] npf):
    c_tqgetv(<char *>statvar, n1, n2, &n3, &npf[0], <void **>&ceq)
    return n3

def ctuple():
    return c_ntup

def cname(int ph):
    cdef bytes phase = c_cnam[ph]
    return phase

def cnel():
    return c_nel;
    




