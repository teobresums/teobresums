import numpy as np
cimport numpy as np
np.import_array()
from libc.stdlib cimport malloc, free
cimport cython

cdef public api tonumpyarray(double* data, long long size) with gil:
    if not (data and size >= 0): raise ValueError
    cdef np.npy_intp dims = size
    #NOTE: it doesn't take ownership of `data`. You must free `data` yourself
    return np.PyArray_SimpleNewFromData(1, &dims, np.NPY_DOUBLE, <void*>data)


cpdef np.ndarray[double, ndim=2, mode = 'c'] pyTEOBResumS(double m1,
                                                         double m2,
                                                         double spin1x,
                                                         double spin1y,
                                                         double spin1z,
                                                         double spin2x,
                                                         double spin2y,
                                                         double spin2z,
                                                         double inclination,
                                                         double polarisation,
                                                         double f_min,
                                                         double dt,
                                                         double LambdaAl2,
                                                         double LambdaAl3,
                                                         double LambdaAl4,
                                                         double LambdaBl2,
                                                         double LambdaBl3,
                                                         double LambdaBl4,
                                                         double distance,
                                                         int    lm,
                                                         dict   waveflags):
    cdef Waveform *hp;
    cdef Waveform *hc;
    cdef TEOBResumFlags flags;

    flags.solver_scheme = waveflags['solver_scheme']
    flags.spin = waveflags['spin']
    flags.tidal = waveflags['tidal']
    flags.RWZ = waveflags['RWZ']
    flags.speedy = waveflags['speedy']
    flags.dynamics = waveflags['dynamics']
    flags.Yagi_fits = waveflags['Yagi_fits']
    flags.multipoles = waveflags['multipoles']
    flags.geometric_units = waveflags['geometric_units']
    flags.set = waveflags['set']

    cdef int i = 0
    cdef unsigned int N
    cdef np.ndarray[double, ndim=2] x

    TEOBResumS(&hp,
              &hc,
              m1,
              m2,
              spin1x,
              spin1y,
              spin1z,
              spin2x,
              spin2y,
              spin2z,
              inclination,
              polarisation,
              f_min,
              dt,
              LambdaAl2,
              LambdaAl3,
              LambdaAl4,
              LambdaBl2,
              LambdaBl3,
              LambdaBl4,
              distance,
              lm,
              &flags)

    N = hp.length
    cdef np.ndarray[double, ndim=1, mode = 'c'] hplus  = tonumpyarray(hp.data,N)
    cdef np.ndarray[double, ndim=1, mode = 'c'] hcross = tonumpyarray(hc.data,N)
    if hp.data: free(hp.data)
    if hp : free(hp)
    if hc.data: free(hc.data)
    if hc: free(hc)
    return np.column_stack((hplus, hcross))

