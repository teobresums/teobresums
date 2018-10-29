import numpy as np
cimport numpy as cnp
cnp.import_array()
from libc.stdlib cimport malloc, free
cimport cython

cdef public api tonumpyarray(double* data, long long size) with gil:
    if not (data and size >= 0): raise ValueError
    cdef cnp.npy_intp dims = size
    #NOTE: it doesn't take ownership of `data`. You must free `data` yourself
    return cnp.PyArray_SimpleNewFromData(1, &dims, cnp.NPY_DOUBLE, <void*>data)


cpdef cnp.ndarray[double, ndim=2, mode = 'c'] pyTEOBResumS(double m1,
                                                         double m2,
                                                         double spin1x,
                                                         double spin1y,
                                                         double spin1z,
                                                         double spin2x,
                                                         double spin2y,
                                                         double spin2z,
                                                         double inclination,
                                                         double coa_phase,
                                                         double f_min,
                                                         double dt,
                                                         double LambdaAl2,
                                                         double LambdaBl2,
                                                         double LambdaAl3,
                                                         double LambdaBl3,
                                                         double LambdaAl4,
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
              coa_phase,
              f_min,
              dt,
              LambdaAl2,
              LambdaBl2,
              LambdaAl3,
              LambdaBl3,
              LambdaAl4,
              LambdaBl4,
              distance,
              lm,
              &flags)

    N = hp.length
    cdef cnp.ndarray[double, ndim=2, mode = 'c'] h = np.zeros((N,2),dtype=np.double)
    h[:,0] = tonumpyarray(hp.data,N)
    hview = hc.data
    h[:,1] = tonumpyarray(hc.data,N)
    if hp.data: free(hp.data)
    if hp : free(hp)
    if hc.data: free(hc.data)
    if hc: free(hc)
    h[np.isnan(h)] = 0.0
    return h

