import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
cimport cython

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

    cdef int i = 0
    cdef unsigned int N = hp.length
    cdef np.ndarray[double, ndim=2] x = np.zeros((N,N), dtype=np.double)
    for i in range(N):
        x[i,0] = hp.data[i]
        x[i,1] = hc.data[i]
    free(hp)
    free(hc)
    return x
