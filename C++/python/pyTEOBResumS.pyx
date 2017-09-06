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
                                                         double sampling_rate,
                                                         double LambdaAl2,
                                                         double LambdaAl3,
                                                         double LambdaAl4,
                                                         double LambdaBl2,
                                                         double LambdaBl3,
                                                         double LambdaBl4,
                                                         double distance,
                                                         int    tidal,
                                                         int    speedy,
                                                         int    RWZ,
                                                         int    dynamics,
                                                         int    lm,
                                                         int    yagi,
                                                         int    solver_scheme):
    cdef Waveform *hp;
    cdef Waveform *hc;
    
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
              sampling_rate,
              LambdaAl2,
              LambdaAl3,
              LambdaAl4,
              LambdaBl2,
              LambdaBl3,
              LambdaBl4,
              distance,
              tidal,
              speedy,
              RWZ,
              dynamics,
              lm,
              yagi,
              solver_scheme)

    cdef int i = 0
    cdef unsigned int N = hp.length
    cdef np.ndarray[double, ndim=2] x = np.zeros((N,N), dtype=np.double)
    for i in range(N):
        x[i,0] = hp.data[i]
        x[i,1] = hc.data[i]
    free(hp)
    free(hc)
    return x
