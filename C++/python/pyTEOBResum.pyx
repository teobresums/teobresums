import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
cimport cython

cpdef np.ndarray[double, ndim=2, mode = 'c'] pyTEOBResum(double m1,
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
                                                         double LambdaBl2,
                                                         double distance,
                                                         int    NQC,
                                                         int    tidal,
                                                         int    speedy,
                                                         int    RWZ,
                                                         int    dynamics,
                                                         int    waveform,
                                                         int    lm,
                                                         int    solver_scheme):
    cdef Waveform *hp;
    cdef Waveform *hc;
    
    TEOBResum(&hp,
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
              LambdaBl2,
              distance,
              NQC,
              tidal,
              speedy,
              RWZ,
              dynamics,
              waveform,
              lm,
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
