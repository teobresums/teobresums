from libcpp cimport bool

cdef extern from "TEOBResumS.h":
    ctypedef struct Waveform:
        double *data
        unsigned int length

    void TEOBResumS(Waveform **hplus,
                    Waveform **hcross,
                    double m1,
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
                    int    solver_scheme);
