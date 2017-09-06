from libcpp cimport bool

cdef extern from "TEOBResumS.h":
    ctypedef struct Waveform:
        double *data
        unsigned int length

    ctypedef struct TEOBResumFlags:
        int solver_scheme
        int spin
        int tidal
        int RWZ
        int speedy
        int dynamics
        int Yagi_fits
        int multipoles
        int geometric_units
        int set


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
                    double dt,
                    double LambdaAl2,
                    double LambdaBl2,
                    double LambdaAl3,
                    double LambdaBl3,
                    double LambdaAl4,
                    double LambdaBl4,
                    double distance,
                    int    lm,
                    TEOBResumFlags *flags
);
