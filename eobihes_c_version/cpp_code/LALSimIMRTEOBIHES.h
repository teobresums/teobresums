#ifndef _LALSimIMRTEOBIHES_h
#define _LALSimIMRTEOBIHES_h

#include <vector>
using namespace::std;
typedef struct tagWaveform
{
    double *data;
    unsigned int length;
}Waveform;

void XLALSimIMRTEOBIHES(Waveform **hplus,          /** h+ return array **/
                        Waveform **hcross,         /** hx return array **/
                        double m1,              /** m1(Msun) **/
                        double m2,              /** m2(Msun) **/
                        double spin1x,          /** dimensionless s1x **/
                        double spin1y,          /** dimensionless s1y **/
                        double spin1z,          /** dimensionless s1z **/
                        double spin2x,          /** dimensionless s2x **/
                        double spin2y,          /** dimensionless s2y **/
                        double spin2z,          /** dimensionless s2z **/
                        double inclination,     /** inclination angle (rad) **/
                        double polarisation,    /** polarisation angle (rad) **/
                        double f_min,           /** starting frequency(Hz) **/
                        double sampling_rate,   /** sampling rate(Hz) **/
                        double LambdaAl2,       /** (tidal deformation of body 1)/(mass of body 1)^5 **/
                        double LambdaBl2,       /** (tidal deformation of body 2)/(mass of body 2)^5 **/
                        double distance,        /** distance(m) **/
                        bool   NQC,             /** NQC corrections flag (BBH only) **/
                        bool   tidal,           /** tidal corrections flag (BNS only) **/
                        bool   speedy,          /** accelerated tails flag **/
                        bool   RWZ,             /** Regge-Wheeler-Zerilli potential (?) **/
                        int    lm,              /** TO BE REMOVED **/
                        int    solver_scheme    /** integration scheme (0:adaptive,1:fixed step) **/
                        );

#endif /* _LALSimIMRTEOBIHES_h */
