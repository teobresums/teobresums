/*
 *  Copyright (C) 2017 Walter Del Pozzo, Gregorio Carullo, Sebastiano Bernuzzi, Alessandro Nagar, Ka Wa Tsang, Philipp Fleig
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include "LALSimIMRTEOB.h"

static int XLALSimIMRTEOB(REAL8TimeSeries **hplus,                     /**< [out] Time-domain waveform h+ */
                          REAL8TimeSeries **hcross,                    /**< [out] Time-domain waveform hx */
                          const REAL8 phi0,                            /**< Orbital phase at f_ref (rad) */
                          const REAL8 deltaT,                          /**< Sampling interval (s) */
                          const REAL8 m1_SI,                           /**< Mass of companion 1 (kg) */
                          const REAL8 m2_SI,                           /**< Mass of companion 2 (kg) */
                          const REAL8 f_min,                           /**< Starting GW frequency (Hz) */
                          const REAL8 f_ref,                           /**< Reference frequency (Hz) */
                          const REAL8 distance,                        /**< Distance of source (m) */
                          const REAL8 thetaJ,                          /**< Angle between J0 and line of sight (z-direction) [inclination] */
                          const REAL8 S1z,                             /**< Aligned-spin parameter of companion 1 */
                          const REAL8 S2z,                             /**< Aligned-spin parameter of companion 2 */
                          const REAL8 quadparam1,                      /**< quadrupole deformation parameter of body 1 (dimensionless, 1 for BH) */
                          const REAL8 quadparam2,                      /**< quadrupole deformation parameter of body 2 (dimensionless, 1 for BH) */
                          const REAL8 lambda1,                         /**< (tidal deformation of body 1)/(mass of body 1)^5 */
                          const REAL8 lambda2,                         /**< (tidal deformation of body 2)/(mass of body 2)^5 */
                          const LALSimInspiralTestGRParam *extraParams /**< linked list containing the extra testing GR parameters */
)
{
    /***************************************
     * call the input preprocessing function
     * call the setup function
     * call the driver/integrator function
     * call the NQC corrections
     * call the post-merger
     * set hplus and hcross and return
     ***************************************/
    return 0;
}
