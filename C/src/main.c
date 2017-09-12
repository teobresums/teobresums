/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
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


#include "TEOBResumS.h"

int main (int argc, char* argv[])
{
  /* Input parameters */
  if (argc < 1) {
    par_commandline_parse("--help", 1);
    exit(OK);
  }
  if (argc == 1) {
    TEOBResumSSetParameters(argv, argc, INPUT_FILE, PR);
  }
  if (argc > 1) {
    TEOBResumSSetParameters(argv, argc, COMMAND_LINE, PR);
  }
   
  /* Alloc data structures */
  // todo: use a single structure with several Waveforms (h+,hx) and multipoles
  Waveform *hplus; 
  Waveform *hcross;
        
  /* Call main routine */
  int err = TEOBResumS(&hplus,
		       &hcross,
		       m1,
		       m2,
		       0.0,
		       0.0,
		       chi1,
		       0.0,
		       0.0,
		       chi2,
		       inclination,
		       polarisation,
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
		       &flags);
  /* Handle errors */
  // todo: list of errors...
  if (err!=OK) errorexit("Main function returned error");

  /* Output */
  // todo: proper output routines:
  // - get the data structure
  // - output each waveform that has been flagged
  /*
  std::FILE* f = std::fopen(output, "w");
  int i        = 0;
  int N        = hplus->length;
  for (i=0;i<N;i++)
    {
      std::fprintf(f, "%f\t%e\t%e\n", i*dt, hplus->data[i], hcross->data[i]);
    }
  std::fclose(f);
  free(hplus->data);
  free(hcross->data);
  free(hplus);
  free(hcross);
  */

}
    


    return OK;
}
