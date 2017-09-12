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
  /** Input parameters */
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

  /** Alloc data structures */
  int size = par_get_i("size"); /* note: size can vary */
  Waveform *hpp; /* (h+,hx) */
  Waveform *hlm; /* h_lm */
  Waveform_alloc (&hpp, size, "hpp");
  Waveform_lm_alloc (&hlm, size, "hlm");

  /** Call main routine */
  int err = TEOBResumS(&hpp, &hlm);
  
  /** Handle errors */
  // todo: list of errors...
  if (err!=OK) errorexit("Main function returned error");

  /** Output */
  Waveform_lm_output (hpp);
  Waveform_lm_output (hlm);

  /** Free memory */
  Waveform_free (hpp);
  Waveform_lm_free (hlm);
  
  return OK;
}
