To compile this code:

$ make

To run this code open the par_file.par and choose your own setup, then:

$ ./EOB_ihes.out par_file.par output_dir

How to use the param_file:

NQC                            // Turn on NQC corrections. BBH only.
tidal                          // Turn on tidal effects.   BNS only.
RWZ                            // Switch to Regge-Wheeler-Zerilli normalization. It just changes the normalization of the waveform by a constant factor.
speedy                         // Turn on a faster computation of the tail effects, TO BE CHECKED.
dynamics                       // Write the dynamics in an output file.
waveform                       // Write the waveform in an output file.
lm                             // PROBABLY UNUSED CHECK 
dt                             // Sets the uniform time step for the evolution.
solver_scheme                  // Select the solver scheme. 0 is fully adaptive, 1 is adaptive till the LSO, 2 is fixed time step



