# Development LOG 

LOG file for TEOBResumS `development` branch

Last updated: 2018.03.07

## Main 

|File                           | Description               | Notes                                                         |
|-------------------------------|---------------------------|---------------------------------------------------------------|
| TEOBResumSUtils.c             | Utilities                 | 2017.11.07 Compiles                                           |
| TEOBResumSPars.c              | parameters and pars I/O   | 2017.11.07 Verified correctness (with standalone test)        |
| TEOBResumSFits.c              | Fits (final BH, NQC, etc) | 2017.11.07 Compiles                                           |
| TEOBResumSInitialConditions.c | Initial conditions        | 2018.02.13 Compiles<br>
                                                              2018.02.26 Implemented post-post-circular conditions          | 
| TEOBResumSMetric.c            | Metric functions          | 2018.02.09 Compiles, verified correctness (with C++ master)   |     
| TEOBResumSDynamics            | r.h.s. for dynamics       | 2018.02.12 Compiles                                           | 
| TEOBResumSFlux.c              | Flux                      | 2017.11.07 Compiles                                           | 
| TEOBResumSWaveform.c          | Waveform                  | 2017.12.02 Compiles                                           |
| TEOBResumS.c                  | Main                      | 2017.03.07 Compiles, made several changes also in other *.c   |
[Debug LOG. Please keep update with main development milestones]

## TODOLIST

 * Check runtime errors

 * Verify corretness, compatre with Matlab and C++ on master

 * Several parts needs checking, please do
 ```
 fgrep "//CHECKME" *
 fgrep "//FIXME" *
 ```

 * LR computation

 * LSO computation only for nonspining BBH

 * Stopping conditions need to be improved

 * Ringdown attachment needs to be improved








