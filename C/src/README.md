# Development LOG 

LOG file for TEOBResumS `development` branch

Last updated: 2018.03.07

## Brief introduction for developers

### Parameters

 * Parameters are amnaged with specific and minimal routine based on libconfig.
 * Pars are stored in a data structure and can be read/written in all the place of the code.
 * For safety, change and write pars only in the main. Do otherwise only if strictly necessary.
 * Some pars are stored in some specific data structure (see below) for quick access. Be careful and do not abuse of this.  
 * The code is initialized with a parameter file only.
 * Default pars are listed in `par/default.par`. Do not remove or overwrite this file; other example files can be found in the same dir.

### Data structures

Special and minimal data structures are defined together with routines to allocate, free and push the storage.

 * `Dynamic` Variables and storage for the dynamics, including some quick-access pars.
 * `Waveform` Waveform data type for h+, hx.
 * `Waveform_lm` Multipolar waveform data type.
 * `Waveform_lm_t` Multipolar waveform data type at given time, comes at handy for given applications.
 * If you need to start to understand what is in the code, the reference files should be `TEOBResumS.h`

### 


## Development timeline

|File                           | Description               | Notes                                                         |
|-------------------------------|---------------------------|---------------------------------------------------------------|
| TEOBResumSUtils.c             | Utilities                 | 2017.11.07 Compiles                                           |
| TEOBResumSPars.c              | parameters and pars I/O   | 2017.11.07 Verified correctness (with standalone test)        |
| TEOBResumSFits.c              | Fits (final BH, NQC, etc) | 2017.11.07 Compiles                                           |
| TEOBResumSInitialConditions.c | Initial conditions        | 2018.02.13 Compiles                                           |
|                                                          || 2018.02.26 Implemented post-post-circular conditions          | 
| TEOBResumSMetric.c            | Metric functions          | 2018.02.09 Compiles, verified correctness (with C++ master)   |     
| TEOBResumSDynamics            | r.h.s. for dynamics       | 2018.02.12 Compiles                                           | 
| TEOBResumSFlux.c              | Flux                      | 2017.11.07 Compiles                                           | 
| TEOBResumSWaveform.c          | Waveform                  | 2017.12.02 Compiles                                           |
| TEOBResumS.c                  | Main                      | 2017.03.07 Compiles, made several changes also in other *.c   |
[Debug LOG. Please keep update with main development milestones]

## TODOLIST

 * Check runtime errors

 * Verify corretness, compatre with Matlab and C++ on master

 * Several parts needs checking; to see what please do
 ```
 fgrep "//CHECKME" *
 fgrep "//FIXME" *
 fgrep "//TODO" *
 ```

 * LR computation

 * LSO computation only for nonspining BBH

 * Stopping conditions need to be improved

 * Ringdown attachment needs to be improved

 * Improve code documentation, especially in `TEOBResumS.h`







