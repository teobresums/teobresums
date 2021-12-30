# TEOBResumS C code

## Compilation & Usage

See the main [README](../README.md) or the [Wiki](https://bitbucket.org/eob_ihes/teobresums/wiki/Home)

## Brief introduction for developers

### Parameters

In v3.0 GIOTTO

 * Simplified handling of EOBParameters and interface
 * Many 0/1 options have now become "yes"/"no". Similarly, in place of enums we use the respective string. 
   For example, `use_lambda234_fits = 1` has become `use_lambda234_fits = "YAGI13"`

In v2.0

 * Parameters are managed with the struct type `EOBParameters`
 * For backwards compatibility and to work with parfiles we maintained the pars database using `libconfig` (will be abandoned in future versions)
 * No need of `par/default.par`

In v1.0

 * Parameters are managed with specific and minimal routines based on `libconfig`
 * Pars are stored in a data structure and can be read/written in all the places of the code using the proper routines.
 * For safety, change and write pars only in the main. Do otherwise only if strictly necessary.
 * Some pars are stored in some specific data structure (see below) for quick access. Be careful and do not abuse of this.  
 * The code is initialized with a parameter file. Should be kept the only option.
 * Default pars are listed in `par/default.par`. Do not remove or overwrite this file; other example files can be found in the same dir.

### Data structures

Special and minimal data structures are defined together with routines to allocate, free and push the storage.

 * `Dynamics` Variables and storage for the dynamics.
 * `SpinDynamics` Variables and storage for the dynamics, including some quick-access pars.
 * `Waveform` Waveform data type for h+, hx.
 * `Waveform_lm` Multipolar waveform data type.
 * `Waveform_lm_t` Multipolar waveform data type at given time, comes at handy for given applications.

and similar types for frequency-domain waveforms.

### How to start coding?

If you need to start with the code, the reference file is `TEOBResumS.h`











