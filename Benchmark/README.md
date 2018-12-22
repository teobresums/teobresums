# Benchmark

## Tests

Parfile for several relevant test cases are collected in the dir

```
$ ../C/par/
```

To add a test simply add parfiles there. A basic but complete set of tests comprises

 * equal/unequal masses (q=1/q~=1)
 * tides/no tides (bns/bbh)
 * spins/no spins
 * equal/unequal spins
 * post-adiabatic/full ODE evolution
 
## How to run

Type

```
$ make build # prepare directory and compile exe
$ make run # run the test

```

or simply

```
$ make # build and run
```

Notes:

 * The benchmark directory is named after the last git revision, a subdirectory with the date is created.
 * Information about the git revision are stored in a text file; the whole source code is backupped. 
