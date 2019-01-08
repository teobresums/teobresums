# Benchmark

## Usage 

 * Before starting to change the code, run th ebenchmark and store the results
 * If you are done with your work, or if you are stuck, run it again
 * Compare the results

In most cases you expect to have not broken/changed the code when
updating say, from revision 0001 to revision 0002, during your xmas
holiday project. So to be sure, just diff the benchmark result: 

```
$ diff bench_0001/181225/bbh_q1_s0s0_10Hz/waveform.txt diff bench_0002/181225/bbh_q1_s0s0_10Hz/waveform.txt 
```

In other cases you do expect a difference because you fixed the
code. So you will need to inspect the data and see if you have done a
good work.

PLEASE USE THIS BENCHMARK, it is easy and can safe you/us lot of time.

## What to run?

Parfiles for several relevant test cases are collected in the dir

```
$ ../C/par/
```

Those parfile can be added in the list `FILES` at the beginning of the makefile in this directory.  
 
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
 * Information about the git revision are stored in a text file; the whole source code is backupped for safety.


