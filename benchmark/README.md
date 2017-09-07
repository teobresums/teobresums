# Benchmark

## Tests

| ID | q | S1z | S2z | Lambda1 | Lambda2 | Notes |
|---|---|---|---|---|---|---|
| BBH00 | 1 |0 | 0 | 0 | 0 | r0=45, Geom.Units, 22 mode, parfile |
| BBH01 | 20 | 0 | 0 | 0 | 0 | r0=45, Geom.Units, 22 mode, parfile |
| BBH02 | 1 | 0.9 | 0.9 | 0 | 0 | r0=45, Geom.Units, 22 mode, parfile |
| BBH03 | 20 | 0.9 | 0.9 | 0 | 0 | r0=45, Geom.Units, 22 mode, parfile |
| BNS00 | 1 | 0 | 0 | 1531 | 1531 | r0=45, Geom.Units, 22 mode, parfile |
| BNS01 | 1.5 | 0 | 0 | 864 | 7021 | r0=45, Geom.Units, 22 mode, parfile |
| BNS02 | 1 | 0.4 | 0.4 | 1531 | 1531 | r0=45, Geom.Units, 22 mode, parfile |
| BNS03 | 1.5 | 0.4 | 0.4 | 864 | 7021 | r0=45, Geom.Units, 22 mode, parfile |

To add tests simply add parfiles (and update the table above) in the directory
```
$ ./parfiles/
```

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

The benchmark directory is named after the last git revision, a subdirectory with the date is created. Information about the git revision are stored in a text file; the whole source code is backupped.
