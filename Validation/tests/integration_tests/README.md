# Integration tests

Tests to check the global behavior of `TEOBResumS`:
* `integration_tests.py` contains a number of sanity checks that can be automatically run via
    ```
    pytest -v integration_tests.py
    ```
    and either fail or pass.

* `parspace_and_plots.py` contains a number of useful routines and sanity checks that produce plots and require human intervention.
Between these are a number of routines useful to explore the parameter space of BBH, BNS and BHNS binaries and ensure that the waveform generator does not fail.
