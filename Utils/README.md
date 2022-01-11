# Python utilies

Various scripts to help running the `TEOBResumS` executable with the parfiles.

## Usage

* Write all example parfiles: 
    ```
    python3 parfile.py
    ```

* Generate, run and plot parfiles based on GW150914.par:
    ```
    python3 multipleruns.py --fname GW150914.par --path ./ --gen_pars --run --multiprocess --nproc 2 --plot_wf
    ```

* Determine NQC iteratively with n=5 iterations: 
    ```
    python3 iterateNQC.py -i GW150914.par -n 5
    ```
