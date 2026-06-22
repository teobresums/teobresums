# Python interface


## Installation
### PyPI
To install the python TEOBResumS extension run

```
pip install teobresums
```

Note: to install from Test PyPI, use 
```
pip install --extra-index-url https://testpypi.python.org/pypi teobresums
```
instead

### From source

If you want to build the extension from source, do

```
python setup.py build_ext --inplace
```
or type `make` for a local build, and
```
pip install .
```
to install in your entire environment.

Note that macros can be enabled during compilations via environment variables.
For instance, to run the code in debug or verbose mode, do:
```
TEOB_DEBUG=1 TEOB_VERBOSE=1 pip install .
```

To set user-defined functions, instead, do:
```
TEOB_USER_FUNCTIONS=1 pip install .
```

Additional options can be inspected directly from the setup.py

## Others

Other contents of this folder are

 * `CWrap/` : Tests for the python wrapper
 * `Examples/` : Examples to run `TEOBResumS` through the python interface 
 

