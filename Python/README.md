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
python setup.py install
```
to install in your entire environment.

## Others

Other contents of this folder are

 * `CWrap/` : Tests for the python wrapper
 * `Examples/` : Examples to run `TEOBResumS` through the python interface 
 

