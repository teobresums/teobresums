To compile this code:

./bootstrap
./configure --prefix=<your destination> CXXFLAGS="-std=c++11 -Ofast -Wall -fPIC"
make
make install

to build and install the python wrapper

cd python 
python setup.py build_ext --inplace 

Once the installation is finished, add 
<your destination> to the path
To run
TEOBResum <options>
