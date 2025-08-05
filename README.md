# luscherzeta
A utility package to compute the Luscher zeta function.

## Installation
This package depends on the gsl library. References to gsl might need to be added to the Makefile under the `LIBS` and `INC` environment variables.

Download and build the repository:
```
git clone git@github.com:asmeyer2012/luscherzeta.git <your_directory>
cd zeta
make
```
This creates a shared object library `zeta/libczeta.so` that is referenced in the python wrapper `zeta/czeta.py`. This reference may need to be changed to an absolute path based on the installation path.

## Usage

There is a class implementation in python that wraps the C++ code. This is the file `zeta/czeta.py`, which has only a few handles. Here is the basic usage syntax in python:
```
from zeta.czeta import czeta
s_vector = [0.3, 0.1, 0] ## 3-vector of floats
boost_gamma = 1.1 ## float
l, m = 1, -1 ## angular momentum quanta (positive_integer, integer)
q2 = 0.5 ## dimensionless center of mass momentum

object_czeta = czeta()
object_czeta.set_svec_gamma( *s_vector, boost_gamma)
object_czeta.set_lm( l, m)
print( object_czeta.evaluate( q2))
```

