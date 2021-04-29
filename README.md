
# The Discrete Lehmann Representation Library **libdlr**

This library provides an implementation of the discrete
Lehmann representation for functions depending on imaginary time,
in particular single-particle Green's functions.

Author: Jason Kaye (2021)

Contributors: Hugo UR Strand (2021)

## Usage

While the library is implemented in Fortran it also provides a stable C API
that can be called from C/C++ as well as a Python wrapper using Numpy/Scipy.

For example programs in Fortran see the `./test/` directory.

## Build instructions

Requirements
- BLAS and LAPACK
- CMake > 3.12

Optional requirements for the Python module `pydlr`:
- Python, Numpy, CFFI

To enable the python module and tests pass the additional flag `-Dwith_python=ON` to cmake below.

To build the library and the test programs run

```
mkdir cbuild
cd cbuild
FC=icc cmake ..
make
```

This will use the Intel Fortran compiler `icc` to build the shared library `libdlr` and the test programs in the folder `./cbuild`. (If the Fortran compiler is not specified CMake will try to use any available compliant compiler it can find.)


## Example programs

To run the Fortran test programs, go to the test folder. Each .f90 test program has some
parameters that can be adjusted at the top.

For the time being, look through the test codes to see how
subroutines are used, along with the documentation of subroutines in the
source files in the src folder for explanations of variables and
subroutines.

## Imaginary time grid

Throughout the code, imaginary time (tau) grid points that are in
(0.5,1) are stored in a peculiar manner. Namely, suppose tau in (0.5,1).
Then instead of storing t, we store the number tau' = tau-1. In other words,
we store the negative distance of tau to 1, rather than tau itself. To
recover the correct tau points (for example from the imaginary time DLR
grid dlrit), simply leave any points in [0,1/2] U {1} alone, and add 1 to any
points in (1/2,1).

The reason for this has to do with maintaining full relative accuracy in
floating point arithmetic. To evaluate the kernel K(tau,omega), we
sometimes need to compute the value 1-tau. In particular, we sometimes
need to compute 1-tau for tau very close to 1. If we work with tau
directly, there is an immediate loss of accuracy due to catastrophic
cancellation. If we instead compute tau' to full relative accuracy and
work with it directly rather than with tau, for example by exploiting
symmetries of K(tau,omega) to avoid evaluating 1-tau, we maintain full
relative accuracy. Note, in particular, that tau' cannot be computed by
first computing tau and then subtracting 1; this would defeat the
purpose and lose the digits the use of tau' is intended to maintain.

To avoid confusion, users should therefore use provided subroutines,
which hide this technical complexity, to carry out all operations. In a
situation, for example, in which the user wants to provide a point tau
in (0,1/2) at which to evaluate a DLR, there are two options: (i)
compute tau' to full relative accuracy, and provide this according to
the instructions in the relevant subroutines, thereby maintaining full
relative accuracy in calculations, or (ii) simply convert the provided
tau point to tau' by subtracting 1, or preferably, process all input
points through the subroutine xxxxxx, which will carry out the
conversion for the user. (ii) is the simpler option, but may result in a
mild loss of accuracy in calculations.
