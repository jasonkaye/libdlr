Background
==========

A brief background on the discrete Lehmann representation (DLR) will appear here soon.


For now, we suggesting reading the original reference on the DLR:

Discrete Lehmann representation of imaginary time Green's functions.
Jason Kaye, Kun Chen, Olivier Parcollet. arXiv:2107.13094. 2021.


For the Fortran library, libdlr, the sample programs contained in the folder libdlr/demo show standard examples of usage, and should serve as a useful starting point. The demos are thoroughly commented. For more information on specific subroutines, please refer to the libdlr (Fortran library) section of this documentation.

The following is a list of the demos, with brief descriptions:

sc_it.f90 - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR imaginary time grid.

sc_it_fit.f90 - Recover DLR of Green's function with semi-circular spectral density from noisy samples on a uniform imaginary time grid.

sc_mf.f90 - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR Matsubara frequency grid.

syk_it.f90 - Solve the SYK equation (Dyson equation with SYK self-energy) by an imaginary time domain method.

syk_mf.f90 - Solve the SYK equation (Dyson equation with SYK self-energy) by a more standard method, which computes the self-energy in the imaginary time domain, and solves the Dyson equation in the Matsubara frequency domain.



For more information on the Python module, pydlr, please refer to the quick-start section of this documentation.




If this library helps you to create software or publications, please let us know, and cite our repository

https://github.com/jasonkaye/libdlr

and the preprint referred to above.
