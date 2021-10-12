
Fortran examples
================

For the Fortran library, ``libdlr``, the sample programs contained in the folder ``./libdlr/demo`` show standard examples of usage, and should serve as a useful starting point. The demos are thoroughly commented. For more information on specific subroutines, please refer to the :ref:`API Documentation ``libdlr```.

The following is a list of the demos, with brief descriptions:

DLR from Matsubara frequency
----------------------------

`sc_mf.f90`_ - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR Matsubara frequency grid.

DLR from imaginary time 
-----------------------

`sc_it.f90`_ - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR imaginary time grid.

Fit noisy imaginary time data
-----------------------------

`sc_it_fit.f90`_ - Recover DLR of Green's function with semi-circular spectral density from noisy samples on a uniform imaginary time grid.

Solve SYK in imaginary time
---------------------------

`syk_it.f90`_ - Solve the SYK model (Dyson equation with SYK self-energy) by an imaginary time domain method.

Solve SYK in Matsubara frequency and imaginary time
---------------------------------------------------

`syk_mf.f90`_ - Solve the SYK model (Dyson equation with SYK self-energy) by a more standard method, which computes the self-energy in the imaginary time domain, and solves the Dyson equation in the Matsubara frequency domain.

.. _sc_it.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_it.f90#L24
.. _sc_it_fit.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_it_fit.f90#L24
.. _sc_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_mf.f90#L24
.. _syk_it.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/syk_it.f90#L24
.. _syk_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/syk_mf.f90#L24
