
.. _Fortran examples:

Fortran examples
================

Sample programs for the Fortran library ``libdlr`` are contained in the folders ``./libdlr/demo`` and ``./libdlr/test``. Both show standard examples of usage, and can serve as a useful starting point, however we recommend looking at the ``demo`` folder first. If you do not find the example you need there, then look through the programs in the ``test`` folder. These also serve as compilation test programs and are typically designed to test specific operations using analytically solvable examples. Both the demos and tests are thoroughly commented. For more information on specific subroutines, please refer to the :ref:`API Documentation ``libdlr```. Below, we provide a list of all demo and test programs, with brief descriptions.

We note, in particular, the C test program `ha_it.c`_, which provides an
example of the C interface included with ``libdlr``.

The Python examples section of this documentation provides a more
step-by-step tutorial of the use of the DLR in standard situations than
is given here, but the Fortran demos and tests should be sufficient for
users who have read at least one of the papers listed in the
:ref:`Background` section.

Demo programs
~~~~~~~~~~~~~

DLR from Matsubara frequency
----------------------------

`sc_mf.f90`_ - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR Matsubara frequency grid.

DLR from imaginary time 
-----------------------

`sc_it.f90`_ - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR imaginary time grid.

Fit noisy imaginary time data
-----------------------------

`sc_it_fit.f90`_ - Recover DLR of Green's function with semi-circular spectral density from noisy samples on a uniform imaginary time grid.

SYK model in imaginary time
---------------------------

`syk_it.f90`_ - Solve the SYK model (Dyson equation with SYK self-energy) by an imaginary time domain method.

SYK model in Matsubara frequency
--------------------------------

`syk_mf.f90`_ - Solve the SYK model (Dyson equation with SYK self-energy) by a more standard method, which computes the self-energy in the imaginary time domain, and solves the Dyson equation in the Matsubara frequency domain.

From time to frequency and back
-------------------------------

`timetofreqtotime.f90`_ - For a Green's function with a delta function spectral
density, form the DLR expansion from samples in imaginary time, evaluate at the
Matsubara frequency DLR nodes, form another DLR expansion from these values, and
evaluate back in imaginary time, making sure the original Green's function is recovered. 

Test programs
~~~~~~~~~~~~~

DLR from Matsubara frequency
----------------------------

`ha_mf.f90`_ - Recover DLR of Green's function with sum-of-delta-functions spectral density from its samples on the DLR Matsubara frequency grid.

DLR from imaginary time 
-----------------------

`ha_it.f90`_ - Recover DLR of Green's function with sum-of-delta-functions spectral density from its samples on the DLR imaginary time grid.

C example: DLR from imaginary time
----------------------------------

`ha_it.c`_ - Same as ha_it.f90, but in C. This program demonstrates the
use of the C interface included with ``libdlr``.

Convolution in imaginary time
-----------------------------

`conv_exp.f90`_ - Convolve two Green's functions represented by their
values on the DLR imaginary time grid, obtaining the result on the same
grid.

Dyson equation with fixed self-energy
-------------------------------------

`dyson_sc.f90`_ - Solve the Dyson equation with a fixed self-energy both
in imaginary time and in Matsubara frequency.

Matrix-valued Dyson equation with fixed self-energy
---------------------------------------------------

`dyson_mat.f90`_ - Solve matrix-valued Dyson equation with a fixed self-energy in
imaginary time.

L^2 inner products
------------------

`ip_exp.f90`_ - Compute the L^2 inner product of two Green's functions
in imaginary time.


.. _sc_it.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_it.f90#L24
.. _sc_it_fit.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_it_fit.f90#L24
.. _sc_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_mf.f90#L24
.. _syk_it.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/syk_it.f90#L24
.. _syk_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/syk_mf.f90#L24
.. _ha_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/test/ha_mf.f90#L24
.. _ha_it.f90: https://github.com/jasonkaye/libdlr/blob/main/test/ha_it.f90#L24
.. _ha_it.c: https://github.com/jasonkaye/libdlr/blob/main/test/ha_it.c#L23
.. _conv_exp.f90: https://github.com/jasonkaye/libdlr/blob/main/test/conv_exp.f90#L24
.. _dyson_sc.f90: https://github.com/jasonkaye/libdlr/blob/main/test/dyson_sc.f90#L24
.. _dyson_mat.f90: https://github.com/jasonkaye/libdlr/blob/main/test/dyson_mat.f90#L24
.. _ip_exp.f90: https://github.com/jasonkaye/libdlr/blob/main/test/ip_exp.f90#L24
