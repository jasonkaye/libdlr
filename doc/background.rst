
.. _Background:

Background
==========

For a brief background on the discrete Lehmann representation (DLR), we suggest reading the original reference on the DLR:

Discrete Lehmann representation of imaginary time Green's functions. Jason Kaye, Kun Chen, Olivier Parcollet. arXiv:2107.13094. 2021.

For more information on the Python module, pydlr, please refer to the quick-start section of this documentation.


Citation information
--------------------

If this library helps you to create software or publications, please let us know, and cite our repository

https://github.com/jasonkaye/libdlr

and the preprint referred to above.


Directory of examples
---------------------

For the Fortran library, ``libdlr``, the sample programs contained in the folder ``./libdlr/demo`` show standard examples of usage, and should serve as a useful starting point. The demos are thoroughly commented. For more information on specific subroutines, please refer to the `libdlr Fortran API`_.

The following is a list of the demos, with brief descriptions:

`sc_it.f90`_ - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR imaginary time grid.

`sc_it_fit.f90`_ - Recover DLR of Green's function with semi-circular spectral density from noisy samples on a uniform imaginary time grid.

`sc_mf.f90`_ - Recover DLR of Green's function with semi-circular spectral density from its samples on the DLR Matsubara frequency grid.

`syk_it.f90`_ - Solve the SYK equation (Dyson equation with SYK self-energy) by an imaginary time domain method.

`syk_mf.f90`_ - Solve the SYK equation (Dyson equation with SYK self-energy) by a more standard method, which computes the self-energy in the imaginary time domain, and solves the Dyson equation in the Matsubara frequency domain.

.. _sc_it.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_it.f90#L24
.. _sc_it_fit.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_it_fit.f90#L24
.. _sc_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/sc_mf.f90#L24
.. _syk_it.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/syk_it.f90#L24
.. _syk_mf.f90: https://github.com/jasonkaye/libdlr/blob/main/demo/syk_mf.f90#L24

Note on libdlr format for representing imaginary time points
------------------------------------------------------------

libdlr works with scaled imaginary time and real frequency coordinates,
in which an imaginary time tau is in [0,1] and a real frequency is in
[-Lambda,Lambda], for Lambda the dimensionless energy cutoff parameter.
Throughout the code, imaginary time grid points in (0.5,1) are stored in
a peculiar manner. Namely, suppose tau in (0.5,1). Then instead of
storing tau directly, we store the number tau* = tau-1.  In other words,
we store the negative distance of tau to 1, rather than tau itself. We
call this the relative format. To recover the correct tau points, one
simply leaves any points in [0,1/2] U {1} alone, and adds 1 to any
points in (1/2,1); the subroutine rel2abs performs this conversion.

The reason for this has to do with maintaining full relative accuracy in
floating point arithmetic. To evaluate the kernel K(tau,omega), we
sometimes need to compute the value 1-tau for tau very close to 1. If we
work with tau directly, there is a loss of accuracy due to catastrophic
cancellation, which begins to appear in extreme physical regimes and at
very high requested accuracies. If we instead compute tau* to full relative accuracy and
work with it directly rather than with tau -- for example by exploiting
symmetries of K(tau,omega) to avoid ever evaluating 1-tau -- we can
maintain full relative accuracy.

This apparently annoying characteristic of libdlr is simply the price of
maintaining all the accuracy that is possible to maintain in floating
point arithmic. It is the right thing to do. But it is also largely
ignoreable if the loss of accuracy is not noticeable in your application
(this will be the case for the vast majority of users). Just use the
provided demos to get started, and follow these guidelines:

(1) Use provided subroutines, which hide
this technical complexity, to carry out all operations.

(2) In a situation in which you want to provide a point tau
in (1/2,1) at which to evaluate a DLR, there are two options:
(i) compute tau* to full relative accuracy, and provide this according to
the instructions in the relevant subroutines, thereby maintaining full
relative accuracy in calculations, or (ii) if you don't care about the
(usually minor) digit loss which comes from ignoring this subtlety, simply convert the provided
tau point to tau* using the subroutine abs2rel, which carries out the
conversion. Since the point will have started its life in the absolute
format, converting it to relative format cannot recover full relative
accuracy, but it still needs to be converted in order to be compatible
with libdlr subroutines.

(3) If you happen to want to evaluate a Green's function on an
equispaced grid on [0,1] in imaginary time, use the subroutine eqpts_rel
to generate this grid in the relative format.
