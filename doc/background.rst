
.. _Background:

Background
==========

Imaginary time formalism
------------------------

The `imaginary time formalism <https://en.wikipedia.org/wiki/Imaginary_time>`_ is a tool used in many-body quantum theory to determine the equilibrium properties of systems at finite temperature. In this case the response functions of the system are functions of imaginary time :math:`\tau` that is defined on :math:`\tau \in [0, \beta]`, where :math:`\beta` is the inverse temperature.

In particular the `single-particle Green's function <https://en.wikipedia.org/wiki/Green%27s_function_(many-body_theory)>`_ :math:`G` of a system is a function in imaginary time, and can be written as the time ordered expectation value

.. math::

   G_{ab}(\tau) = - \langle \mathcal{T} c_a(\tau) c_b^\dagger(0) \rangle, 

where :math:`c^\dagger_b(0)` creates a particle in state :math:`b` at time :math:`0` and :math:`c_a(\tau)` destroys a particle in state :math:`a` at time :math:`\tau`. The Green's function -- and other derived quantities -- can be extended to the interval :math:`\tau \in [-\beta, \beta]` using the (anti-)periodicity property

.. math::
   
   G_{ab}(-\tau) = \xi G_{ab}(\beta - \tau),

where :math:`\xi = \pm 1` for bosonic and fermionic particles, respectively.


Discrete Lehmann Representation
-------------------------------

The discrete Lehmann representation (DLR) provides and efficient method of representing imaginary time functions numerically. The theory and methodology behind is based on the analysis of an integral kernel (`the analytic continuation kernel <https://en.wikipedia.org/wiki/Numerical_analytic_continuation>`_) combined with a clever linear algebra factorization (`the interpolative decomposition <https://en.wikipedia.org/wiki/Interpolative_decomposition>`_). 

For an in depth background on the discrete Lehmann representation (DLR), we suggest looking at the two papers cited just below.


.. _citations:

Citation information
--------------------

If this library helps you to create software or publications, please let
us know, and cite

- our repository: `<https://github.com/jasonkaye/libdlr>`_
- `libdlr: Efficient imaginary time calculations using the discrete Lehmann representation, Jason Kaye, Kun Chen, and Hugo U.R. Strand, arXiv:2110.06765. <https://arxiv.org/abs/2110.06765>`_
- `Discrete Lehmann representation of imaginary time Green's functions, Jason Kaye, Kun Chen, and Olivier Parcollet, Phys. Rev. B 105, 235115, 2022. <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.235115>`_ [`arXiv:2107.13094 <https://arxiv.org/abs/2107.13094>`_]


``libdlr`` imaginary time point format
--------------------------------------

``libdlr`` works with scaled imaginary time and real frequency coordinates,
in which an imaginary time :math:`\tau` is in :math:`\tau \in [0,1]` and a real frequency :math:`\omega` is in
:math:`\omega \in [-\Lambda,\Lambda]`, for :math:`\Lambda` the dimensionless energy cutoff parameter.

Throughout the code, imaginary time grid points in :math:`(0.5,1)` are stored in
a peculiar manner. Namely, suppose :math:`\tau \in (0.5,1)`. Then instead of
storing :math:`\tau` directly, we store the number :math:`\tau^* = \tau-1`.  In other words,
we store the negative distance of :math:`\tau` to 1, rather than tau itself. We
call this the relative format. To recover the correct tau points, one
simply leaves any points in :math:`[0,1/2] \cup {1}` alone, and adds 1 to any
points in :math:`(1/2,1)`; the subroutine ``rel2abs`` performs this conversion.

The reason for this has to do with maintaining full relative accuracy in
floating point arithmetic. To evaluate the kernel :math:`K(\tau,\omega)`, we
sometimes need to compute the value :math:`1-\tau` for :math:`\tau` very close to 1. If we
work with tau directly, there is a loss of accuracy due to catastrophic
cancellation, which begins to appear in extreme physical regimes and at
very high requested accuracies. If we instead compute :math:`\tau^*` to full relative accuracy and
work with it directly rather than with :math:`\tau`, for example by exploiting
symmetries of :math:`K(\tau,\omega)` to avoid ever evaluating :math:`1-\tau`, we can
maintain full relative accuracy.

This apparently annoying characteristic of ``libdlr`` is simply the price of
maintaining all the accuracy that is possible to maintain in floating
point arithmic. But it is largely
ignoreable if the loss of accuracy is not noticeable in your application
(this will be the case for the vast majority of users). Just use the
provided demos to get started, and follow these guidelines:

1. Use provided subroutines, which hide this technical complexity, to carry out all operations.

2. In a situation in which you want to provide a point :math:`\tau \in (1/2,1)`
   at which to evaluate a DLR, there are two options:

   - compute :math:`\tau^*` to full relative accuracy, and provide this according to
     the instructions in the relevant subroutines, thereby maintaining full
     relative accuracy in calculations, or
   - if you don't care about the
     (usually minor) digit loss which comes from ignoring this subtlety, simply convert the provided
     :math:`\tau`-point to :math:`\tau^*` using the subroutine ``abs2rel``, which carries out the
     conversion. Since the point will have started its life in the absolute
     format, converting it to relative format cannot recover full relative
     accuracy, but it still needs to be converted in order to be compatible
     with ``libdlr`` subroutines.

3. If you happen to want to evaluate a Green's function on an
   equispaced grid on :math:`[0,1]` in imaginary time, use the subroutine ``eqpts_rel``
   to generate this grid in the relative format.

A more in-depth discussion of the relative format is given in Appendix C
of `our paper on libdlr <https://arxiv.org/abs/2110.06765>`_.
