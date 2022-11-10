
``libdlr``: Imaginary time calculations using the Discrete Lehmann Representation
=================================================================================

The ``libdlr`` library provides Fortran routines for efficient imaginary time calculations of single particle Green's functions using the **Discrete Lehmann Representation** (DLR), as well as a stand-alone Python module ``pydlr`` implementing the same functionality. For more information on the DLR and the imaginary time formalism, please see the :ref:`Background`.

Take a look at our recent preprint for an introduction to the library:

- `libdlr: Efficient imaginary time calculations using the discrete Lehmann representation, Jason Kaye, Kun Chen, and Hugo U.R. Strand, Comput. Phys. Commun. 280, 108458, 2022. <https://www.sciencedirect.com/science/article/pii/S0010465522001771>`_ [`arXiv:2110.06765 <https://arxiv.org/abs/2110.06765>`_]
  
We note that ``libdlr`` provides an easy to use C interface. An example
program is described in the :ref:`Fortran examples` section. There is
also a Julia package, `Lehmann.jl <https://github.com/numericaleft/Lehmann.jl>`_, which implements similar
functionality to ``libdlr``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   background.rst
   fortran_examples.rst
   python_examples.ipynb
   libdlr.rst
   pydlr.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
