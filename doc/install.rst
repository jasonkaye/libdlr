.. index:: installation

.. highlight:: bash

Installation
============

Requirements
------------

``libdlr`` (Fortran library)

- Fortran compiler
- BLAS / LAPACK
- CMake

  
``pydlr`` (Python module)

- Python 3.x
- Numpy
- Scipy
- ``nosetest`` (for running the tests)

Build
-----

Using ``CMake``::
  
   git clone https://github.com/jasonkaye/libdlr.git
   cd libdlr
   mkdir cbuild
   cd cbuild

Include python module (optional)::
  
   cmake -Dwith_python=ON ..

Build documentation (optional)::
  
   cmake -Dwith_python=ON -DBUILD_DOCS=ON ..

Compile and build::
  
   make

Run tests::

   make test

Install
-------

Set ``-DCMAKE_INSTALL_PREFIX="/some/path/to/install/to"`` when calling ``cmake`` and run::
  
  make install
