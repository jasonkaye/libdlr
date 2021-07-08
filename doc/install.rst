.. index:: installation

.. highlight:: bash

Installation
============

Requirements
------------

``dlrcode`` (Fortran library)

- Fortran compiler
- Blas / Lapack
- ``id_dist`` library (bundled with ``dlrcode``)
- CMake

  
``pydlr`` (Python module)

- Python 3.x
- Numpy
- Scipy
- ``nosetest`` (for running the tests)

Build
-----

Using ``CMake``::
  
   git clone https://github.com/jasonkaye/dlrcode.git
   cd dlrcode
   mkdir cbuild
   cd cbuild

Only static Fortran library::
  
   cmake -DBUILD_SHARED_LIBS=OFF ..

Include python module::
  
   cmake -Dwith_python=ON ..

Build documentation::
  
   cmake -Dwith_python=ON -DBUILD_DOCS=ON ..

Compile and build::
  
   make

Install
-------

Set ``-DCMAKE_INSTALL_PREFIX="/some/path/to/install/to"`` when calling ``cmake`` and run::
  
  make install
