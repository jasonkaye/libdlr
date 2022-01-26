.. index:: installation

.. highlight:: bash

Installation
============

Python module ``pydlr``
-----------------------

The Python module ``pydlr`` is available on the `Python Package Index (PyPi) <https://pypi.org/project/pydlr/>`_ and can be installed by running::

  pip3 install pydlr

If you do not have root privileges add the flag ``--user`` to install to your home directory. Note that ``pydlr`` requires Python version 3 as well as the standard Python modules ``numpy`` and ``scipy``.

Developer installation
^^^^^^^^^^^^^^^^^^^^^^

If you want to contribute to the development of  ``pydlr`` consider doing a "developer installation" by cloning the repo and set your ``PYTHONPATH`` environment variable manually::

   cd /my/install/directory
   git clone https://github.com/jasonkaye/libdlr.git
   export PYTHONPATH=/my/install/directory/libdlr:$PYTHONPATH
   python -c "import pydlr; print(pydlr.__file__)"

The last command is a test that the ``pydlr`` module can be found by python and should print::

  /my/install/directory/libdlr/__init__.py

where ``/my/install/directory`` is your system-specific install location.

For a permanent setup, consider adding the ``export ...`` command above to the initialization script of your shell, e.g. to ``~/.bashrc``, ``~/.profile``, or similar.

The ``pydlr`` unit-tests can be run using::

  cd /my/install/directory/libdlr/pydlr/test
  python -m unittest

Fortran library ``libdlr``
--------------------------

The Fortran library ``libdlr`` has the following requirements/dependencies

- Fortran compiler
- BLAS / LAPACK
- CMake >=v3.18

Build
^^^^^

``libdlr`` is using the `CMake build system <https://cmake.org/>`_.

To build ``libdlr`` from source run::
  
   git clone https://github.com/jasonkaye/libdlr.git
   cd libdlr
   mkdir cbuild
   cd cbuild
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make
   make test
   make install

Build documentation (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To build the documentation locally add the CMake flags::
  
   cmake -DCMAKE_BUILD_TYPE=Release -Dwith_python=ON -DBUILD_DOCS=ON ..

Custom install location
^^^^^^^^^^^^^^^^^^^^^^^
   
To install ``libdlr`` to a non standard location, e.g. in your home directory, add the CMake flag::

  -DCMAKE_INSTALL_PREFIX="/my/install/directory"

when calling ``cmake``.

Custom Blas/Lapack
^^^^^^^^^^^^^^^^^^

CMake will attempt to auto-detect the available `BLAS libraries <https://cmake.org/cmake/help/latest/module/FindBLAS.html>`_ and `LAPACK libraries <https://cmake.org/cmake/help/latest/module/FindLAPACK.html>`_.

If the auto-detection fails and you have libraries installed in custom locations try adding your custom location in the CMake flag::

  -DCMAKE_LIBRARY_PATH="/custom/path/to/blas/lapack/"

Linking with ``libdlr``
^^^^^^^^^^^^^^^^^^^^^^^

With ``libdlr`` correctly installed on your system, linking with the library requires only the link flag ``-ldlr``.

When using a custom install location, the path also needs to be specified::

  -L/my/install/directory/lib -ldlr

The same path has also to be present in your ``LD_LIBRARY_PATH`` environment variable at runtime.
  
Support
^^^^^^^

If you experience issues installing ``libdlr`` or ``pydlr`` please consider openinig an issue on `our GitHub repository <https://github.com/jasonkaye/libdlr/issues>`_.
