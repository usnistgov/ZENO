===========
Source code
===========

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

The source code, which can be found at https://github.com/usnistgov/ZENO,
is written in ``C++`` and requires a compiler that supports
the ``C++11`` standard; recent versions of ``g++`` have been found to
work.

Dependencies and compilation
----------------------------

The only essential external library is the SPRNG random number library
(which requires compilation). This can be obtained from:

-  http://sprng.org

When configuring the SPRNG library, you will have the option to include
support for Fortran and MPI.  Neither of these are required for ZENO (even
if you plan to use MPI with ZENO), so you may build SPRNG without these.
   
Other required libraries are zlib and a threading library.  These should
be provided by your Operating System distribution.

Compiling the code requires cmake version 3.1 or greater.

To compile, you should first make a build directory.  Then change to the
build directory and run ``cmake`` on the ``zeno/cpp`` directory.  You will
need to set the ``SPRNG_INCLUDE_DIR`` and ``SPRNG_LIBRARY`` variables based
on where you installed SPRNG.  For example:

::
   
   mkdir zeno/cpp/build
   cd zeno/cpp/build
   cmake -DSPRNG_INCLUDE_DIR=/opt/sprng5/include -DSPRNG_LIBRARY=/opt/sprng5/lib/libsprng.a ..

You should then be able to build the executable ``zeno`` by simply
typing ``make``.  Optionally, you may also install the executable into
the directory specified in cmake by typing ``make install``.

MPI support
-----------

MPI support is included, but is optional. If MPI libraries are installed
on your system, you may be able to build zeno with MPI support by running
cmake with

``-DUSE_MPI=ON``

Cmake will attempt to automatically find the MPI libraries on your system.
If this does not work, you will need to manually set the various cmake MPI
variables to the values appropriate for your installation.

Modifying the code
------------------

If you modify the source code, some external utilities may be required.
``Gengetopt`` is required to modify command-line parameters, while
``Bisonc++`` and ``Flexc++`` are required to modify the input file
format. These can be obtained from the sources:

-  https://www.gnu.org/software/gengetopt

-  https://fbb-git.github.io/bisoncpp

-  https://fbb-git.github.io/flexcpp

Simple check
------------

Once the code has been compiled, you can perform a quick self-test by
typing ``make check``. This will run ``zeno`` 
on various test cases and compare the output against the
output from a correctly built version. Floating-point values will be
allowed some tolerance to account for differences in compilers, machine
precision, etc.

.. _sec:runcode:

