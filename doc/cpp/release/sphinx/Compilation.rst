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

The only essential external libraries are the nanoflann Nearest
Neighbor library (header-only, does not require compilation) and the
SPRNG random number library (requires compilation). These can be
obtained from the sources:

-  https://github.com/jlblancoc/nanoflann

-  http://sprng.org

You will need to set ``NANOFLANN_DIR`` and ``SPRNG_DIR`` at the top of
the ``Makefile`` to point to the locations of these libraries,
respectively.

You should then be able to build the executable ``zeno`` by simply
typing ``make``.

MPI support
-----------

MPI support is included, but is optional. If MPI libraries are installed
on your system, you may be able to build the MPI-enabled executable
``zeno-mpi`` by simply typing ``make mpi``. If this does not work, you
will need to change ``MPI_CXX``, ``MPI_CXXFLAGS``, and ``MPI_LDFLAGS``
in the ``Makefile`` to values appropriate for your installation.

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
typing ``make check``. This will run ``zeno`` and, if it exists,
``zeno-mpi`` on various test cases and compare the output against the
output from a correctly built version. Floating-point values will be
allowed some tolerance to account for differences in compilers, machine
precision, etc.

.. _sec:runcode:

