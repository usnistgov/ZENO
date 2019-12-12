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

Dependencies
------------

The only essential external library is the SPRNG :cite:`Mascagni2000` random number library.
This requires compilation, which is described in the next section.

Other required libraries are zlib and a threading library.  These should
be provided by your Operating System distribution.

Building SPRNG
~~~~~~~~~~~~~~

SPRNG can be obtained from http://sprng.org.  The required version is 5.0.

When configuring the SPRNG library, you will have the option to include
support for Fortran and MPI.  Neither of these are required for ZENO (even
if you plan to use MPI with ZENO), so you may build SPRNG without these.

Since the core ZENO functionality is built as a shared library, SPRNG must be
configured to be compiled with the ``-fPIC`` flag.

These example commands should build and install SPRNG in your home directory:

::

   tar xjf sprng5.tar.bz2
   mkdir sprng5-build && cd sprng5-build
   CXXFLAGS=-fPIC ../sprng5/configure --with-fortran=no --with-mpi=no --prefix=$HOME/sprng5-install
   make && make install


Building ZENO
-------------

Compiling the code requires cmake version 3.1 or greater.

To compile, you should first make a build directory.  Then change to the
build directory and run ``cmake`` on the ``ZENO/cpp`` directory.  You will
need to set the ``SPRNG_INCLUDE_DIR`` and ``SPRNG_LIBRARY`` variables based
on where you installed SPRNG.

You should then be able to build the executable ``zeno`` by simply
typing ``make``.  Optionally, you may also install the executable into
the directory specified in cmake by typing ``make install``.

MPI support
~~~~~~~~~~~

MPI support is included, but is optional. If MPI libraries are installed
on your system, you may be able to build zeno with MPI support by running
cmake with ``-DUSE_MPI=ON``.

Cmake will attempt to automatically find the MPI libraries on your system.
If this does not work, you will need to manually set the various cmake MPI
variables to the values appropriate for your installation.

Simple check
~~~~~~~~~~~~

Once the code has been compiled, you can perform a quick self-test by
typing ``make check``. This will run ``zeno`` 
on various test cases and compare the output against the
output from a correctly built version. Floating-point values will be
allowed some tolerance to account for differences in compilers, machine
precision, etc.

Example
~~~~~~~

These example commands should build and install ZENO in your home directory
(without MPI support), assuming SPRNG has also been installed in your home
directory:

::

   git clone https://github.com/usnistgov/ZENO.git
   mkdir zeno-build && cd zeno-build
   cmake -DSPRNG_INCLUDE_DIR=$HOME/sprng5-install/include -DSPRNG_LIBRARY=$HOME/sprng5-install/lib/libsprng.a -DCMAKE_INSTALL_PREFIX=$HOME/zeno-install ../ZENO/cpp
   make && make install
   make check


Modifying the code
------------------

If you modify the source code, some external utilities may be required.
``Gengetopt`` is required to modify command-line parameters, while
``Bisonc++`` and ``Flexc++`` are required to modify the input-file
format. These can be obtained from the sources:

-  https://www.gnu.org/software/gengetopt

-  https://gitlab.com/fbb-git/bisoncpp

-  https://gitlab.com/fbb-git/flexcpp

You will then need to enable regeneration of the comand-line parser and
input-file parser with the cmake options ``-DGEN_CMD_PARSER=ON`` and
``-DGEN_BOD_PARSER=ON``.

.. _sec:runcode:

