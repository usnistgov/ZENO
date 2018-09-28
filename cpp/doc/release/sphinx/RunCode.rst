============
Command-line 
============

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

The code is run using:

``./zeno [OPTIONS]``

For example, the command to run ZENO for 1e6 random walks and 1e5
interior samples on an object described in the file ``obj.bod`` with the
results being written to the file ``out.csv`` is:

``./zeno -i obj.bod --num-walks=1000000 --num-interior-samples=100000 
--csv-output-file out.csv``

Command-line options may also be stored in a config file.  The example above
could also have been run by creating a file ``config.cfg`` containing

::
   
  input-file = obj.bod
  num-walks = 1000000
  num-interior-samples = 100000
  csv-output-file = out.csv

and running

``./zeno -c config.cfg``

If an option given in a config file is also given on the command line, the
command-line value overrides the config-file value.

Detailed descriptions of all options are below.

Required command-line options
-----------------------------

Input file
~~~~~~~~~~

The input or ``.bod`` file must contain a list of geometric primitives that
define the shape of the object. See :ref:`defineobj` for
details on the content of this file. Note that it can also contain
optional quantities as specified in :ref:`optinputs`.
The input file is specified via the command-line argument

``-i <input file name>``

.. _exterior-calculation-1:

Exterior calculation
~~~~~~~~~~~~~~~~~~~~

Either the number of exterior walks, the maximum relative standard
deviation of the capacitance, or the maximum relative standard
deviation of the mean electric polarizability must be specified. If
one of the relative standard deviations are specified, then the
calculation will continue performing walks until the relative standard
deviation is less than the specified value. One of the following
options is required for running the exterior calculation.

``--num-walks=<number of walks>``

``--max-rsd-capacitance=<maximum relative standard deviation of capacitance>``

``--max-rsd-polarizability=<maximum relative standard deviation of mean electric polarizability>``

.. _interior-calculation-1:

Interior calculation
~~~~~~~~~~~~~~~~~~~~

Either the number of interior samples or the maximum relative standard
deviation of the volume must be specified. If the relative standard
deviation is specified, then the calculation will continue performing
walks until the relative standard deviation is less than the specified
value. One of the following options is required for running the
interior calculation.

``--num-interior-samples=<number of interior samples>``

``--max-rsd-volume=<maximum relative standard deviation of volume>``

.. _sec:cmdline:

Description of command-line options
-----------------------------------

+---------+-------------------------------------+-----------------------+
| ``-h``, | ``--help``                          | Print help and exit   |
+---------+-------------------------------------+-----------------------+
| ``-V``, | ``--version``                       | Print version and     |
|         |                                     | exit                  |
+---------+-------------------------------------+-----------------------+
| ``-c``, | ``--config-file=string``            | Config file name      |
+---------+-------------------------------------+-----------------------+
| ``-i``, | ``--input-file=string``             | Input file name       |
|         |                                     |                       |
+---------+-------------------------------------+-----------------------+
|         | ``--csv-output-file=string``        | Write output in CSV   |
|         |                                     | format to the         |
|         |                                     | specified file in     |
|         |                                     | addition to           |
|         |                                     | displaying the        |
|         |                                     | regular output        |
+---------+-------------------------------------+-----------------------+
|         | ``--num-walks=int``                 | Number of exterior    |
|         |                                     | walks to perform      |
+---------+-------------------------------------+-----------------------+
|         | ``--num-interior-samples=int``      | Number of interior    |
|         |                                     | samples to take       |
+---------+-------------------------------------+-----------------------+
|         | ``--max-rsd-capacitance=double``    | Perform exterior      |
|         |                                     | walks until the       |
|         |                                     | relative standard     |
|         |                                     | deviation of the      |
|         |                                     | capacitance drops     |
|         |                                     | below this value.     |
|         |                                     | Relative standard     |
|         |                                     | deviation is defined  |
|         |                                     | as                    |
+---------+-------------------------------------+-----------------------+
|         | ``--max-rsd-polarizability=double`` | Perform exterior      |
|         |                                     | walks until the       |
|         |                                     | relative standard     |
|         |                                     | deviation of the mean |
|         |                                     | electric              |
|         |                                     | polarizability drops  |
|         |                                     | below this value.     |
|         |                                     | Relative standard     |
|         |                                     | deviation is defined  |
|         |                                     | as (Standard          |
|         |                                     | Deviation/Mean)       |
|         |                                     | :math:`\times` 100%   |
+---------+-------------------------------------+-----------------------+
|         | ``--max-rsd-volume=double``         | Take interior samples |
|         |                                     | until the relative    |
|         |                                     | standard deviation of |
|         |                                     | the volume drops      |
|         |                                     | below this value.     |
|         |                                     | Relative standard     |
|         |                                     | deviation is defined  |
|         |                                     | as                    |
+---------+-------------------------------------+-----------------------+
|         | ``--min-num-walks=int``             | Minimum number of     |
|         |                                     | exterior walks to     |
|         |                                     | perform when using    |
|         |                                     | max-relative standard |
|         |                                     | deviation stopping    |
|         |                                     | conditions            |
+---------+-------------------------------------+-----------------------+
|         | ``--min-num-interior-samples=int``  | Minimum number of     |
|         |                                     | interior samples to   |
|         |                                     | take when using       |
|         |                                     | max-relative standard |
|         |                                     | deviation stopping    |
|         |                                     | conditions            |
+---------+-------------------------------------+-----------------------+
|         | ``--max-run-time=double``           | Max time (in seconds) |
|         |                                     | that the program is   |
|         |                                     | allowed to run.  If   |
|         |                                     | this time is reached, |
|         |                                     | the computation will  |
|         |                                     | be stopped and the    |
|         |                                     | results computed so   |
|         |                                     | far will be displayed |
+---------+-------------------------------------+-----------------------+
|         | ``--num-threads=int``               | Number of threads to  |
|         |                                     | use (default=Number   |
|         |                                     | of logical cores)     |
+---------+-------------------------------------+-----------------------+
|         | ``--seed=INT``                      | Seed for the random   |
|         |                                     | number generator      |
|         |                                     | (default=Randomly     |
|         |                                     | set)                  |
+---------+-------------------------------------+-----------------------+
|         | ``--surface-points-file-string``    | Name of file for      |
|         |                                     | writing the surface   |
|         |                                     | points from exterior  |
|         |                                     | calculation           |
+---------+-------------------------------------+-----------------------+
|         | ``--interior-points-file=string``   | Name of file for      |
|         |                                     | writing the interior  |
|         |                                     | sample points         |
+---------+-------------------------------------+-----------------------+
|         | ``--print-counts``                  | Print statistics      |
|         |                                     | related to counts of  |
|         |                                     | hit points            |
+---------+-------------------------------------+-----------------------+
|         | ``--print-benchmarks``              | Print detailed RAM    |
|         |                                     | and timing            |
|         |                                     | information           |
+---------+-------------------------------------+-----------------------+

.. raw:: latex

   \addtocounter{table}{-1}

.. _input-file-1:

