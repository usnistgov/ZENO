Validation Reproducibility
==========================

The data tables in the Validation section of the documentation can be reproduced from scratch using the ``run_validation.py`` script.  This script performs three tasks: generate the surrogate ground truth data, run the correctness tests, and generate the data tables.  These tasks are described in more detail in the following sections.

The tasks in the ``run_validation.py`` script each perform their function by calling several other scripts.  A user could perform a subset of a task by calling the other scripts directly.  The syntax for eacxh script is given by the ``--help`` option.

Generate Surrogate Ground Truth
-------------------------------

This task generates data sets using a very large number of samples for certain gemoetries that do not have an analytic ground truth available.  These data sets can be used as a surrogate ground truth for purposes such as testing convergence and checking for regressions.  It is therefore recommended to reuse the surrogate ground truth data from the original version of the program, rather than generating new, to avoid the possibility of gradual drift in the program output between versions.  However, a user may wish to regenerate the surrogate gound truth data in order to verify the reported values.

Note that this task could potentially take several days of computation, depending on the system.

This task is commented out by default in the ``run_validation.py`` script.

Run Correctness Tests
---------------------

This task runs the correctness tests on the test data sets.  Each data set is tested with a lower and a higher number of samples and several different random number seeds.  The results of these tests can be used to validate the correctness of the program and get a sense of its convergence with different numbers of samples.  Note that for tests without an analytic ground truth available, correctness can only be checked by testing for regressionss, not by testing if the values are actually correct.

Generate Data Tables
--------------------

This task generates the data tables for the Validation section of the documentation from the ground truth data (surrogate and analytic) and the correctness test results.  It computes the required statistics, converts the tables to RST format at the appropriate precisions, copies the tables to the appropriate location in the documentation source tree, and uses Sphinx to regenerate the documentation.  Note that this final step will fail if Sphinx is not installed.

