Correctness tests
-----------------

The correctness tests can be reproduced by running the following scripts.

``./generate_correctness_tests.py``

Generate a Zeno ``.cfg`` file for each test.  The tests are based on the ``.bod`` files in ``correctness/bod_files``.  The ``.cfg`` files are written to ``correctness/config_files``.

``./run_correctness_tests.py``

Run the tests specified by the ``.cfg`` files in ``correctness/config_files``.  The results of the tests are written to ``.csv`` files in ``correctness/raw_output``.

``./combine_correctness_tests.py``

Combine the ``.csv`` files in ``correctness/raw_output`` for runs of the same test case with different seeds into a single ``.csv`` file.  The combined ``.csv`` files are written to ``correctness/combined_output``.

Surrogate ground truth
----------------------

The surrogate ground truth can be reproduced in the same way as the correctness tests.  Simply follow the directions above, replacing ``correctness`` with ``surrogate_ground_truth``.  Be aware that running the surrogate ground truth tests can take several days, depending on your hardware.
