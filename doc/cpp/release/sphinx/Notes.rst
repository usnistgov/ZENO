================
Historical Notes
================

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

The development of the current version (v5) of the code was motivated by
the need to modernize the code base and to significantly speedup the
computation. ZENO up to version 3.x was written in Fortran77. These
versions were developed at Stevens Institute of Technology by Dr. Marc
Mansfield and coworkers; they can be downloaded from
https://web.stevens.edu/zeno For version 4, the code was ported to
Fortran 2008. For version 5, ZENO is implemented from scratch in C++.
This new version takes advantage of parallelism to deliver up to four
orders of magnitude speedup compared to the Fortran versions as
described in Ref. :cite:`Juba2016`. Note that the algorithms
in ZENO have been extensively tested previously, e.g.,
Refs. :cite:`Mansfield2008,Mansfield2001`.

