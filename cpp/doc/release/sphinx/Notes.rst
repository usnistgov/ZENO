===================
Version Information
===================

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

Version 3
---------

3.x : Written in Fortran77 and developed at Stevens Institute of Technology by Dr. Marc
Mansfield and coworkers, this version can be downloaded at https://web.stevens.edu/zeno

Version 4
---------

4.0 : Takes version 3 and ports it to Fortran 2008

Version 5
---------

5.0.x : Implemented from scratch in C++, makes use of k-d trees and implements parallelism to improve performance
as documented in Ref. :cite:`Juba2016`

5.1.x : Replaces the k-d trees with axis-aligned bounding box data structure and includes 
the ability to read voxels in addition to spheres as input


