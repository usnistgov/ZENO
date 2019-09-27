==========
Validation
==========

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle


The algorithms in ZENO have been extensively tested previously, e.g.,
Refs. :cite:`Mansfield2008,Mansfield2001`. However, in order to 
validate the implementation of the current version of the code, 
as well as provide examples for users, five different systems were considered:
(1) two touching spheres of radius 1, (2) two touching spheres: one of
radius 1 and one of :math:`1/4`, (3) a cube, (4) a cuboid of 1 x 2 x 3,
(5) a “polymer" composed of 20 beads, (6) the protein lysozyme, 
(7) a torus with a minor radius of 1 and a major radius of 4 that
is represented by 416 spheres
(8) two touching spheres: one of radius 1 and one of 4 that has
been voxelized such that 32 voxels have a length of 1. 
For lysozyme, the pdb file 1LYD (see Refs. :cite:`pdbonline,pdb`) was used, and
a 5 :math:`\unicode{xC5}` sphere was placed at the center of each alpha
carbon of each amino acid. This procedure was previously used in
Ref. :cite:`Kang2004`. All ``.bod`` files used in testing
can be found in ``src/cpp/SelfTests``.

To validate the accuracy of ZENO, we use analytic values where possible 
to establish a true value, or ground truth. 
For systems composed of two spheres or a torus, the analytic values
are taken from Ref. :cite:`Mansfield2001`.
For the cube and cuboid, only analytic values
are known for the volume and eigenvalues of the gyration tensor.
In all other cases, ground truths do not exist so
the output associated with 1e11 walks and 1e11
interior samples was used as a surrogate.
In thses cases, the ground truth is taken as
the mean from three runs, and the accuracy is determined using the
expanded uncertainty (see below) from those three runs. Note that the
:math:`t` statistic is modified (:math:`t=4.302653`).
In cases where the representation will generate error, i.e. a torus
represented by a collection of spheres and two spheres that have been
voxelized, the results are compared to both analytic values where the
uncertainity is due to both 
longer runs of ZENO where the uncertainity is due to the algorithm
and the representation and the algorithm.

For each property, the following quantities were calculated: the mean
(:math:`\bar{y}=\sum_{i} y_{i}/N`), difference between ground truth
and mean (:math:`\Delta=|\mu-\bar{y}|`), relative
difference (:math:`{(\Delta / \mu) \cdot 100\,\%}`), standard deviation
(:math:`s=\sum_{i} (y_{i}-\bar{y})^{2}/(N-1)`), standard uncertainty
(:math:`s/\sqrt{N}`), expanded uncertainty (:math:`st/\sqrt{N}`), and
relative expanded uncertainty
(:math:`st/(\mu\sqrt{N}) \cdot 100\,\%`). For the equations, the
ground truth is represented by :math:`\mu`, the number of samples is
:math:`N`, and the value of observation :math:`i` of a property is
:math:`y_{i}`. The expanded uncertainty was computed using a 95 %
confidence interval (:math:`t=2.009575`). Each of the properties was
found to have a normal distribution. As a result, the mean and
standard deviation should define the distribution.

Results from the aforementioned tests can be found in the following
tables. Each test was run a total of 50 different times for 1e6 walks
and 1e6 interior samples, as well as 1e7 walks and 1e7 interior
samples. Note that these two quantities should be chosen carefully as
they affect the uncertainty of the calculation; choosing larger values
reduces the expanded uncertainty.


Two touching spheres of radius 1
--------------------------------

Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/two_spheres_1_1_1e+06.rst

Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/two_spheres_1_1_1e+07.rst


Two touching spheres one of radius 1 and one of radius 1/4
----------------------------------------------------------

Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/two_spheres_1_4_1e+06.rst

Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/two_spheres_1_4_1e+07.rst


Cube
----

Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/cube_1_1e+06.rst

Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/cube_1_1e+07.rst


Cuboid of 1 x 2 x 3
-------------------

Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/cuboid_1_2_3_1e+06.rst

Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/cuboid_1_2_3_1e+07.rst


Polymer
-------

Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/polymer_1e+06.rst

Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/polymer_1e+07.rst


Protein
-------

Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/1LYD_1e+06.rst

Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/1LYD_1e+07.rst


Torus with a minor radius of 1 and major radius of 4 sphere representation
--------------------------------------------------------------------------

Comparison with same representation; Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/spheres_torus_1_4_16_1e+06_surrogate.rst

Comparison with same representation; Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/spheres_torus_1_4_16_1e+07_surrogate.rst

Comparison to analytic values such that representation generates error; Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/spheres_torus_1_4_16_1e+06_analytic.rst

Comparison to analytic values such that representation generates error; Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/spheres_torus_1_4_16_1e+07_analytic.rst


Two touching spheres one of radius 1 and one of radius 4 voxelized
------------------------------------------------------------------

Comparison with same representation; Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/voxels_two_spheres_1_4_32_1e+06_surrogate.rst

Comparison with same representation; Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/voxels_two_spheres_1_4_32_1e+07_surrogate.rst

Comparison to analytic values such that representation generates error; Number of walks = 1e6; Number of interior samples = 1e6

.. include:: validation_tables/voxels_two_spheres_1_4_32_1e+06_analytic.rst

Comparison to analytic values such that representation generates error; Number of walks = 1e7; Number of interior samples = 1e7

.. include:: validation_tables/voxels_two_spheres_1_4_32_1e+07_analytic.rst

