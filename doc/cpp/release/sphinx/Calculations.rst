============
Calculations
============

.. _Calculations:

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

The ZENO code is composed of two types of calculations: exterior and
interior.

Exterior calculation
--------------------

The exterior calculation focuses on the computation of electrical
properties including the capacitance, the electric polarizability
tensor, and the intrinsic conductivity. Once the electrical properties
are known, the hydrodynamic properties, including the hydrodynamic
radius and the intrinsic viscosity, can be precisely estimated by
invoking an electrostatic-hydrodynamic analogy as detailed in
Refs. :cite:`Douglas1995,Douglas1994,Hubbard1993`. Other
related properties are also determined.

To compute the aforementioned properties for an object requires the
solution of Laplace’s equation outside the object with appropriate
boundary conditions. This is efficiently accomplished by using a Monte
Carlo method, which involves (1) creating a launch sphere that
encloses the object, (2) launching random walks from the surface of
the launch sphere, and (3) determining the fate of such walks—if they
hit the object or go to infinity. These walks are exterior to the
object, hence the name for the calculation. Each random walk is
generated using a method called Walk on Spheres. This algorithm
requires generating a sphere for each step in the random walk. The
center of this sphere is located at the end of the current random
walk; the radius of the sphere is determined by finding the shortest
distance between the center of the sphere and the object. Finally, the
step in the walk is taken by randomly choosing a point on the surface
of the sphere. The process is then repeated. Since the size of spheres
will progressively get smaller as the object is approached, a cutoff
distance, known as the skin thickness, is required. Without a cutoff
distance, the algorithm would continue, at least theoretically,
indefinitely. As this is reminiscent of Zeno’s paradox of Tortoise and
Achilles, the code is named in Zeno’s honor. For more details on this
method refer to
Refs. :cite:`Douglas1995,Mansfield2008,Mansfield2001`.

Interior calculation
--------------------

The interior calculation determines the volume and the gyration tensor
for an object using a Monte Carlo method. Specifically, this calculation
involves generating random points within the same launch sphere as in
the exterior calculation. The location of these points can then be used
to approximate all of the relevant properties. For example, the volume
of the object is estimated by the fraction of points inside the object
multiplied by the the volume of the launch sphere. The interior
calculation is given its name since the points in the interior of the
object are essential for computing the properties.
