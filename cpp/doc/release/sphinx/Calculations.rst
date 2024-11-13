============
Calculations
============

.. _Calculations:

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

The ZENO code is composed of three types of calculations: exterior, 
interior, and virial.

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

Virial calculation
--------------------

The virial calculation evaluates the virial coefficients, which appear in the power-series expansion of the pressure as a function density. 
The calculations formally yield the gas-phase virial coefficients, but these can be interpreted as coefficients of a series for the osmotic pressure 
as a function of solute concentation, where the solutes interact according to a solvent-averaged pair potential. Unlike the exterior and interior 
calculations, which are performed for a single object, the virial calculations involve multiple objects: calculation of the second virial 
coefficient involves two objects, third virial coefficient three objects, and so on.  Coefficients :math:`B_N` to any order :math:`N` may be 
calculated by ZENO, but the required computational effort grows rapidly with :math:`N`, such that practically a limit is found at about :math:`N = 7`.  
Even then, usually just :math:`B_2` and maybe :math:`B_3` are of primary interest.

The virial calculation employs the Mayer-sampling Monte Carlo algorithm.  This entails random sampling of configurations of the objects relative 
to each other, weighted by their importance to the averages being collected (see Ref. :cite:`Bansal2022` for details).  If the objects have 
internal degrees of freedom (e.g., conformations of a chain molecule), these are sampled as well. The algorithm yields the virial coefficient 
relative to a known reference system of hard spheres. The diameter of these spheres should be selected to optimize the calculation, based on a 
rough measure of the objects of interest (e.g., radius of gyration of a chain molecule). The calculations become inefficient (but are still in 
principle correct) when applied to objects that are severely non-spherical (e.g., a long, rigid rod). The objects may interact with each other 
according to several choices of interatomic potential; also, intramolecular potentials (e.g., bond bending or stretching) may be specified.
