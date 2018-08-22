==========
Input file
==========

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

The input file, also known as a ``.bod`` file, contains the description
of the object and some additional input parameters.

.. _defineobj:

Defining the object
-------------------

The object of interest must be described by a collection of spheres,
which may or may not be overlapping. The shape of the object is
defined in the ``.bod`` file. This file should contain lines of

``SPHERE x y z r``

where ``x``, ``y``, and ``z`` are the :math:`x`, :math:`y`, and
:math:`z` positions of the center of the sphere, respectively. ``r``
is the radius. For example, a ``.bod`` file that contains the
following describes an object composed of two spheres: one of radius 2
at :math:`x=0`, :math:`y=0`, and :math:`z=1` and one of radius 3 at
:math:`x=0`, :math:`y=0`, and :math:`z=-1`.

``SPHERE 0 0 1 2``
``SPHERE 0 0 -1 3``

.. _optinputs:

Optional inputs
---------------

Launch radius
~~~~~~~~~~~~~

+-------------------+-----------------------------------+
| Command:          | ``rlaunch double``                |
+-------------------+-----------------------------------+
| Explanation:      | Sets the radius, which is radius  |
|                   | of the sphere from which random   |
|                   | walks are launched. The radius    |
|                   | must be large enough to enclose   |
|                   | the entire object.                |
+-------------------+-----------------------------------+
| Default value:    | The smallest radius that encloses |
|                   | the smallest axis-aligned         |
|                   | bounding-box of the object.       |
+-------------------+-----------------------------------+
| Example:          | ``rlaunch 20`` means that the     |
|                   | launch radius is 20.              |
+-------------------+-----------------------------------+

Skin thickness
~~~~~~~~~~~~~~

+-------------------+-----------------------------------+
| Command:          | ``st double``                     |
+-------------------+-----------------------------------+
| Explanation:      | Sets the skin thickness. A random |
|                   | walker is assumed to have hit the |
|                   | surface of the object if the      |
|                   | distance between the surface and  |
|                   | the walker is less than the skin  |
|                   | thickness.                        |
+-------------------+-----------------------------------+
| Default value:    | 1e-6 times the launch radius      |
+-------------------+-----------------------------------+
| Example:          | ``st 0.01`` means that the skin   |
|                   | thickness is 0.01.                |
+-------------------+-----------------------------------+

Units for length
~~~~~~~~~~~~~~~~

+-------------------+-----------------------------------+
| Command:          | ``hunits double string``          |
+-------------------+-----------------------------------+
| Explanation:      | Specifies the units for the       |
|                   | length for all objects.           |
+-------------------+-----------------------------------+
| Options:          | The string can take the following |
|                   | values:                           |
|                   |                                   |
|                   | -  ``m`` (meters)                 |
|                   |                                   |
|                   | -  ``cm`` (centimeters)           |
|                   |                                   |
|                   | -  ``nm`` (nanometers)            |
|                   |                                   |
|                   | -  ``A`` (Angstroms)              |
|                   |                                   |
|                   | -  ``L`` (generic or unspecified  |
|                   |    length units)                  |
+-------------------+-----------------------------------+
| Default value:    | 1 ``L``                           |
+-------------------+-----------------------------------+
| Example:          | ``hunits 10 cm`` means that a     |
|                   | length of 1 for an object is      |
|                   | equivalent to 10 cm.              |
+-------------------+-----------------------------------+

Temperature
~~~~~~~~~~~

+-------------------+-----------------------------------+
| Command:          | ``temp double string``            |
+-------------------+-----------------------------------+
| Explanation:      | Specifies the temperature, which  |
|                   | is used for computing the         |
|                   | diffusion coefficient.            |
+-------------------+-----------------------------------+
| Options:          | The string can take the following |
|                   | values:                           |
|                   |                                   |
|                   | -  ``C`` (Celsius)                |
|                   |                                   |
|                   | -  ``K`` (Kelvin)                 |
+-------------------+-----------------------------------+
| Default value:    | None                              |
+-------------------+-----------------------------------+
| Example:          | ``temp 20 C`` means that the      |
|                   | temperature is                    |
|                   | 20\ :math:`^\circ`\ C.            |
+-------------------+-----------------------------------+

Mass
~~~~

+-------------------+-----------------------------------+
| Command:          | ``mass double string``            |
+-------------------+-----------------------------------+
| Explanation:      | Specify the mass of the object,   |
|                   | which is used for computing the   |
|                   | intrinsic viscosity in            |
|                   | conventional units and the        |
|                   | sedimentation coefficient.        |
+-------------------+-----------------------------------+
| Options:          | The string can take the following |
|                   | values:                           |
|                   |                                   |
|                   | -  ``Da`` (Daltons)               |
|                   |                                   |
|                   | -  ``kDa`` (kiloDaltons)          |
|                   |                                   |
|                   | -  ``g`` (grams)                  |
|                   |                                   |
|                   | -  ``kg`` (kilograms)             |
+-------------------+-----------------------------------+
| Default value:    | None                              |
+-------------------+-----------------------------------+
| Example:          | ``mass 2 g`` means that the mass  |
|                   | of the object is 2 grams.         |
+-------------------+-----------------------------------+

Solvent viscosity
~~~~~~~~~~~~~~~~~

+-------------------+-----------------------------------+
| Command:          | ``viscosity double string``       |
+-------------------+-----------------------------------+
| Explanation:      | Specify the solvent viscosity,    |
|                   | which is used for computing the   |
|                   | diffusion coefficient, the        |
|                   | friction coefficient, and the     |
|                   | sedimentation coefficient.        |
+-------------------+-----------------------------------+
| Options:          | The string can take the following |
|                   | values:                           |
|                   |                                   |
|                   | -  ``p`` (poise)                  |
|                   |                                   |
|                   | -  ``cp`` (centipoise)            |
+-------------------+-----------------------------------+
| Default value:    | None                              |
+-------------------+-----------------------------------+
| Example:          | ``viscosity 2 cp`` means that the |
|                   | solvent has a viscosity of 2      |
|                   | centipoise.                       |
+-------------------+-----------------------------------+

Buoyancy factor
~~~~~~~~~~~~~~~

+-------------------+-----------------------------------+
| Command:          | ``bf double``                     |
+-------------------+-----------------------------------+
| Explanation:      | Specify the buoyancy factor,      |
|                   | which is used for computing the   |
|                   | sedimentation coefficient.        |
+-------------------+-----------------------------------+
| Default value:    | None                              |
+-------------------+-----------------------------------+
| Example:          | ``bf 2`` means that the buoyancy  |
|                   | factor is 2.                      |
+-------------------+-----------------------------------+

