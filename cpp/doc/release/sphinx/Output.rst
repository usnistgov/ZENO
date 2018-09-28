======
Output
======

.. role:: raw-latex(raw)
   :format: latex
..

.. raw:: latex

   \maketitle

Depending on the quantity, the requirements for calculations are
different. Some quantities require the exterior calculation, some the
interior, and some both. Additionally, some quantities are direct
outputs, while others are indirect—requiring only the direct outputs
coupled with algebraic expressions. Finally, some special quantities
require optional input such as the temperature. Each of these
dependencies are listed.

Capacitance
-----------

+--------------+-----------------------------------+
| Equation:    | :math:`C = tR`                    |
+--------------+-----------------------------------+
| Calculation: | Direct exterior                   |
+--------------+-----------------------------------+
| Explanation: | :math:`t` is the fraction of      |
|              | random walks that hit the object  |
|              | as opposed to go to infinity, and |
|              | :math:`R` is the radius of the    |
|              | launch sphere.                    |
+--------------+-----------------------------------+
| Uncertainty: | Determined by directly estimating |
|              | the variance of :math:`t` and     |
|              | then applying propagation of      |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length                            |
+--------------+-----------------------------------+

Electric polarizability tensor
------------------------------

+--------------+-----------------------------------+
| Equation:    | :math:`\mathbf{\alpha}` is        |
|              | expressed in                      |
|              | Ref. :cite:`Mansfield2001`.       |
|              |                                   |
+--------------+-----------------------------------+
| Calculation: | Direct exterior                   |
+--------------+-----------------------------------+
| Explanation: | Several different counting        |
|              | variables from the exterior       |
|              | calculation are combined to       |
|              | evaluate this quantity.           |
+--------------+-----------------------------------+
| Uncertainty: | Determined by directly estimating |
|              | the variances of :math:`t`,       |
|              | :math:`u`, :math:`v`, and         |
|              | :math:`w` from                    |
|              | Ref. :cite:`Mansfield2001`        |
|              | and then applying propagation of  |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length cubed                      |
+--------------+-----------------------------------+

Eigenvalues of electric polarizability tensor
---------------------------------------------

+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Explanation: | The eigenvalues of the previously |
|              | computed electric polarizability  |
|              | tensor :math:`\mathbf{\alpha}`    |
|              | are determined.                   |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length cubed                      |
+--------------+-----------------------------------+

Mean electric polarizability
----------------------------

+--------------+-----------------------------------+
| Equation:    | :math:`\langle \alpha \rangle =   |
|              | \mathrm{Tr}(\mathbf{\alpha})/3`   |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Explanation: | The trace of the previously       |
|              | computed electric polarizability  |
|              | tensor :math:`\mathbf{\alpha}` is |
|              | computed and then divided by      |
|              | three.                            |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length cubed                      |
+--------------+-----------------------------------+

Intrinsic conductivity
----------------------

+--------------+-----------------------------------+
| Equation:    | :math:`[\sigma]_\infty = \langle  |
|              | \alpha \rangle/V`                 |
+--------------+-----------------------------------+
| Explanation: | :math:`\mathbf{\alpha}` is the    |
|              | electric polarizability tensor;   |
|              | :math:`V` is the volume.          |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior and interior    |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | None                              |
+--------------+-----------------------------------+

Volume
------

+--------------+-----------------------------------+
| Equation:    | :math:`V= p \frac{4}{3} \pi R^{3}`|
|              |                                   |
+--------------+-----------------------------------+
| Explanation: | :math:`p` is the fraction of      |
|              | points inside the object;         |
|              | :math:`R` is the radius of the    |
|              | launch sphere.                    |
+--------------+-----------------------------------+
| Calculation: | Direct interior                   |
+--------------+-----------------------------------+
| Uncertainty: | Determined by directly estimating |
|              | the variance of :math:`p` and     |
|              | then applying propagation of      |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length cubed                      |
+--------------+-----------------------------------+

Gyration tensor
---------------

+--------------+-----------------------------------+
| Equation:    | :math:`\mathbf{S}` is expressed   |
|              | in                                |
|              | Ref. :cite:`Theodorou1985`        |
+--------------+-----------------------------------+
| Explanation: | Several different counting        |
|              | variables from the interior       |
|              | calculation are combined to       |
|              | evaluate this quantity.           |
+--------------+-----------------------------------+
| Calculation: | Direct interior                   |
+--------------+-----------------------------------+
| Uncertainty: | Determined by directly estimating |
|              | the variance of sums and the sums |
|              | of products of interior sample    |
|              | point coordinates, and then       |
|              | applying propagation of           |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length squared                    |
+--------------+-----------------------------------+

Eigenvalues of gyration tensor
------------------------------

+--------------+-----------------------------------+
| Calculation: | Indirect interior                 |
+--------------+-----------------------------------+
| Explanation: | The eigenvalues of the previously |
|              | computed gyration tensor          |
|              | :math:`\mathbf{S}` are            |
|              | determined.                       |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length squared                    |
+--------------+-----------------------------------+

Capacitance of a sphere of the same volume
------------------------------------------

+--------------+----------------------------------------------+
| Equation:    | :math:`C_0 = \left(3V/(4\pi)\right)^{1/3}`   |
+--------------+----------------------------------------------+
| Calculation: | Indirect interior                            |
+--------------+----------------------------------------------+
| Explanation: | :math:`V` is the volume of the object.       |
+--------------+----------------------------------------------+
| Uncertainty: | Determined via propagation of uncertainties. |
+--------------+----------------------------------------------+
| Units:       | Length                                       |
+--------------+----------------------------------------------+

Hydrodynamic radius
-------------------

+--------------+-----------------------------------+
| Equation:    | :math:`R_{h}=q_{R_{h}}C`          |
+--------------+-----------------------------------+
| Explanation: | :math:`q_{R_{h}}\approx 1`, and   |
|              | :math:`C` is the capacitance.     |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties assuming the        |
|              | standard deviation of             |
|              | :math:`q_{R_{h}}` is              |
|              | :math:`0.01`.                     |
+--------------+-----------------------------------+
| Units:       | Length                            |
+--------------+-----------------------------------+

Prefactor relating average polarizability to intrinsic viscosity
----------------------------------------------------------------

+--------------+-----------------------------------+
| Equation:    | :math:`q_\eta` varies slowly with |
|              | shape and is expressed in         |
|              | Ref. :cite:`Mansfield2008`        |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Explanation: | The electric polarizability       |
|              | tensor plus a complicated Padé    |
|              | approximate is used to determine  |
|              | this quantity.                    |
+--------------+-----------------------------------+
| Uncertainty: | :math:`0.015q_\eta`               |
+--------------+-----------------------------------+
| Units:       | None                              |
+--------------+-----------------------------------+

Viscometric radius
------------------

+--------------+-----------------------------------+
| Equation:    | :math:`R_{v}= (3 q_\eta \langle   |
|              | \alpha \rangle/(10 \pi))^{1/3}`   |
+--------------+-----------------------------------+
| Explanation: | :math:`q_\eta` is the prefactor   |
|              | for the intrinsic viscosity, and  |
|              | :math:`\langle \alpha \rangle` is |
|              | the mean polarizability.          |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length                            |
+--------------+-----------------------------------+

Intrinsic viscosity
-------------------

+--------------+-----------------------------------+
| Equation:    | :math:`[\eta]=q_\eta              |
|              | [\sigma]_\infty`                  |
+--------------+-----------------------------------+
| Explanation: | :math:`q_\eta` is a prefactor,    |
|              | and :math:`[\sigma]_\infty` is    |
|              | the intrinsic conductivity.       |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior and interior    |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | None                              |
+--------------+-----------------------------------+

Intrinsic viscosity with mass units
-----------------------------------

+--------------+-----------------------------------+
| Equation:    | :math:`[\eta]_{m}=q_\eta \langle  |
|              | \alpha\rangle/m`                  |
+--------------+-----------------------------------+
| Explanation: | :math:`q_\eta` is the prefactor,  |
|              | :math:`\alpha` is the             |
|              | polarizability tensor, and        |
|              | :math:`m` is the specified mass.  |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length cubed / mass               |
+--------------+-----------------------------------+
| Requirements:| Specified mass.                   |
+--------------+-----------------------------------+

Friction coefficient
--------------------

+--------------+-----------------------------------+
| Equation:    | :math:`f = 6\pi\eta R_h`          |
+--------------+-----------------------------------+
| Explanation: | :math:`\eta` is the solvent       |
|              | viscosity, and :math:`R_h` is the |
|              | hydrodynamic radius.              |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Mass / time                       |
+--------------+-----------------------------------+
| Requirements:| Specified length scale and        |
|              | solvent viscosity.                |
+--------------+-----------------------------------+

Diffusion coefficient
---------------------

+--------------+-----------------------------------+
| Equation:    | :math:`D = k_{B}T/f`              |
+--------------+-----------------------------------+
| Explanation: | :math:`k_{B}` is the Boltzmann    |
|              | constant, :math:`T` is the        |
|              | temperature, and :math:`f` is the |
|              | friction coefficient.             |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Length squared / time             |
+--------------+-----------------------------------+
| Requirements:| Specified length scale, solvent   |
|              | viscosity, and temperature.       |
+--------------+-----------------------------------+

Sedimentation coefficient
-------------------------

+--------------+-----------------------------------+
| Equation:    | :math:`s = mb/f`                  |
+--------------+-----------------------------------+
| Explanation: | :math:`m` is the mass, :math:`b`  |
|              | is the buoyancy factor, and       |
|              | :math:`f` is the friction         |
|              | coefficient.                      |
+--------------+-----------------------------------+
| Calculation: | Indirect exterior                 |
+--------------+-----------------------------------+
| Uncertainty: | Determined via propagation of     |
|              | uncertainties.                    |
+--------------+-----------------------------------+
| Units:       | Time                              |
+--------------+-----------------------------------+
| Requirements:| Specified length scale, solvent   |
|              | viscosity, mass, and buoyancy     |
|              | factor.                           |
+--------------+-----------------------------------+

