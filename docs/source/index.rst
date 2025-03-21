#########
Equations
#########

The incompressible Navier-Stokes equations in three-dimensional cylindrical coordinates are derived `here <https://github.com/NaokiHori/SimpleTCSolver>`_.

Due to the azimuthal homogeneity:

.. math::

    &
    \vel{2}
    \equiv
    0,

    &
    \pder{q}{\gcs{2}}
    \equiv
    0,

the equations can be further simplified, which are given below.

****************************
Incompressibility constraint
****************************

.. math::

    \frac{1}{J}
    \pder{}{\gcs{1}}
    \left(
        \frac{J}{\sfact{1}} \vel{1}
    \right)
    +
    \frac{1}{J}
    \pder{}{\gcs{3}}
    \left(
        \frac{J}{\sfact{3}} \vel{3}
    \right)
    =
    0.

This is enforced by the SMAC method.

************
Mass balance
************

.. math::

   \pder{\density}{t}
   +
   \frac{1}{J}
   \pder{}{\gcs{1}}
   \left(
      \frac{J}{\sfact{1}}
      \density
      \vel{1}
   \right)
   +
   \frac{1}{J}
   \pder{}{\gcs{3}}
   \left(
      \frac{J}{\sfact{3}}
      \density
      \vel{3}
   \right)
   =
   0.

This is solved by means of the volume-of-fluid method.

****************
Momentum balance
****************

.. math::

    \momtemp{1}
    =
    &
    \momadv{1}{1}
    \momadv{3}{1}

    &
    \mompre{1}

    &
    \momdif{1}{1}
    \momdif{3}{1}

    &
    \momdifx,

.. math::

    \momtemp{3}
    =
    &
    \momadv{1}{3}
    \momadv{3}{3}

    &
    \mompre{3}

    &
    \momdif{1}{3}
    \momdif{3}{3}.

This is solved by the energy-conserving scheme.

*******************
Shear-stress tensor
*******************

The shear-stress tensor for Newtonian liquids is defined as

.. math::

    \sst{i}{j}
    \equiv
    \viscosity
    \vgt{i}{j}
    +
    \viscosity
    \vgt{j}{i},

where :math:`\vgt{i}{j}` is a second-order tensor representing the gradient of velocity vector:

.. math::

    \sum_i
    \sum_j
    \vec{e}_i
    \otimes
    \vec{e}_j
    \vgt{i}{j}.

The components are

.. math::

    \begin{pmatrix}
        \vgt{1}{1} & \vgt{2}{1} & \vgt{3}{1} \\
        \vgt{1}{2} & \vgt{2}{2} & \vgt{3}{2} \\
        \vgt{1}{3} & \vgt{2}{3} & \vgt{3}{3} \\
    \end{pmatrix}
    =
    \begin{pmatrix}
        \frac{1}{\sfact{1}}
        \pder{\vel{1}}{\gcs{1}}
        &
        0
        &
        \frac{1}{\sfact{3}}
        \pder{\vel{1}}{\gcs{3}}
        \\
        0
        &
        \frac{1}{J}
        \pder{}{\gcs{1}}
        \left(
            \frac{J}{\sfact{1}}
        \right)
        \vel{1}
        &
        0
        \\
        \frac{1}{\sfact{1}}
        \pder{\vel{3}}{\gcs{1}}
        &
        0
        &
        \frac{1}{\sfact{3}}
        \pder{\vel{3}}{\gcs{3}}
    \end{pmatrix}.

The dynamic viscosity is computed using the cell-faced volume-of-fluid value.

