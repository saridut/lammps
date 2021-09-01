.. index:: pair_style sciortino

pair_style sciortino command
=====================

Accelerator Variants: *none*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style sciortino

Examples
""""""""

.. code-block:: LAMMPS

   pair_style sciortino
   pair_coeff * * AB.sciortino A B

Description
"""""""""""

The *sciortino* style computes a 3-body :ref:`Sciortino <Sciortino>`
potential for the energy E of a system of atoms as

.. math::

   E & =  \sum_i \sum_{j > i} \phi_2 (r_{ij}) +
          \sum_i \sum_{j \neq i} \sum_{k > j}
          \phi_3 (r_{ij}, r_{ik}) \\
  \phi_2(r_{ij}) & =  4 \epsilon_{ij} \left[ (\frac{\sigma_{ij}}{r_{ij}})^{2m_{ij}} -
                    (\frac{\sigma_{ij}}{r_{ij}})^{m_{ij}} \right], 
                    \quad r_{ij} < r_{\mathrm{cut},ij} \\
  \phi_3(r_{ij},r_{ik}) & = \lambda_{ijk} \epsilon_{ijk} 
    \hat{\phi_2}^{(2b)}\left(r_{ij}\right) \hat{\phi_2}^{(2b)}\left(r_{ik}\right)

where

.. math::

   \hat{\phi_2}^{(2b)}\left(r_{ij}\right) = \begin{cases} 1 & r_{ij} \le r_{\mathrm{min}}\\
        -\frac{\phi_2(r_{ij})}{\epsilon_{ij}} & \mathrm{otherwise} \end{cases},

:math:`\phi_2` is a two-body term and :math:`\phi_3` is a
three-body term.  The summations in the formula are over all neighbors J
and K of atom I within a cutoff distance :math:`r_{\mathrm{cut}}`.

Only a single pair_coeff command is used with the *sciortino* style which
specifies a Sciortino potential file with parameters for all
needed elements.  These are mapped to LAMMPS atom types by specifying
N additional arguments after the filename in the pair_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of Sciortino elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine a file AB.sciortino has Sciortino values for
A and B.  If your LAMMPS simulation has 4 atoms types and you want
the first 3 to be A, and the fourth to be B, you would use the following
pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * AB.sciortino A A A B

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three A arguments map LAMMPS atom types 1,2,3 to the A
element in the Sciortino file.  The final B argument maps LAMMPS atom type 4
to the B element in the Sciortino file.  If a mapping value is specified as
NULL, the mapping is not performed.  This can be used when a *sciortino*
potential is used as part of the *hybrid* pair style.  The NULL values
are placeholders for atom types that will be used with other
potentials.

Sciortino files in the *potentials* directory of the LAMMPS
distribution have a ".sciortino" suffix.  Lines that are not blank or
comments (starting with #) define parameters for a triplet of
elements.  The parameters in a single entry correspond to the two-body
and three-body coefficients in the formula above:

* element 1 (the center atom in a 3-body interaction)
* element 2
* element 3
* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* :math: `m`
* :math: `r_{\mathrm{cut}}` (distance units)
* :math:`\lambda`

The :math:`\lambda` parameter is used exclusively for three-body interactions,
the rest are used both for two-body and three-body terms. Both :math:`\lambda`
and :math:`\epsilon` parameters must be defined for all interactions, though usually
they are constants across all triplets. The non-annotated parameters are unitless.
The parameter :math:`r_{\mathrm{min}}` is the minimum of the potential
:math:`\phi_2` for the appropriate pair, this value will be calculated
internally and need not be provided by the user.

The Sciortino potential file must contain entries for all the
elements listed in the pair_coeff command.  It can also contain
entries for additional elements not being used in a particular
simulation; LAMMPS ignores those entries.

For a single-element simulation, using the Sciortino potential do not make
much sense, though it can be used. In such case only a single entry is required
(e.g. AAA).  For a two-element simulation, the file must contain 8
entries (for AAA, AAB, ABA, ABB, BAA, BAB, BBA, BBB), that
specify Sciortino parameters for all permutations of the two elements
interacting in three-body configurations.  Thus for 3 elements, 27
entries would be required, etc.

As annotated above, the first element in the entry is the center atom
in a three-body interaction.  Thus an entry for ABB means a A atom
with 2 B atoms as neighbors.  The parameter values used for the
two-body interaction come from the entry where the second and third
elements are the same.  Thus the two-body parameters for A
interacting with B, comes from the ABB entry.  The three-body
parameters can in principle be specific to the three elements of the
configuration. In the :ref:`Sciortino <Sciortino>`, however, the three-body parameters
are defined in terms of  two sets of pair-wise
parameters, corresponding to the `ij` and `ik` pairs, where `i` is the
center atom. The user must ensure that the correct combining rule is
used to calculate the values of the three-body parameters for
alloys. Note also that the function :math:`\phi_3` contains two exponential
screening factors with parameter values from the ij pair and ik
pairs. 
Since the order of the two neighbors is arbitrary, the three-body parameters for
entries BAB and BBA should be the same.  Similarly, the two-body parameters for
entries ABB and BAA should also be the same.  The parameters :math:`\sigma` and
:math:`\m` need to be specified only for two-body interactions in entries whose
second and third element are same (e.g. ABB). In rest of the entries the values
are not used and can be set to zero. Furthermore, the :math:`r_{\mathrm{cut}}`
need not be specified between like atoms (e.g. AAA or BBB) and can be set to
zero. For like atoms, :math:`\phi_2` is shifted such that it is purely
repulsive; no shifting is performed for unlike atoms.


Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as
described above from values in the potential file.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, 
since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the MANYBODY package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The Sciortino potential files provided with LAMMPS (see the
potentials directory) are parameterized for *lj* :doc:`units <units>`.
You can use the Sciortino potential with any LAMMPS units, but you would need
to create your own Sciortino potential file with coefficients listed in the
appropriate units if your simulation does not use *lj* units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Sciortino:

**(Sciortino)** Sciortino, Eur Phys J E, 40, 3 (2017).
