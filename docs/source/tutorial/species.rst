.. _species:

===========================================
species
===========================================

Some atoms are special:

- atoms with different starting_magnetization
- atoms with different U values
- atoms with special basis set

.. code-block:: python

    atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype = 'U20'))
    atoms.arrays['species'][1] = 'Fe1'
