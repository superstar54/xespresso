.. _pw:

===========================================
pw
===========================================
The :class:`~xespresso.xespresso.Espresso` object is used to perform a PW calculation.


A example of PW calculation of H2 molecule. One need to define the following inputs:

- ASE Atoms.
- kpoints.
- pseudopotentials.
- input_data.
- queue.

.. code-block:: python

    from ase.build import molecule
    from xespresso import Espresso

    # define H2 molecule
    h2 = molecule("H2")
    h2.cell = [8, 8, 8]
    h2.pbc = [True, True, True]
    # set pseudo-potential
    pseudopotentials = {
        "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
    }
    calc = Espresso(
        pseudopotentials=pseudopotentials,
        label="calculations/scf/h2",
        ecutwfc=20,
        occupations="smearing",
        degauss=0.03,
        kpts=(1, 1, 1), # kpoints
        debug=True,
    )
    h2.calc = calc
    e = h2.get_potential_energy()
    print("Energy: {0:1.4f}".format(e))
    assert np.isclose(e, -31.4454)

.. _nscf:

nscf
===========================================


A example of nscf calculation following the above one.

.. code-block:: python

    from xespresso.post.nscf import EspressoNscf
    nscf = EspressoNscf(calc.directory, prefix = calc.prefix,
                    occupations = 'tetrahedra',
                    kpts = (2, 2, 2),
                    debug = True,
                    )
    nscf.run()

List of all Methods
--------------------

.. autoclass:: xespresso.xespresso.Espresso
   :members:
