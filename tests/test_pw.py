from _common_helpers import set_envs
import numpy as np


def test_scf():
    from ase.build import molecule
    from xespresso import Espresso

    set_envs()
    h2 = molecule("H2")
    h2.cell = [8, 8, 8]
    h2.pbc = [True, True, True]
    pseudopotentials = {
        "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
    }
    calc = Espresso(
        pseudopotentials=pseudopotentials,
        label="calculations/scf/h2",  # 'scf/fe' is the directory, 'fe' is the prefix
        ecutwfc=20,
        occupations="smearing",
        degauss=0.03,
        kpts=(1, 1, 1),
        debug=True,
    )
    h2.calc = calc
    e = h2.get_potential_energy()
    print("Energy: {0:1.4f}".format(e))
    assert np.isclose(e, -31.4454)


def test_relax():
    from ase.build import molecule
    from xespresso import Espresso

    set_envs()
    atoms = molecule("H2")
    atoms.center(5)
    atoms.pbc = [True, True, True]
    pseudopotentials = {"H": "H.pbe-rrkjus_psl.1.0.0.UPF"}
    calc = Espresso(
        label="calculations/relax/h2",
        pseudopotentials=pseudopotentials,
        calculation="relax",
        ecutwfc=20,
        kpts=(1, 1, 1),
        debug=True,
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print("Energy = {0:1.4f} eV".format(e))
    assert np.isclose(e, -31.4552)
