from _common_helpers import set_envs
import numpy as np


def test_nscf():
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
    # ===============================================================
    # start nscf calculation
    from xespresso.post.nscf import EspressoNscf

    nscf = EspressoNscf(
        calc.directory,
        prefix=calc.prefix,
        occupations="tetrahedra",
        kpts=(2, 2, 2),
        debug=True,
    )
    nscf.run()
