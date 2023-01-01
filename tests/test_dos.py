from _common_helpers import set_envs
import numpy as np


def test_dos():
    from ase.build import molecule
    from xespresso import Espresso
    from xespresso.post.dos import EspressoDos

    set_envs()
    atoms = molecule("CO")
    atoms.center(3)
    atoms.pbc = [True, True, True]
    pseudopotentials = {
        "O": "O.pbe-n-rrkjus_psl.1.0.0.UPF",
        "C": "C.pbe-n-rrkjus_psl.1.0.0.UPF",
    }
    calc = Espresso(
        label="calculations/scf/co",
        pseudopotentials=pseudopotentials,
        ecutwfc=30,
        occupations="smearing",
        degauss=0.03,
        kpts=(1, 1, 1),
        debug=True,
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print("Energy = {0:1.3f} eV".format(e))
    assert np.isclose(e, -606.94121029)
    # ===============================================================
    # start nscf calculation
    fe = calc.get_fermi_level()
    print("fermi level: ", fe)
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
    # ===============================================================
    dos = EspressoDos(
        parent_directory="calculations/scf/co",
        prefix=calc.prefix,
        Emin=fe - 30,
        Emax=fe + 30,
        DeltaE=0.01,
    )
    dos.run()


def test_dos_analysis():
    from xespresso.dos import DOS
    import matplotlib.pyplot as plt

    # DOS analysis
    dos = DOS(label="calculations/scf/co", prefix="co")
    dos.read_dos()
    dos.plot_dos(Emin=-10, Emax=10, smearing=[0.02, 0.01])
    plt.savefig("images/co-dos.png")
