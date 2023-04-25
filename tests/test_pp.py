from _common_helpers import set_envs
import numpy as np


def test_pp():
    from ase.build import fcc111
    from xespresso import Espresso
    from xespresso.post.pp import EspressoPp

    set_envs()
    atoms = fcc111("Al", (1, 1, 4), vacuum=5)
    atoms.translate([0, 0, 0.1 - min(atoms.positions[:, 2])])
    pseudopotentials = {
        "Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
    }
    calc = Espresso(
        label="calculations/scf/al-fcc111",
        pseudopotentials=pseudopotentials,
        ecutwfc=30,
        occupations="smearing",
        degauss=0.03,
        kpts=(8, 8, 1),
        debug=True,
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print("Energy = {0:1.3f} eV".format(e))
    # assert np.isclose(e, -606.94121029)
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
        kpts=(8, 8, 1),
        debug=True,
    )
    nscf.run()
    # ===============================================================
    pp = EspressoPp(
        parent_directory=calc.directory,
        prefix=calc.prefix,
        plot_num=11,
        fileout="potential.cube",
        iflag=3,
        output_format=6,
        debug=True,
    )
    pp.run()
    wf = calc.get_work_function()
    print(f"work function: {wf} eV")


# def test_pp_analysis():
#     from xespresso.dos import DOS
#     import matplotlib.pyplot as plt

#     # DOS analysis
#     dos = DOS(label="calculations/scf/co", prefix="co")
#     dos.read_dos()
#     dos.plot_dos(Emin=-10, Emax=10, smearing=[0.02, 0.01])
#     plt.savefig("images/co-dos.png")
