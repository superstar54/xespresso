from ase.build import bulk, fcc111
from ase.visualize import view
from xespresso import Espresso
import matplotlib.pyplot as plt


atoms = fcc111("Al", size=(1, 1, 2), vacuum=4.0)
pseudopotentials = {
    "Al": "Al.pbe-n-rrkjus_psl.1.0.0.UPF",
}
calc = Espresso(
    pseudopotentials=pseudopotentials,
    label="scf/al",
    ecutwfc=40,
    occupations="smearing",
    degauss=0.01,
    kpts=(4, 4, 1),
    debug=True,
)
atoms.calc = calc
e = atoms.get_potential_energy()
calc.read_results()
e = calc.results["energy"]
fermi_energy = calc.get_fermi_level()
print("Energy: {0:1.3f}".format(e))
# ===============================================================================
bandpath = atoms.cell.bandpath()
calc = Espresso(
    calculation="bands",
    label="scf/al",
    pseudopotentials=pseudopotentials,
    ecutwfc=40,
    occupations="smearing",
    degauss=0.01,
    kpts=bandpath,
    debug=True,
)
atoms.calc = calc
calc.run(atoms)
# =================================================
bs = calc.band_structure()
bs.reference = fermi_energy
bs.plot()
