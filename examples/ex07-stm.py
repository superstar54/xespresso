import pickle
from ase.io.cube import read_cube_data
from ase.dft.stm import STM
from ase.build import bulk, fcc111
from ase.visualize import view
from xespresso import Espresso
from xespresso.dos import DOS
from xespresso.tools import get_nbnd
import matplotlib.pyplot as plt


atoms = fcc111("Al", size=(1, 1, 2), vacuum=4.0)
view(atoms)
print(atoms)
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
)
atoms.calc = calc
e = atoms.get_potential_energy()
calc.read_results()
e = calc.results["energy"]
fermi = calc.get_fermi_level()
print("Energy: {0:1.3f}".format(e))
# ===============================================================
# post calculation
calc.post(
    package="pp",
    plot_num=5,
    sample_bias=0.0735,
    iflag=3,
    output_format=6,
    fileout="stm1_1eV.cube",
)
# ========================
# Creates: 2d.png, 2d_I.png, line.png, dIdV.png
ldos, atoms = read_cube_data("scf/al/stm1_1eV.cube")
with open("scf/al/datas.pickle", "wb") as f:
    pickle.dump([ldos, 1.0, atoms.cell], f)
stm = STM("scf/al/datas.pickle")
z = 8.0
bias = 1.0
c = stm.get_averaged_current(bias, z)
x, y, h = stm.scan(bias, c, repeat=(3, 5))
# print(x, y, h)
plt.gca(aspect="equal")
plt.contourf(x, y, h, 40)
plt.colorbar()
plt.savefig("2d.png")
