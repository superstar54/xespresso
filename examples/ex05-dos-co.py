from xespresso.post.projwfc import EspressoProjwfc
from xespresso.post.dos import EspressoDos
from xespresso.post.nscf import EspressoNscf
from ase.build import molecule
from xespresso import Espresso
from xespresso.dos import DOS
import matplotlib.pyplot as plt

atoms = molecule("CO")
atoms.center(5)
atoms.pbc = [True, True, True]
pseudopotentials = {
    "O": "O.pbe-n-rrkjus_psl.1.0.0.UPF",
    "C": "C.pbe-n-rrkjus_psl.1.0.0.UPF",
}
queue = None

calc = Espresso(
    label="scf/co",
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
fe = calc.get_fermi_level()
print("fermi level: ", fe)
# ===============================================================
# start nscf calculation
nscf = EspressoNscf(
    calc.directory,
    prefix=calc.prefix,
    occupations="tetrahedra",
    kpts=(2, 2, 2),
    debug=True,
)
nscf.run()
# ===============================================================
# dos
dos = EspressoDos(
    parent_directory="scf/co",
    prefix=calc.prefix,
    Emin=fe - 30,
    Emax=fe + 30,
    DeltaE=0.01,
)
dos.run()
#
# DOS analysis
dos = DOS(label="scf/co", prefix="co")
dos.read_dos()
dos.plot_dos(Emin=-10, Emax=10, smearing=[0.02, 0.01])
plt.savefig("images/co-dos.png")

# ===========================================
# pdos
projwfc = EspressoProjwfc(parent_directory="scf/co", prefix="co", DeltaE=0.01)
projwfc.run()
# PDOS analysis
dos = DOS(label="scf/co", prefix="co")
dos.read_pdos()
dos.plot_pdos(Emin=-10, Emax=10, smearing=[0.02, 0.01], legend=True)
plt.savefig("images/co-pdos.png")
