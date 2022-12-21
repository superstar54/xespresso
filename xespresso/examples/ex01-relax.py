from ase.build import molecule
from xespresso import Espresso

atoms = molecule("H2")
atoms.center(5)
atoms.pbc = [True, True, True]
pseudopotentials = {"H": "H.pbe-rrkjus_psl.1.0.0.UPF"}
calc = Espresso(
    label="relax/h2",
    pseudopotentials=pseudopotentials,
    calculation="relax",
    ecutwfc=20,
    kpts=(1, 1, 1),
    debug=True,
)
atoms.calc = calc
e = atoms.get_potential_energy()
print("Energy = {0:1.3f} eV".format(e))
