from ase.build import molecule
from xespresso import Espresso

h2 = molecule("H2")
h2.cell = [8, 8, 8]
h2.pbc = [True, True, True]
pseudopotentials = {
    "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
}
calc = Espresso(
    pseudopotentials=pseudopotentials,
    label="scf/h2",  # 'scf/fe' is the directory, 'fe' is the prefix
    ecutwfc=40,
    occupations="smearing",
    degauss=0.03,
    kpts=(1, 1, 1),
    debug=True,
)
h2.calc = calc
e = h2.get_potential_energy()
print("Energy: {0:1.3f}".format(e))

"""
Energy: -31.730
"""
