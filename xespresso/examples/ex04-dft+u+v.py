"""
Example using 1) AFM 2) DFT+U+V
"""
from ase.build import bulk
from ase.io import read
from xespresso import Espresso
import numpy as np

atoms = bulk("Fe", cubic=True)
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
atoms.arrays["species"][1] = "Fe1"
input_ntyp = {
    "starting_magnetization": {
        "Fe": 0.5,
        "Fe1": -0.5,
    },
    "Hubbard_U": {"Fe": 4.3, "Fe1": 4.3},
}
input_data = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "input_ntyp": input_ntyp,
    "hubbard_v": {"(1,1,1)": 4.0, "(3,3,1)": 1.0},  # Hubbard_V
}
pseudopotentials = {
    "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Fe1": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
}
queue = {"nodes": 1, "ntasks-per-node": 4, "partition": "debug", "time": "0:10:00"}
calc = Espresso(
    pseudopotentials=pseudopotentials,
    label="scf/fe+u+v",
    input_data=input_data,
    kpts=(4, 4, 4),
    queue=queue,
)
atoms.calc = calc
e = atoms.get_potential_energy()
print("Energy {0:1.3f}".format(e))

"""
Energy -10451.104
"""
