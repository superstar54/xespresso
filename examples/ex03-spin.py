from ase.build import bulk
from xespresso import Espresso
import numpy as np

atoms = bulk("Fe", cubic=True)
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
atoms.arrays["species"][0] = "Fe"
atoms.arrays["species"][1] = "Fe1"
print(atoms.arrays["species"])
input_ntyp = {
    "starting_magnetization": {
        "Fe": 1.0,
        "Fe1": -1.0,
    }
}
pseudopotentials = {
    "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "Fe1": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
}
calc = Espresso(
    pseudopotentials=pseudopotentials,
    label="scf/fe-afm",
    ecutwfc=40,
    occupations="smearing",
    degauss=0.02,
    nspin=2,
    input_data={"input_ntyp": input_ntyp},
    kpts=(4, 4, 4),
    debug=True,
)
atoms.calc = calc
e = atoms.get_potential_energy()
print("Energy  {0:1.3f}".format(e))

"""
Energy  -6737.189
"""
