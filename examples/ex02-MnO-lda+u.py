from ase.build import bulk, sort
from ase.io import read
from ase.visualize import view
from xespresso import Espresso
import numpy as np

atoms = read("datas/MnO.cif")
# view(atoms)
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
# set species Mn1, Mn2, ...
for i in range(4):
    atoms.arrays["species"][2 * i] = atoms.arrays["species"][2 * i] + "%s" % (i + 1)
    atoms.arrays["species"][2 * i + 1] = atoms.arrays["species"][2 * i + 1] + "%s" % (
        i + 1
    )
# here I add a new species O1, just for testing
atoms.arrays["species"][15] = "O1"
# here you can set any parameters which support "(i), i=1,ntyp" in qw,
# for example  'starting_magnetization', 'Hubbard_U', 'charge', ...
input_ntyp = {
    "starting_magnetization": {
        "Mn1": 1.0,
        "Mn2": -1.0,
        "Mn3": -1.0,
        "Mn4": 1.0,
        "O1": 1.0,
    },
    "Hubbard_U": {
        "Mn1": 5.7504,
        "Mn2": 5.7606,
        "Mn3": 5.7490,
        "Mn4": 5.7551,
        "O1": 3.0,
    },
    "starting_charge": {"Mn1": 1},
}

input_data = {
    "ecutwfc": 70.0,
    "ecutrho": 840.0,
    "occupations": "smearing",
    "degauss": 0.01,
    "nspin": 2,
    "lda_plus_u": True,
    "input_ntyp": input_ntyp,
}
queue = {
    "nodes": 1,
    "ntasks-per-node": 4,
    "account": "dcb",
    "partition": "all",
    "time": "0:10:00",
}
pseudopotentials = {
    "Mn1": "Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF",
    "Mn2": "Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF",
    "Mn3": "Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF",
    "Mn4": "Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF",
    "O": "O.pbesol-n-rrkjus_psl.1.0.0.UPF",
    "O1": "O.pbesol-n-rrkjus_psl.1.0.0.UPF",
    "Sr": "Sr.pbesol-spn-rrkjus_psl.1.0.0.UPF",
}
# print(queue)
calc = Espresso(
    pseudopotentials=pseudopotentials,
    package="pw",
    parallel="-nk 2 -nt 4 -nd 144",  # parallel parameters
    queue=queue,
    label="scf/mno",
    input_data=input_data,
    kpts=(6, 6, 6),
)
atoms.set_calculator(calc)
e = atoms.get_potential_energy()

print("Energy: ", e)
