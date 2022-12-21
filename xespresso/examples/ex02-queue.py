"""
edit file ~/.espressorc, and add lines need for HPC, for example:

"
module load QuantumESPRESSO
"

"""
from ase.build import bulk
from xespresso import Espresso

# build Fe
atoms = bulk("Fe")
# input parameters for pw
input_data = {
    "verbosity": "high",
    "ecutwfc": 40.0,
    "ecutrho": 302.0,
    "occupations": "smearing",
    "smearing": "gaussian",
    "degauss": 0.03,
    #
    "mixing_beta": 0.3,
    "conv_thr": 1.0e-8,
}
pseudopotentials = {
    "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
}
queue = {
    "nodes": 1,
    "ntasks-per-node": 4,
    "account": "dcb",
    "partition": "all",
    "time": "0:10:00",
}
calc = Espresso(
    pseudopotentials=pseudopotentials,
    queue=queue,
    label="scf/fe",  # 'scf/fe' is the directory, 'fe' is the prefix
    input_data=input_data,
    kpts=(6, 6, 6),
)
atoms.calc = calc
e = atoms.get_potential_energy()
print("Energy: {0:1.3f}".format(e))

"""
Energy: -3368.435
"""
