from ase.io import vasp
from xespresso.uscsd import Uscsd

atoms = vasp.read_vasp("./STO.vasp")
pseudopotentials = {
    "Sr": "Sr.pbesol-spn-rrkjus_psl.1.0.0.UPF",
    "Ti": "Ti.pbesol-spn-rrkjus_psl.1.0.0.UPF",
    "O": "O.pbesol-n-rrkjus_psl.1.0.0.UPF ",
}

parameters_hp = {"nq1": 1, "nq2": 1, "nq3": 1, "find_atpert": 1, "conv_thr_chi": 0.001}

input_ntyp = {
    # 'starting_magnetization': {'Ti': 1},
    "Hubbard_U": {"Ti": 3.0},
}
queue = {
    "nodes": 1,
    "ntasks-per-node": 12,
    "partition": "debug",
    "time": "00:10:00",
    "account": "mr26",
    "constraint": "mc",
    "cpus-per-task": 1,
    "ntasks-per-core": 1,
}

input_data = {
    "calculation": "scf",
    "ecutwfc": 40,
    "ecutrho": 320,
    "occupations": "smearing",
    "lda_plus_u": True,
    "degauss": 0.01,
    "verbosity": "high",
    "input_ntyp": input_ntyp,
    # 'nspin' :2,
}
test = Uscsd(
    label="Sto_test/xespresso",
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(3, 3, 3),
    input_data=input_data,
    parameters_hp=parameters_hp,
    Hubbard_site="Ti",
    max_iter=6,
    sc_type="relax",
    relax_type="vc-relax",
    ethr_relax=1.0e-8,
    ethr_scf=1.0e-8,
    ethr_scf2=1.0e-8,
    parallel="-nk 1",
    queue=queue,
    queue_hp=queue,
    parallel_hp="-nk 1",
)
test.run_sc()
print("U values:", test.UVscsd)
print("Energy:", test.etot)
print("Structure:", test.opt_structures)
