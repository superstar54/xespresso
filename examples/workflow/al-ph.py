from ase.build import bulk
from pseudo import rrkjus_psl
from xespresso import Espresso
from xespresso.workflow.phonon import Phonon

#
queue = {"nodes": 2, "ntasks-per-node": 20, "partition": "debug", "time": "00:10:00"}
paras = {
    "pseudopotentials": rrkjus_psl,
    "calculation": "relax",
    "ecutwfc": 40.0,
    "ecutrho": 320.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "kpts": (8, 8, 8),
    "queue": queue,
    "debug": True,
    "parallel": "-nk 4",
}

atoms = bulk("Li")
calc = Espresso(label="vc-relax/li", **paras)
atoms.calc = calc
e = atoms.get_potential_energy()
print("energy: ", e)
#
paras.update({"calculation": "scf", "conv_thr": 1e-10})
atoms = calc.results["atoms"]
calc = Espresso(label="scf/li", **paras)
atoms.calc = calc
e = atoms.get_potential_energy()
print("energy: ", e)
# Phonon calculations at gamma
input = {
    "tr2_ph": 1.0e-14,
    "epsil": True,
}
#
# phonon calculation on a (444) uniform grid of q-points
input = {
    "tr2_ph": 1.0e-12,
    "ldisp": True,
    "nq1": 4,
    "nq2": 4,
    "nq3": 4,
}
calc.post(queue=queue, package="ph", parallel="-nk 4", **input)
#
input = {
    "asr": "simple",
    "dos": True,
    "nk1": 12,
    "nk2": 12,
    "nk3": 12,
}
calc.post(queue=queue, package="matdyn", **input)
# oer = Phonon(atoms,
#                 label = 'ZnO',
#                 calculator = calc,
#                 view = True,
#                 )
# oer.run()
