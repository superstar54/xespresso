from ase.build import bulk, surface, fcc111
from pseudo import mypseudo
from xespresso import Espresso
from xespresso.workflow.oer import (
    OER_bulk,
    OER_surface,
    fix_layers,
    OER_pourbaix,
    OER_site,
)
from copy import deepcopy
from qedatas.qe_mol import G_h2o, G_h2

# I imported the Gibbs free energies of H2O and H2 from 'qedatas` module,
# which is a module build by myself to save mine datas.
# Please provide your own enegies here.

molecule_energies = {
    "H2O": G_h2o,
    "H2": G_h2,
}
# Parameters for QE calculation
queue = {
    "nodes": 2,
    "ntasks-per-node": 20,
    "partition": "bdw",
    "mem-per-cpu": "5G",
    "time": "23:59:00",
    "config": ".xespresso-intel-2020b",
}
parameters = {
    "pseudopotentials": mypseudo,
    "calculation": "relax",
    "ecutwfc": 40.0,
    "ecutrho": 320.0,
    "occupations": "smearing",
    "degauss": 0.01,
    "kpts": (4, 4, 1),
    "queue": queue,
    "debug": True,
    "parallel": "-nk 4",
}

# Build the Pt (001) surface
atoms = bulk("Pt")
atoms = surface(atoms, (0, 0, 1), 2, vacuum=1)
atoms.pbc = [True, True, True]
atoms = atoms * [2, 2, 1]
atoms.positions[:, 2] -= min(atoms.positions[:, 2]) - 0.01
atoms.cell[2][2] = max(atoms.positions[:, 2]) + 10
atoms = fix_layers(atoms, (0, 0, 1), 1.0, [0, 2])
#
calc = Espresso(
    label="relax/001-221-Pt",
    **parameters,
)

atoms.calc = calc
calc.run(atoms=atoms)
calc.read_results()
atoms_opt = calc.results["atoms"]
del atoms_opt.info["species"]

# ==========================================================================
# parameter for QE
parameters = {
    "pseudopotentials": mypseudo,
    "calculation": "relax",
    "ecutwfc": 40.0,
    "ecutrho": 320.0,
    "occupations": "smearing",
    "degauss": 0.01,
    "kpts": (4, 4, 1),
    "queue": queue,
    "debug": True,
    "parallel": "-nk 5",
}
# build the OER calculator
site = -1
oer = OER_site(
    atoms_opt,
    label="oer/Pt-001-ontop",
    site_type="ontop",
    site=-1,
    height=2.0,
    calculator=parameters,
    molecule_energies=molecule_energies,
)
oer.run()
