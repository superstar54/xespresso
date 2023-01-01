from ase.io import read
from ase.build import bulk
from ase.visualize import view
from pseudo import mypseudo
from xespresso import Espresso
from oer import OER_bulk
from copy import deepcopy

#
queue = {"nodes": 2, "ntasks-per-node": 20, "partition": "debug", "time": "00:10:00"}
calc = {
    "pseudopotentials": mypseudo,
    "calculation": "relax",
    "ecutwfc": 40.0,
    "ecutrho": 320.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "kpts": (4, 4, 1),
    "queue": queue,
    "debug": True,
    "parallel": "-nk 5",
}
atoms = bulk("Pt")
oer = OER_bulk(
    atoms,
    label="Pt",
    indices=[(0, 0, 1), (1, 1, 0), (1, 1, 1)],
    nlayer=1,
    fix=[0, 2],
    calculator=calc,
    view=True,
)
# oer.build_surface(indices = [(0, 0, 1), (0, 1, 0), (0, 0 ,1)])
# oer.build_surface()
# oer.view()
oer.run()
