from ase.io import read
from ase.build import bulk, surface, fcc111
from ase.visualize import view
from pseudo import mypseudo
from xespresso import Espresso
from oer import OER_bulk, OER_surface, fix_layers, OER_pourbaix
from copy import deepcopy
#
queue = {'nodes': 2, 'ntasks-per-node': 20, 'partition': 'debug', 'time': '00:10:00'}
calc = {
'pseudopotentials': mypseudo,
'calculation':'relax',
'ecutwfc': 40.0,
'ecutrho': 320.0,
'occupations': 'smearing',
'degauss': 0.02,
'kpts': (4, 4, 1),
'queue': queue,
'debug': True,
'parallel': '-nk 5',
}

atoms = bulk('Pt')
atoms = surface(atoms, (0, 0, 1), 3, vacuum=1)
atoms.pbc = [True, True, True]
atoms = atoms*[2, 2, 1]
atoms.positions[:, 2] -= min(atoms.positions[:, 2]) - 0.01
atoms.cell[2][2] = max(atoms.positions[:, 2]) + 10
atoms = fix_layers(atoms, (0, 0, 1), 1.0, [0, 2])
# sites_dict = {'fcc': [2, 2, 2.2], 'ontop': [0, 0, 2.2]}
# atoms.info['sites'] = sites_dict
oer = OER_pourbaix(atoms, label = 'Pt/001-221-Pt', prefix = '001-221-Pt', adsorbates = [], calculator = calc)
oer.run()
