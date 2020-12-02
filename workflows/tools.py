'''
0 bulk
1 clean surfaces, different termination
2 add OH* and O*, different coverage
3 build Pourbaix diagram
4 choose surface 
5 add O*, OH* and OO*
6 free energy diagram
7 report

'''

import os
from ase.geometry import get_layers
from ase.constraints import FixAtoms
from ase.io.espresso import construct_namelist
from xespresso import Espresso
from xespresso.xio import build_atomic_species_str
from ase.dft.bandgap import bandgap
import pickle
import multiprocessing
import numpy as np



class OER():
    def __init__():
        pass
    def build_surface():
        '''
        surfaces and terminations
        '''
        pass
    def pourbaix_diagram():
        '''
        '''
        pass
    def free_energy_diagram():
        '''
        '''
        pass
    def coverage():
        '''
        '''
        pass
    def add_adsorbate():
    '''
    '''
    #----------------------------------
    ooh = Atoms('O2H', positions = [[0, 0, 0], [1.4, 0, 0], [1.4, 0, 1.0]])
    ooh.rotate('z', np.pi/4)
    mols = {
    'o' : Atoms('O'),
    'oh' : Atoms('OH', positions = [[0, 0, 0], [0, 0, 1.0]]),
    'ooh' : ooh,
    }
    #
    jobs = {}
    maxz = max(atoms.positions[:, 2])
    indm = [atom.index for atom in atoms if atom.z > maxz - 1.0 and atom.symbol != 'O'][0]
    for job, mol in mols.items():
        # print(job, mol)
        ads = mol.copy()
        natoms = atoms.copy()
        ads.translate(atoms[indm].position - ads[0].position + [0, 0, 1.9])
        natoms = natoms + ads
        jobs[job] = natoms
    return jobs

def fix_layers(atoms, miller = (0, 0, 1), tol = 1.0, n = [0, 4]):
    '''
    '''
    layers = get_layers(atoms, miller, tol)[0]
    index = [j for j in range(len(atoms)) if layers[j] in range(n[0], n[1])]
    constraint = FixAtoms(indices=index)
    atoms.set_constraint(constraint)
    return atoms

