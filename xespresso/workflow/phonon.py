import os
import numpy as np
from copy import deepcopy
from time import sleep

from ase.atoms import Atoms
from ase.build import surface
from ase.geometry import get_layers
from ase.constraints import FixAtoms
from ase.formula import Formula
from ase.visualize import view

from xespresso import Espresso
from xespresso.tools import mypool, fix_layers, dipole_correction
from xespresso.workflow.base import Base
from xespresso.xlog import XLogger

import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor


class PLogger(XLogger):
    """Class for handling all text output."""
    def __init__(self, ):
        XLogger.__init__(self,)
    def logo(self):
        self()
        self(' //===\\  ||   ||   //==\\  ')
        self(' ||   ||  ||   ||  ||    || ')
        self(' ||===//  ||===||  ||    || ')
        self(' ||       ||   ||  ||    ||  ')
        self(' ||       ||   ||   \\==//   ')
        self()


#view(mols.values())
class Phonon(Base):
    def __init__(self, atoms, label = '.', prefix = None, calculator = None, view = False):
        Base.__init__(self, atoms, label = label, prefix = prefix ,calculator=calculator, view = view)
        self.set_logger(PLogger)
    def run(self, ):
        # vc-relax
        self.calculator.update({'calculation': 'vc-relax'})
        self.pool_atoms({'vc-relax': self.atoms})
        # scf
        self.atoms = self.results['vc-relax']['atoms']
        self.calculator.update({'calculation': 'scf'})
        self.pool_atoms({'scf': self.atoms})
        # gamma
        self.results['scf']
    def run_atoms(self, job, atoms):
        self.log('-'*60)
        self.log('Submit job {0}'.format(job))
        self.log.print_atoms(atoms)
        dip = dipole_correction(atoms, edir = 3)
        calc = Espresso(
                label = os.path.join(self.label, job),
                **self.calculator,
                **dip,
                )
        calc.parameters['input_data']['prefix'] = job
        # self.log('    Commnad: {0}'.format(calc.command))
        atoms.calc = calc
        # atoms.get_potential_energy()
        calc.run(atoms = atoms, restart = 1)
        calc.read_results()
        self.results[job] = deepcopy(calc.results)
        calc.clean()
        # scf
        
        # matdyn
        return job, self.results[job]['energy']
