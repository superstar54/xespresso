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
from xespresso.xlog import XLogger

import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor

class Base():
    def __init__(self, atoms, label = '.', prefix = None, calculator = None, view = False):
        self.set_label(label, prefix)
        self.atoms = atoms
        self.calculator = calculator
        self.view = view
        self.images = {}
        self.children = {}
        self.results = {}
        # self('%s'%self.__class__.__name__)
    def set_logger(self, logger):
        self.log = logger()
        self.log.fd = os.path.join(self.label, '%s.txt'%self.prefix)
        self.log('directory: %s'%self.directory)
        self.log('prefix: %s'%self.prefix)
        self.log('-'*60)
        self.log('Atomic structure:')
        self.log.print_atoms(self.atoms)
        # self.log('Calculator parameters:')
        # self.log.print_dict(self.calculator)
    def set_label(self, label, prefix):
        self.label = label
        if not prefix:
            self.directory = os.path.split(label)[0]
            self.prefix = os.path.split(label)[1]
        else:
            self.directory = label
            self.prefix = prefix
        if not os.path.exists(label):
            os.makedirs(label)
    def view_images(self, ):
        view(self.images.values())
    def run(self):
        self.build_children()
        self.run_children()
    def run_children(self):
        from random import random
        max_workers = len(self.children)
        with ThreadPoolExecutor(max_workers = max_workers) as executor:
            futures = []
            for job, child in self.children.items():
                sleep(random())
                futures.append(executor.submit(child.run))
            for future in as_completed(futures):
                print(future.result())
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
        return job, self.results[job]['energy']
    def pool_atoms(self, jobs):
        from random import random
        max_workers = len(jobs)
        with ThreadPoolExecutor(max_workers = max_workers) as executor:
            futures = []
            for job, atoms in jobs.items():
                sleep(random())
                futures.append(executor.submit(self.run_atoms, job, atoms))
            for future in as_completed(futures):
                print(future.result())