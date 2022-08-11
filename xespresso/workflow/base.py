import os
import numpy as np
from copy import deepcopy
from time import sleep
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
import multiprocessing
from xespresso import Espresso
from xespresso.tools import mypool, fix_layers, dipole_correction


class Base():
    def __init__(self, atoms, label='.', prefix=None, calculator=None, view=False, debug=False):
        self.set_label(label, prefix)
        self.atoms = atoms
        self.calculator = calculator
        self.view = view
        self.debug = debug
        self.images = {}
        self.children = {}
        self.results = {}
        # self('%s'%self.__class__.__name__)

    def set_logger(self, logger):
        self.log = logger()
        # self.log.fd = os.path.join(self.label, '%s.txt'%(self.prefix))
        self.log.fd = os.path.join(self.label, '%s.%so' %
                                   (self.prefix, self.name))
        self.log('-'*60)
        self.log('Class:  %s\n' % self.__class__.__name__)
        self.log('directory: %s' % self.directory)
        self.log('prefix: %s' % self.prefix)
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

    def view_images(self, images=None):
        """
        """
        from ase.visualize import view
        if not images:
            images = self.images
        print('%s, Number of images: %s' % (self.prefix, len(images)))
        if len(images) == 0:
            return
        view(images.values())

    def run(self):
        self.build_children()
        self.run_children()

    def run_children(self):
        from random import random
        max_workers = len(self.children)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for job, child in self.children.items():
                sleep(random())
                futures.append(executor.submit(child.run))
            for future in as_completed(futures):
                print(future.result())

    def run_children_multiprocessing(self, showInfo=False):
        '''
        '''
        from random import random
        from time import sleep
        pool = multiprocessing.Pool(len(self.children))
        results = []
        for job, child in self.children.items():
            sleep(random())
            results.append(pool.apply_async(child.run))
        for r in results:
            r.get()
        pool.close()
        pool.join()

    def cell_relax_espresso(self, job, atoms):
        """
        """
        self.log('-'*60)
        self.log('Submit job {0}'.format(job))
        self.log.print_atoms(atoms)
        calculator = deepcopy(self.calculator)
        calculator.update({'calculation': 'vc-relax',
                           'prefix': job, })
        calc = Espresso(
            label=os.path.join(self.label, job),
            **calculator,
        )
        atoms.calc = calc
        calc.run(atoms=atoms)
        calc.read_results()
        self.results[job] = deepcopy(calc.results)
        return job, self.results[job]['energy']

    def geo_pdos_zpe(self, job, atoms):
        self.geo_relax_espresso(job, atoms)
        atoms = self.results[job]['atoms']
        self.vib_zpe(job, atoms)

    def get_kpts(self, atoms):
        """
        Generate k-point by density or unit cell.
        """
        a, b, c, alpha, beta, gamma = atoms.get_cell_lengths_and_angles()
        kpts = (int(30/a), int(30/b), int(30/c))
        return kpts

    def geo_relax_espresso(self, job, atoms):
        """
        """
        self.log('-'*60)
        self.log('Run geometry relaxtion: {0}'.format(job))
        self.log.print_atoms(atoms)
        calculator = deepcopy(self.calculator)
        calculator.update({'calculation': 'relax',
                           'prefix': job, })
        dip = dipole_correction(atoms, edir=3)
        calculator.update(dip)
        calc = Espresso(
            label=os.path.join(self.label, job),
            **calculator,
        )
        atoms.calc = calc
        calc.run(atoms=atoms)
        calc.read_results()
        energy = calc.results['energy']
        self.results[job] = deepcopy(calc.results)
        self.calcs[job] = calc
        return job, energy

    def pdos(self, job, atoms):
        """
        calculate the pdos
        """
        self.log('-'*60)
        self.log('Run pdos calculation: {0}'.format(job))
        calc = self.calcs[job]
        fe = calc.get_fermi_level()
        atoms = calc.results['atoms']
        calculator = deepcopy(self.calculator)
        kpts = [calculator['kpts'][0]*2, calculator['kpts']
                [1]*2, calculator['kpts'][2]*2]
        if 'parallel' not in calculator:
            calculator['parallel'] = '-nk 1'
        calc.nscf(queue=calculator['queue'], parallel=calculator['parallel'],
                  occupations='tetrahedra', kpts=kpts)
        calc.nscf_calculate()
        calc.read_results()
        calc.post(queue=calculator['queue'], package='dos',
                  Emin=fe - 30, Emax=fe + 30, DeltaE=0.01)
        calc.post(queue=calculator['queue'], package='projwfc',
                  Emin=fe - 30, Emax=fe + 30, DeltaE=0.01)

    def pool_atoms(self, jobs, func):
        from random import random
        if self.debug:
            for job, atoms in jobs.items():
                func(job, atoms)
            return 0
        max_workers = len(jobs)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for job, atoms in jobs.items():
                sleep(random())
                futures.append(executor.submit(func, job, atoms))
            for future in as_completed(futures):
                print(future.result())

    def build_children(self):
        pass

    def run_command(self, command=None):
        import subprocess
        from ase.calculators.calculator import CalculationFailed
        if command is None:
            command = self.command
        print('Running %s' % command)
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err
        errorcode = proc.wait()
        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Command "{}" failed in '
                   '{} with error code {}'.format(command,
                                                  path, errorcode))
            raise CalculationFailed(msg)
        print('Done')
