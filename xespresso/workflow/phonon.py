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

from phonopy import Phonopy
import numpy as np

import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
import subprocess


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
        self.supercells = {}
        self.set_of_forces = []
    def run(self):
        pass
    def phonopy(self, 
                supercell_matrix = [[1,0,0], [0,1,0], [0,0,1]], 
                primitive_matrix = [[1,0,0], [0,1,0], [0,0,1]],
                distance = 0.01):
        '''
        '''
        self.build_phonon(supercell_matrix = supercell_matrix, 
                         primitive_matrix = primitive_matrix,
                         distance = distance)
        self.build_supercells()
        self.pool_atoms(self.supercells)
        self.post()
    def build_phonon(self, atoms = None, 
                       supercell_matrix = None, 
                       primitive_matrix = None,
                       distance = None):
        if atoms is None:
            atoms = self.atoms
        phonon = Phonopy(atoms,
                        supercell_matrix,
                        primitive_matrix = primitive_matrix)
        phonon.generate_displacements(distance = distance)
        self.log("[Phonopy] Atomic displacements:")
        disps = phonon.get_displacements()
        for d in disps:
            self.log("[Phonopy] %d %s" % (d[0], d[1:]))
        self.phonon = phonon
    def build_supercells(self):
        supercells = self.phonon.get_supercells_with_displacements()
        # Force calculations by calculator
        i = 0
        for scell in supercells:
            atoms = Atoms(symbols=scell.get_chemical_symbols(),
                        scaled_positions=scell.get_scaled_positions(),
                        cell=scell.get_cell(),
                        pbc=True)
            self.supercells[str(i).zfill(3)] = atoms
            i += 1
        self.nsupercells = len(self.supercells)
        self.log('')
        self.log('Number of supercells: %s'%(i))
    def post(self):
        # read forces
        for i in range(self.nsupercells):
            job = str(i).zfill(3)
            forces = self.results[job]['forces']
            drift_force = forces.sum(axis=0)
            self.log(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force))
            # Simple translational invariance
            for force in forces:
                force -= drift_force / forces.shape[0]
            self.set_of_forces.append(forces)
        #
        self.phonon.produce_force_constants(forces=self.set_of_forces)
        self.log('')
        self.log("[Phonopy] Phonon frequencies at Gamma:")
        for i, freq in enumerate(self.phonon.get_frequencies((0, 0, 0))):
            self.log("[Phonopy] %3d: %10.5f THz" %  (i + 1, freq)) # THz

        # DOS
        self.phonon.set_mesh([41, 41, 41])
        self.phonon.set_total_DOS(tetrahedron_method=True)
        self.log('')
        self.log("[Phonopy] Phonon DOS:")
        self.omega, self.phdos = self.phonon.get_total_DOS()
        self.phonon_energies = 0.00414*self.omega
        for omega, dos in np.array(self.phonon.get_total_DOS()).T:
            self.log("%15.7f%15.7f" % (omega, dos))
    def opt(self, **kwargs):
        # vc-relax
        self.log('')
        self.log('Optmizatiion:')
        self.old_calculator = self.calculator.copy()
        self.calculator.update({'calculation': 'vc-relax'})
        self.calculator.update(kwargs)
        self.pool_atoms({'vc-relax': self.atoms})
        self.atoms = self.results['vc-relax']['atoms']
        self.calculator = self.old_calculator
        # scf
        # self.calculator.update({'calculation': 'scf'})
        # self.pool_atoms({'scf': self.atoms})
    def run_atoms(self, job, atoms):
        self.log('-'*60)
        self.log('Submit job {0}'.format(job))
        self.log.print_atoms(atoms)
        calc = Espresso(
                label = os.path.join(self.label, job),
                **self.calculator,
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
    
    def dfpt(self, job, queue):
        '''
        '''
        atoms = self.atoms
        self.log('-'*60)
        self.log('Submit job {0}'.format(job))
        self.log.print_atoms(atoms)
        calc = Espresso(
                label = os.path.join(self.label, job),
                **self.calculator,
                )
        calc.parameters['input_data']['prefix'] = job
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        # dos
        calc.post(queue = queue,
            package = 'ph',
            tr2_ph = 1e-14,
            ldisp = True,
            nq1 = 4,
            nq2 = 4, 
            nq3 = 3,
            )
        calc.post(queue = queue,
            fildyn = 'matdyn',
            package = 'q2r',
            zasr = 'simple',
            flfrc = '%s.fc'%job,
            )
        calc.post(queue = queue,
            package = 'matdyn',
            asr = 'simple',
            dos = True,
            flfrc = '%s.fc'%job,
            fldos = '%s.phdos'%job,
            nk1 = 40, 
            nk2 = 40,
            nk3 = 30,
            )
        calc.plot_phdos(fldos = '%s'%job,
                output = 'images/%s-phdos.png'%job)
        pass
    
    def read_phdos(self, fldos = None, ax = None, output = None):
        '''
        '''
        import matplotlib.pyplot as plt
        if fldos is None:
            fldos = self.prefix
        phdos = np.loadtxt(self.directory+'/%s.phdos' % fldos)
        self.phdos = phdos
    def plot_phdos(self, Emin = -5, Emax = 20, ax = None, output = None):
        '''
        '''
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(figsize = (6, 3))
            # ax = plt.gca()
        xindex = (self.omega > Emin) & (self.omega < Emax)
        ax.plot(self.omega[xindex], self.phdos[xindex], linewidth=0.7)
        ax.set_xlabel('Frequency (THz)')
        ax.set_ylabel('Phonon DOS (a.u.)')
        plt.tight_layout()
        if output:
            plt.savefig('%s' %output)
        return ax