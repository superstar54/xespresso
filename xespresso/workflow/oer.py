'''
0 bulk
1 clean surfaces, different termination
2 add OH* and O*, different coverage
3 build Pourbaix diagram
4 choose surface 
5 different OER site
6 add O*, OH* and OOH*
7 free energy diagram
8 report


'''
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

#----------------------------------
mols = {
'O' : Atoms('O'),
'OH' : Atoms('OH', positions = [[0, 0, 0], [0.6, 0.7, 0.15]]),
'OOH' : Atoms('OOH', positions = [[0, 0, 0], [0.8, 0.9, 0.0], [1.4, 1.6, 0.45]]),
}

#view(mols.values())
class OER_base(Base):
    def __init__(self, atoms, label = '.', prefix = None, calculator = None, view = False):
        Base.__init__(self, atoms, label = label, prefix = prefix ,calculator=calculator, view = view)
        self.set_logger(OERLogger)
    
class OER_bulk(OER_base):
    def __init__(self, atoms, label = '.', prefix = None, surfaces = {}, indices = [], nlayer = 4, fix = [0, 3], tol = 1.0, termination = None, calculator = None, view = False):
        '''
        Only one calculation, vc-relax for the bulk
        '''
        OER_base.__init__(self, atoms, label = label, prefix = prefix ,calculator=calculator, view = view)
        self.children = surfaces
        self.indices = indices
        self.nlayer = nlayer
        self.fix = fix
        self.tol = tol
    def build_children(self, indices = None, nlayer = None, fix = None, tol = None):
        '''
        surfaces and terminations

        '''
        if not indices:
            indices = self.indices
        if not nlayer:
            nlayer = self.nlayer
        if not fix:
            fix = self.fix
        if not tol:
            tol = self.tol
        self.log('-'*60)
        self.log('Build surface:\n')
        count = 0
        for indice in indices:
            count += 1
            self.log('{0:15s}: {1}'.format('%s. Indice'%count, indice))
            atoms = surface(self.atoms, indice, nlayer + 2, vacuum=1)
            atoms.pbc = [True, True, True]
            layers = get_layers(atoms, (0, 0, 1), tol)[0]
            terminations = self.get_terminations(atoms, layers)
            for ter, ind in terminations.items():
                self.log('{0:15s}: {1}'.format('    termination', ter))
                surf = atoms.copy()
                index = [j for j in range(len(surf)) if layers[j] in range(ind, ind - nlayer, -1)]
                surf = surf[index]
                surf.positions[:, 2] -= min(surf.positions[:, 2]) - 0.01
                surf.cell[2][2] = max(surf.positions[:, 2]) + 15
                surf = fix_layers(surf, (0, 0, 1), tol, fix)
                #
                prefix = '%s%s%s-%s'%(indice[0], indice[1], indice[2], ter)
                self.images[prefix] = surf
                surf = OER_pourbaix(surf, label = os.path.join(self.label, prefix), prefix = prefix, calculator = self.calculator, view = self.view)
                self.children[prefix] = surf
            self.log()
        self.log('Total number of surfaces: {0}\n'.format(len(self.children)))
    def get_terminations(self, slab, layers):
        '''
        determine how many terminations for a surface index
        for a slab with SrO, TiO2, SrO, TiO2 ... ordering, there will be two terminations
        'SrO' and 'TiO2' 
        Search from the top to the bottom
        '''
        nlayer = max(layers) + 1
        terminations = {}
        for i in range(nlayer - 2, 0, -1):
            atoms = slab[layers == i]
            formula = atoms.get_chemical_formula(mode = 'hill')
            if formula not in terminations:
                terminations[formula] = i
        return terminations
    
class OER_pourbaix(OER_base):
    def __init__(self, atoms, label = '', prefix = None, sites_dict = None, calculator = None, coverages = [0, 1], adsorbates = ['O', 'OH'], view = False):
        '''
        Calculate surface Pourbaix diagram.
        coverages: list
            A list of float, e.g. [0, 0.5, 1]
        adsorbates: list
            The adsorbate used to build the coverage, support ['O', 'OH']
        sites_dict: dict
            In OER_pourbaix, site is position. e.g. {'fcc': [0.5, 0.5, 3]}. 
            If None, the ontop site will be found and used.
        '''
        OER_base.__init__(self, atoms, label = label, prefix = prefix, calculator = calculator, view = view)
        if not sites_dict:
            self.sites_dict = self.get_sites_dict(self.atoms)
        else:
            self.sites_dict = sites_dict
        self.coverages = coverages
        self.adsorbates = [''] + adsorbates
        self.images = {}
        self.results = {}
        self.children = {'0-O-0-OH': atoms}
    def run(self):
        # coverage
        self.build_coverage()
        if self.view:
            self.view_images()
            return
        #
        self.run_coverage()
        fig, ax = plt.subplots(1, 1)
        ax, children = self.plot_pourbaix_diagram(self.results, ax = ax)
        filename = os.path.join(self.label, '%s-pourbaix_diagram.png'%self.prefix)
        fig.savefig(filename)
        self.log('-'*60)
        self.log('Save pourbaix diagram to: {0}'.format(filename))
        self.log('Stable surfaces from pourbaix diagram:')
        for prefix in children:
            surf = self.images[prefix]
            self.images[prefix] = surf
            surf = OER_surface(surf, label = os.path.join(self.label, prefix), prefix = prefix, calculator = self.calculator)
            self.log('{0:15s}: {1}'.format('    Stable surface', prefix))
            self.children[prefix] = surf
        self.log()
        self.log('Total number of stable surfaces: {0}\n'.format(len(self.children)))
        # OER
        self.run_children()
        
    def run_coverage(self):
        self.pool_atoms(self.images)
        self.log('-'*60)
        self.log('Image energies (eV):')
        for job, result in self.results.items():
            self.log('{0:10s}  Energy: {1:1.3f}'.format(job, result['energy']))
    def build_coverage(self, coverages = None, sites_dict = None):
        '''
        add 1ML O* and OH*
        could be:
        1 ML O*
        1 ML OH*
        0.5 ML O* and 0.5 ML OH*
        '''
        import itertools
        if coverages is None:
            coverages = self.coverages
        if sites_dict is None:
            sites_dict = self.sites_dict
        self.sites = []
        site_mols = []
        self.log('-'*60)
        self.log('Build coverages: ')
        self.log('Surface sites: ')
        self.log('    {0:10s} {1:4s}'.format('Symbol', 'Position'))
        count = 0
        for symbol, position in sites_dict.items():
            # self.log('    %s %s'%(symbol, len(positions)))
            count += 1
            self.log('    {0:3d}   {1:10s} {2}'.format(count, symbol, position))
            self.sites.append(position)
            site_mols.append(self.adsorbates)
        self.nsite = len(self.sites)
        self.log('Number of sites: {0}'.format(self.nsite))
        combinations = list(itertools.product(*site_mols))
        self.log('Coverages: {0}'.format(self.coverages))
        for combination in combinations:
            noh = combination.count('OH')
            no = combination.count('O')
            coverage = (no + noh)/self.nsite
            flag = [abs(coverage - x) for x in coverages]
            if min(flag) > 0.1: continue
            formula = '%s-O-%s-OH'%(no, noh)
            if formula in self.images: continue
            ads = Atoms('')
            for i in range(len(combination)):
                if not combination[i]: continue
                pos = self.sites[i]
                mol = mols[combination[i]].copy()
                mol.translate(pos- mol[0].position + [0, 0, 1.5])
                ads = ads + mol
            natoms = self.atoms.copy()
            self.images[formula] = natoms + ads
            self.log('Structure: {0}'.format(formula))
        self.log('Total number of structures: {0}'.format(len(self.images)))
    
    
    def plot_pourbaix_diagram(self, results, ax = None, E_O = None, E_OH = None):
        '''
        mh = potential
        '''
        from qedatas import E_h2o, E_h2
        E_h2o = E_h2o + 0.56 - 0.67
        E_h2 = E_h2 + 0.27 - 0.41
        if not E_O:
            E_O = E_h2o - E_h2
        if not E_OH:
            E_OH = E_h2o - E_h2/2.0
        potential = np.linspace(0, 2, 20)
        if ax is None:
            ax = plt.gca()
        E0 = results['0-O-0-OH']['energy']
        for job, result in results.items():
            E = result['energy']
            no = int(job.split('-')[0])
            noh = int(job.split('-')[2])
            dE = E - E0 - noh*(E_OH + potential) - no*(E_O + 2*potential)
            ax.plot(potential, dE, label = job)
            self.results[job]['pourbaix_energy'] = dE
        ax.legend()
        ax.set_xlabel('Potential (V)')
        ax.set_ylabel('$\Delta$ G (eV)')
        StableSurfaces = []
        for i in range(20):
            Emin = 1000
            jobmin = ''
            for job, result in self.results.items():
                E = result['pourbaix_energy']
                if E[i] < Emin:
                    Emin = E[i]
                    jobmin = job
            StableSurfaces.append(jobmin)
        StableSurfaces = list(set(StableSurfaces))
        return ax, StableSurfaces
    def get_sites_dict(self, atoms, miller = (0, 0, 1), tol = 1.0):
        '''
        '''
        layers = get_layers(atoms, miller, tol)[0]
        last_layer = max(layers)
        index = [j for j in range(len(atoms)) if layers[j] ==  last_layer]
        sites_dict = {}
        count = 0
        for i in index:
            symbol = atoms[i].symbol
            if symbol not in ['O', 'N']:
                count += 1
                if symbol not in sites_dict:
                    sites_dict['%s-%s'%(count, symbol)] = atoms[i].position
        return sites_dict

class OER_surface(OER_base):
    def __init__(self, atoms, label = '', prefix = None, sites_dict = None, calculator = None):
        '''
        Calculate surface Pourbaix diagram.
        sites_dict: dict
            the active site, should be the atom in the structure. e.g. {'Pt': [2]}. 
            If None, the ontop atom will be found and used.
            This is different from the site in OER_pourbaix.
        '''
        OER_base.__init__(self, atoms, label = label, prefix = prefix, calculator = calculator)
        self.atoms = atoms
        if not sites_dict:
            self.sites_dict = self.get_oer_sites_dict(self.atoms)
        self.children = {}
    def run(self):
        self.build_children()
        self.run_children()
        return 0
    def build_children(self, ):
        '''
        '''
        self.log('-'*60)
        self.log('OER sites:')
        for prefix, sites in self.sites_dict.items():
            if len(sites) == 0: continue
            ind = sites[0]
            self.images[prefix] = self.atoms
            child = OER_site(self.atoms, label = os.path.join(self.label, prefix), site = ind, prefix = prefix, calculator = self.calculator)
            self.log('{0:15s}: {1}'.format('    sites', prefix))
            self.children[prefix] = child
            self.log()
        self.log('Total number of OER sites: {0}\n'.format(len(self.children)))
    def get_oer_sites_dict(self, atoms, miller = (0, 0, 1), tol = 1.0):
        '''
        metal: Ti, La -> OH -> O -> OOH -> O2
        O:          O -> OOH -> O2 -> OH 
        H:          OH -> O -> OOH -> O2
        '''
        layers = get_layers(atoms, miller, tol)[0]
        last_layer = max(layers)
        index = [j for j in range(len(atoms)) if layers[j] ==  last_layer]
        sites_dict = {'O': [], 'H': []}
        for i in index:
            symbol = atoms[i].symbol
            if symbol not in sites_dict:
                sites_dict[symbol] = [i]
            else:
                sites_dict[symbol].append(i)
        return sites_dict

class OER_site(OER_base):
    def __init__(self, atoms, label = '', prefix = None, site = None, calculator = None):
        '''
        site: int
            The atom index for the active site.
        '''
        OER_base.__init__(self, atoms, label = label, prefix = prefix, calculator = calculator)
        self.site = site
        self.site_symbol = atoms[self.site].symbol
        self.children = {}
    def run(self):
        #
        self.build_oer_adsorbate()
        #
        self.pool_atoms(self.children)
        self.log('-'*60)
        self.log('Image energies (eV):')
        self.log('{0:10s}  Energy: {1:1.3f}'.format(self.site_symbol, self.atoms.calc.results['energy']))
        for job, result in self.results.items():
            self.log('{0:10s}  Energy: {1:1.3f}'.format(job, result['energy']))
        #
        self.get_free_energy()
        self.log('-'*60)
        self.log('Gibbs free energies (eV): ')
        for e in self.free_energies:
            self.log('        {0:1.2f}'.format(e))
        self.log('{0:10s}: {1:1.2f} eV'.format('Over potential', self.over_potential))
        # view(self.children.values())
        return 'Over potential', self.over_potential
    def build_oer_adsorbate(self, ):
        '''
        '''
        self.log('-'*60)
        self.log('Adsoption site: index ({0})  symbol({1})'.format(self.site, self.site_symbol))
        self.log('OER adsorbate:')
        ind = self.site
        atoms = self.atoms
        if self.site_symbol not in ['H', 'O']:
            for prefix, mol in mols.items():
                # print(prefix, mol)
                self.log('    %s'%prefix)
                ads = mol.copy()
                natoms = atoms.copy()
                ads.translate(atoms[ind].position - ads[0].position + [0, 0, 1.9])
                natoms = natoms + ads
                self.children[prefix] = natoms
        # O
        if self.site_symbol == 'O':
            atoms = self.atoms.copy()
            del atoms[[ind]]
            self.log('    O2')
            self.children['O2'] = atoms.copy()
            for prefix, mol in mols.items():
                # print(prefix, mol)
                if prefix == 'O': continue
                self.log('    %s'%prefix)
                ads = mol.copy()
                natoms = atoms.copy()
                ads.translate(self.atoms[ind].position - ads[0].position)
                natoms = natoms + ads
                self.children[prefix] = natoms
        # OH
        if self.site_symbol == 'H':
            atoms = self.atoms.copy()
            dis = atoms.get_distance(ind, range(len(atoms)))
            indo = dis.index(min(dis))
            # O
            self.log('    O')
            natoms = atoms.copy()
            del natoms[[ind]]
            self.children['O'] = natoms
            # O2
            self.log('    O2')
            natoms = atoms.copy()
            del natoms[[ind, indo]]
            self.children['O2'] = natoms
            # OOH
            self.log('    OOH')
            ads = mols['OOH'].copy()
            natoms = atoms.copy()
            ads.translate(self.atoms[indo].position - ads[0].position)
            natoms = natoms + ads
            self.children['OOH'] = natoms
        self.log('Total number of OER adsorbate: {0}\n'.format(len(self.children)))
    def get_free_energy(self, results = None, E_O = None, E_OH = None, E_OOH = None, E_H = None):
        '''
        mh = potential
        '''
        from qedatas import E_h2o, E_h2, E_o2
        E_h2o = E_h2o + 0.56 - 0.67
        E_h2 = E_h2 + 0.27 - 0.41
        E_o2 = 0.10 - 0.64
        if not E_O:
            E_O = E_h2o - E_h2
        if not E_OH:
            E_OH = E_h2o - E_h2/2.0
        if not E_OOH:
            E_OOH = 2*E_h2o - 3*E_h2/2.0
        if not E_H:
            E_H = E_h2/2.0
        E0 = self.atoms.calc.results['energy']
        fig, ax = plt.subplots()
        if self.atoms[self.site].symbol not in ['H', 'O']:
            E_OH = self.results['OH']['energy'] - E0 - E_OH
            E_O = self.results['O']['energy'] - E0 - E_O
            E_OOH = self.results['OOH']['energy'] - E0 - E_OOH
            E_O2 = 4.92
            self.free_energies = [E_OH, E_O - E_OH, E_OOH - E_O, E_O2 - E_OOH]
            self.plot_free_energy([0, E_OH, E_O, E_OOH, E_O2], ['*', 'OH', 'O', 'OOH', 'O2'], ax = ax)
        if self.atoms[self.site].symbol == 'O':
            E_OOH = self.results['OOH']['energy'] - E0 - E_OH
            E_O2 = self.results['O2']['energy'] - E0 + E_O
            E_OH = self.results['OH']['energy'] - E0 - E_H
            E_O = 4.92
            self.free_energies = [E_OOH, E_O2 - E_OOH, E_OH - E_O2, E_O - E_OH]
            self.plot_free_energy([0, E_OOH, E_O2, E_OH, E_O], ['O', 'OOH', 'O2', 'OH', 'O'], ax = ax)
        if self.atoms[self.site].symbol == 'H':
            E_O = self.results['OH']['energy'] - E0 + E_H
            E_OOH = self.results['O']['energy'] - E0 - E_O
            E_O2 = self.results['O2']['energy'] - E0 + E_OH
            E_OH = 4.92
            self.free_energies = [E_O, E_OOH - E_O, E_O2 - E_OOH, E_OH - E_O2]
            self.plot_free_energy([0, E_O, E_OOH, E_O2, E_OH], ['OH', 'O', 'OOH', 'O2', 'OH'], ax = ax)
        fig.savefig(os.path.join(self.label, '%s-free-energy.png'%self.prefix))
        self.over_potential = max(self.free_energies) - 1.23
            
    def plot_free_energy(self, E, labels, ax = None):
        if ax is None:
            ax = plt.gca()
        ax.plot([0, 1], [0, 0], '-')
        ax.text(0 + 0.3, E[0] + 0.5, labels[0])
        for i in range(1, 5):
            ax.plot([i, i], [E[i - 1], E[i]], '-')
            ax.plot([i, i+1], [E[i], E[i]], '-')
            ax.text(i + 0.3, E[i] + 0.5, labels[i])
        ax.set_xlabel('Step')
        ax.set_ylabel('$\Delta$ G (eV)')
        return ax
        

class OERLogger(XLogger):
    """Class for handling all text output."""
    def __init__(self, ):
        XLogger.__init__(self,)
    def logo(self):
        self()
        self('  //==\\   //====  //====\\  ')
        self(' ||    ||  ||      ||    || ')
        self(' ||    ||  ======  |\====// ')
        self(' ||    ||  ||      || \\    ')
        self('  \\==//   \\====  ||  \\   ')
        self()

