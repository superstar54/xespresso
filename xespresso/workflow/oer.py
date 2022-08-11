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
from ase.atoms import Atoms
from ase.io import write
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from xespresso import Espresso
from xespresso.tools import mypool, fix_layers, dipole_correction
from xespresso.workflow.base import Base
from xespresso.xlog import XLogger

import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor

# ----------------------------------


# view(mols.values())
class OER_base(Base):
    def __init__(self, atoms, name='oer', label='.', prefix=None, calculator=None, molecule_energies=None, view=False, debug=False):
        Base.__init__(self, atoms, label=label, prefix=prefix,
                      calculator=calculator, view=view, debug=debug)
        self.name = name
        self.set_logger(OERLogger)
        self.molecule_energies = molecule_energies
        self.mols = {
            'O': Atoms('O'),
            'OH': Atoms('OH', positions=[[0, 0, 0], [0.6, 0.7, 0.15]]),
            'OOH': Atoms('OOH', positions=[[0, 0, 0], [0.8, 0.9, 0.0], [1.4, 1.6, 0.45]]),
        }
        '''
        from qedatas import E_h2o, E_h2, G_o2
        E_h2o = E_h2o + 0.56 - 0.67
        E_h2 = E_h2 + 0.27 - 0.41
        G_o2 = 0.10 - 0.64
        '''

    def get_sites_dict(self, atoms, excludes=[]):
        '''
        get the surface site for adsorption
        1) for OER_pourbaix, ['O', 'N'] sites are generally excluded.
        2) for OER_site
                metal: Ti, La -> OH -> O -> OOH -> O2
                O:          O -> OOH -> O2 -> OH 
                H:          OH -> O -> OOH -> O2

        Return
        sites: dict
            e.g. {'Pt': [0.0, 0.0, 5.0]}
        '''
        from scipy.spatial.distance import squareform, pdist
        slabs = AseAtomsAdaptor.get_structure(atoms)
        asf_slabs = AdsorbateSiteFinder(slabs)
        ads_sites = asf_slabs.find_adsorption_sites(distance=0.0)
        sites = {}
        count = 0
        # view(atoms)
        for pos in ads_sites['ontop']:
            poss = atoms.get_positions()
            poss = np.append(poss, pos).reshape(-1, 3)
            ind = np.argsort(squareform(pdist(poss))[-1])[1] - 1
            # print(ind)
            symbol = atoms[ind].symbol
            if symbol in excludes:
                continue
            sites['site-%s-%s' % (count, symbol)] = pos + np.array([0, 0, 2.0])
        return sites

    def vib_zpe(self, job, atoms):
        from ase.vibrations import Vibrations
        self.log('-'*60)
        self.log('Run ZEP calculation: {0}'.format(job))
        label = os.path.join(self.label, job)
        label = os.path.join(label, 'vib')
        # copy file
        calculator = deepcopy(self.calculator)
        dip = dipole_correction(atoms, edir=3)
        calculator.update(dip)
        calculator.update({'calculation': 'scf',
                           'tstress': True,
                           'tprnfor': True,
                           'outdir': '../',
                           'prefix': '%s' % job,
                           'startingpot': 'file',
                           'startingwfc': 'file',
                           'etot_conv_thr': 1e-6,
                           'disk_io': 'none',
                           })
        # better to restart from previous geo_relax calculation
        calc = Espresso(label=label,
                        **calculator,
                        )
        atoms.calc = calc
        if job[-1] == 'O':
            indices = [-1]
        elif job[-2:] == 'OH' and job[-3:] != 'OOH':
            indices = [-1, -2]
        elif job[-3:] == 'OOH':
            indices = [-1, -2, -3]
        elif job[-2:] == 'O2':
            indices = [-1, -2]
        else:
            indices = []
            # print('%!!!')
            return job, 0
        vib = Vibrations(atoms, name=os.path.join(label, job), indices=indices)
        vib.run()
        vib_energies = np.real(vib.get_energies())
        zpe = 0.
        for energy in vib_energies:
            zpe += 0.5 * energy
        self.zpes[job] = zpe
        return job, zpe

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
                [1]*2, calculator['kpts'][2]]
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


class OER_bulk(OER_base):
    def __init__(self, atoms, label='.', prefix=None, surfaces={}, indexs=[], min_slab_size=6.0, min_vacuum_size=15.0,  nlayer=4, fix=[0, 3], tol=1.0, termination=None, pourbaix_input={}, surface_input={}, calculator=None, molecule_energies=None, view=False, debug=False):
        '''
        Generate surface slab model from bulk structure

        To do: adding vc-relax for the bulk
        '''
        OER_base.__init__(self, atoms, label=label, prefix=prefix,
                          calculator=calculator, view=view, molecule_energies=molecule_energies)
        self.name = 'bulk'
        self.children = {}
        self.surfaces = surfaces
        self.indexs = indexs
        self.min_slab_size = min_slab_size
        self.min_vacuum_size = min_vacuum_size
        self.nlayer = nlayer
        self.fix = fix
        self.tol = tol
        self.pourbaix_input = pourbaix_input
        self.surface_input = surface_input

    def get_slabs(self, atoms=None, index=None):
        """
        """
        from math import ceil
        if not atoms:
            atoms = self.atoms
        struct = AseAtomsAdaptor.get_structure(self.atoms)
        struct = SpacegroupAnalyzer(
            struct).get_conventional_standard_structure()
        slabgen = SlabGenerator(struct,
                                miller_index=index,
                                min_slab_size=self.min_slab_size,
                                min_vacuum_size=self.min_vacuum_size,
                                center_slab=False)
        slabs = {}
        for n, slab in enumerate(slabgen.get_slabs()):
            slab = AseAtomsAdaptor.get_atoms(slab)
            self.set_ase_cell(slab)
            a, b, c, alpha, beta, gamma = slab.get_cell_lengths_and_angles()
            supercell = [ceil(5/a), ceil(5/b), ceil(5/c)]
            # print(a, b, c, supercell)
            slab = slab*supercell
            slabs[n] = slab.copy()
        return slabs

    def set_ase_cell(self, surf):
        from numpy.linalg import norm
        a1, a2, a3 = surf.cell
        surf.set_cell([a1, a2,
                       np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /
                       norm(np.cross(a1, a2))**2])
        #
        a1, a2, a3 = surf.cell
        surf.set_cell([(norm(a1), 0, 0),
                       (np.dot(a1, a2) / norm(a1),
                        np.sqrt(norm(a2)**2 - (np.dot(a1, a2) / norm(a1))**2), 0),
                       (0, 0, norm(a3))],
                      scale_atoms=True)
        surf.wrap()

    def build_children(self, surfaces=None, nlayer=None, fix=None, tol=None):
        '''
        surfaces and terminations

        '''
        if not surfaces:
            surfaces = self.surfaces
        if not nlayer:
            nlayer = self.nlayer
        if not fix:
            fix = self.fix
        if not tol:
            tol = self.tol
        self.log('-'*60)
        self.log('Build surface:\n')
        count = 0
        for index, ters in self.indexs.items():
            count += 1
            self.log('{0:15s}: {1}'.format('%s. index' % count, index))
            # atoms = surface(self.atoms, index, nlayer + 2, tol = self.tol, vacuum=1)
            # atoms.pbc = [True, True, True]
            # print(self.tol)
            # layers = get_layers(atoms, (0, 0, 1), self.tol)[0]
            # print(layers)
            # terminations = self.get_terminations(atoms, layers)
            terminations = self.get_slabs(self.atoms, index)
            print(terminations)
            # if len(ters) == 0:
            # ters = terminations.keys()
            for ter, surf in terminations.items():
                # if ter not in ters: continue
                self.log('{0:15s}: {1}'.format('    termination', ter))
                # surf = atoms.copy()
                # index = [j for j in range(len(surf)) if layers[j] in range(ind, ind - nlayer, -1)]
                # surf = surf[index]
                # surf.positions[:, 2] -= min(surf.positions[:, 2]) - 0.01
                # surf.cell[2][2] = max(surf.positions[:, 2]) + 15
                surf = fix_layers(surf, (0, 0, 1), tol, fix)
                #
                prefix = '%s%s%s-%s' % (index[0], index[1], index[2], ter)
                self.images[prefix] = surf
                surf = OER_pourbaix(surf, label=os.path.join(self.label, prefix), prefix=prefix, calculator=self.calculator,
                                    view=self.view, molecule_energies=self.molecule_energies, **self.pourbaix_input, surface_input=self.surface_input)
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
        from ase.build import minimize_rotation_and_translation
        from ase.calculators.calculator import compare_atoms, equal
        import pprint
        nlayer = max(layers) + 1
        terminations = {}
        natoms = len(slab)
        for i in range(nlayer - 2, 0, -1):
            atoms = slab[layers == i]
            formula = atoms.get_chemical_formula(mode='hill')
            print(formula)
            inds = [j for j in range(natoms) if layers[j]
                    in range(i, i - self.nlayer, -1)]
            atoms = slab[inds]
            atoms.wrap()
            if not terminations:
                terminations['%s-%s' % (formula, i)] = atoms.copy()
            else:
                flag = True
                atoms1 = atoms.copy()
                for ter, atoms2 in terminations.items():
                    if len(atoms1) != len(atoms2):
                        flag = False
                    else:
                        minimize_rotation_and_translation(atoms2, atoms1)
                        if equal(atoms1.arrays, atoms2.arrays, tol=0.1):
                            flag = False
                if flag:
                    terminations['%s-%s' % (formula, i)] = atoms.copy()
                    print(formula, ter)
                    # view([atoms1, atoms2])
        pprint.pprint(terminations)
        # view(terminations.values())
        return terminations


class OER_pourbaix(OER_base):
    def __init__(self, atoms, label='', prefix=None, sites_dict=None, calculator=None, coverages=[0, 1], adsorbates=['O', 'OH'], surface_input={}, molecule_energies=None, view=False, debug=False):
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
        self.name = 'pourbaix'
        OER_base.__init__(self, atoms, label=label, prefix=prefix,
                          calculator=calculator, view=view, molecule_energies=molecule_energies)
        if not sites_dict:
            self.sites_dict = self.get_sites_dict(
                self.atoms, excludes=['O', 'N'])
        else:
            self.sites_dict = sites_dict
        self.coverages = coverages
        self.adsorbates = [''] + adsorbates
        self.surface_input = surface_input
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
        # fig, ax = plt.subplots(1, 1)
        ax, children = self.plot_pourbaix_diagram(self.results)  # , ax = ax)
        filename = os.path.join(
            self.label, '%s-pourbaix_diagram.png' % self.prefix)
        # fig.savefig(filename)
        self.log('-'*60)
        self.log('Save pourbaix diagram to: {0}'.format(filename))
        self.log('Stable surfaces from pourbaix diagram:')
        for prefix in children:
            surf = self.images[prefix]
            self.images[prefix] = surf
            surf = OER_surface(surf, label=os.path.join(self.label, prefix), prefix=prefix,
                               calculator=self.calculator, molecule_energies=self.molecule_energies, **self.surface_input)
            self.log('{0:15s}: {1}'.format('    Stable surface', prefix))
            self.children[prefix] = surf
        self.log()
        self.log('Total number of stable surfaces: {0}\n'.format(
            len(self.children)))
        # OER
        self.run_children()

    def run_coverage(self):
        self.pool_atoms(self.images)
        self.log('-'*60)
        self.log('Image energies (eV):')
        for job, result in self.results.items():
            self.log('{0:10s}  Energy: {1:1.3f}'.format(job, result['energy']))

    def build_coverage(self, coverages=None, sites_dict=None):
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
            self.log('    {0:3d}   {1:10s} {2}'.format(
                count, symbol, position))
            self.sites.append(position)
            site_mols.append(self.adsorbates)
        self.nsite = len(self.sites)
        self.log('Number of sites: {0}'.format(self.nsite))
        if len(sites_dict) == 0:
            return
        combinations = list(itertools.product(*site_mols))
        self.log('Coverages: {0}'.format(self.coverages))
        for combination in combinations:
            noh = combination.count('OH')
            no = combination.count('O')
            coverage = (no + noh)/self.nsite
            flag = [abs(coverage - x) for x in coverages]
            if min(flag) > 0.1:
                continue
            formula = '%s-O-%s-OH' % (no, noh)
            if formula in self.images:
                continue
            ads = Atoms('')
            for i in range(len(combination)):
                if not combination[i]:
                    continue
                pos = self.sites[i]
                mol = self.mols[combination[i]].copy()
                mol.translate(pos - mol[0].position + [0, 0, 1.5])
                ads = ads + mol
            natoms = self.atoms.copy()
            self.images[formula] = natoms + ads
            self.log('Structure: {0}'.format(formula))
        self.log('Total number of structures: {0}'.format(len(self.images)))

    def plot_pourbaix_diagram(self, results, ax=None, G_O=None, G_OH=None):
        '''
        mh = potential
        '''
        E_h2o = self.molecule_energies['H2O']
        E_h2 = self.molecule_energies['H2']
        if not G_O:
            G_O = E_h2o - E_h2
        if not G_OH:
            G_OH = E_h2o - E_h2/2.0
        potential = np.linspace(0, 2, 20)
        # if ax is None:
        # ax = plt.gca()
        E0 = results['0-O-0-OH']['energy']
        for job, result in results.items():
            E = result['energy']
            no = int(job.split('-')[0])
            noh = int(job.split('-')[2])
            dE = E - E0 - noh*(G_OH + potential) - no*(G_O + 2*potential)
            # ax.plot(potential, dE, label = job)
            self.results[job]['pourbaix_energy'] = dE
        # ax.legend()
        # ax.set_xlabel('Potential (V)')
        # ax.set_ylabel('$\Delta$ G (eV)')
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


class OER_surface(OER_base):
    def __init__(self, atoms, label='', prefix=None, activate_species=None, sites_dict=None, calculator=None, molecule_energies=None):
        '''
        sites_dict: dict
            the active site, should be the atom in the structure. e.g. {'Pt': [2]}. 
            If None, the ontop atom will be found and used.
            This is different from the site in OER_pourbaix.
        '''
        self.name = 'surface'
        OER_base.__init__(self, atoms, label=label, prefix=prefix,
                          calculator=calculator, molecule_energies=molecule_energies)
        self.atoms = atoms
        self.activate_species = activate_species
        if not sites_dict:
            self.sites_dict = self.get_sites_dict(self.atoms)
        else:
            self.sites_dict = sites_dict
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
        for prefix, pos in self.sites_dict.items():
            ele = prefix.split('-')[-1]
            if ele not in self.activate_species:
                continue
            self.images[prefix] = self.atoms
            child = OER_site(self.atoms, label=os.path.join(self.label, prefix), position=pos,
                             prefix=prefix, calculator=self.calculator, molecule_energies=self.molecule_energies)
            self.log('{0:15s}: {1}'.format('    sites', prefix))
            self.children[prefix] = child
            self.log()
        self.log('Total number of OER sites: {0}\n'.format(len(self.children)))


class OER_site(OER_base):
    def __init__(self, atoms, label='', prefix=None, site_type='ontop', site=None, height=2.0, calculator=None, molecule_energies=None, mols=None, view=False, debug=False):
        '''
        Workflow:
        0. Read the surface slab
        1. Build the adsorption site
        2. Add O*, OH* and OOH*, and relax the structure
        3. Calculate the dos and pdos
        4. Calculate the ZPE energy using Vibrations module
        5. Calculate the Gibbs free energy
        6. Generate the report for OER overpotential

        atoms:
            The surface slab, from previous optmization calculation.
        model_type: 
            (1) 'ontop', (2) 'bridge', (3) 'hollow'. For 'ontop' site, the 'O' and 'H' vacancy pathway are possible.
        site: 
            index for the adsoprtion. 
        height: 
            height of the O*, OH* and OOH* from the adsorption site.

        Examples:
        >>> oer = OER_site(atoms_opt,
                    label = 'oer/Pt-001-ontop',
                    site_type = 'ontop',
                    site = -1,
                    height=2.0,
                        calculator = parameters, 
                    molecule_energies = molecule_energies,
                    )
        >>> oer.run()

        '''
        self.name = 'site'
        OER_base.__init__(self, atoms, label=label, prefix=prefix, calculator=calculator,
                          molecule_energies=molecule_energies, view=view, debug=debug)
        self.images = {}
        self.site_type = site_type
        self.site = site
        self.height = height
        self.site_symbol = site_type
        if site_type.upper() == 'ONTOP':
            self.site_symbol = atoms[site].symbol
            self.site_position = atoms[site].position + \
                np.array([0, 0, height])
        elif site_type.upper() == 'BRIDGE':
            self.site_position = (
                atoms[site[0]].position + atoms[site[1]].position)/2.0 + np.array([0, 0, height])
        elif site_type.upper() == 'HOLLOW':
            self.site_position = (atoms[site[0]].position +
                                  atoms[site[1]].position +
                                  atoms[site[2]].position)/2.0 + np.array([0, 0, height])
        elif site_type.upper() == 'POSITION':
            self.site_position = site + np.array([0, 0, height])
        if mols:
            self.mols = mols
        self.children = {}
        self.calcs = {}
        self.zpes = {}
        self.free_energies = {}

    def run(self):
        """
        Function build all the main workflow
        """
        #
        self.build_oer_adsorbate()
        if self.view:
            self.view_images(self.children)
            return
        # geometry relaxiton
        self.pool_atoms(self.children, self.geo_relax_espresso)
        self.log('-'*60)
        self.log('Image energies (eV):')
        # self.log('{0: <20}  Energy: {1:1.3f}'.format('Clean Surface', self.atoms.calc.results['energy']))
        for job, result in self.results.items():
            self.log('    {0:<20}  {1:1.3f}'.format(job, result['energy']))
            self.children[job] = result['atoms']
        # pdos
        self.pool_atoms(self.children, self.pdos)
        # vibration calculation to get zero point energy
        self.pool_atoms(self.children, self.vib_zpe)
        self.log('-'*60)
        self.log('ZPE energies (eV):')
        for job, zpe in self.zpes.items():
            self.log('    {0:<20}  {1:1.3f}'.format(job, zpe))
        # calculate the oer over potential
        self.get_free_energy()
        self.write_image()
        self.log('-'*60)
        self.log('Gibbs free energies (eV): ')
        for step, e in zip(self.steps, self.free_energies):
            self.log('    {0:<10}    {1:1.2f}'.format(step, e))
        self.log('{0:10s}: {1:1.2f} eV'.format(
            'Over potential', self.over_potential))
        print('{0:10s}: {1:1.2f} eV'.format(
            'Over potential', self.over_potential))
        # view(self.children.values())
        return 'Over potential', self.over_potential

    def build_oer_adsorbate(self, ):
        """
        Add O*, OH* and OOH*, and relax the structure
        For 'ontop' site, the 'O' and 'H' vacancy pathway are possible.
        """
        self.log('-'*60)
        self.log('Adsoption site: index ({0})  symbol({1})'.format(
            self.site, self.site_symbol))
        self.log('OER adsorbate:')
        atoms = self.atoms
        label = self.label.replace('/', '_')
        # print(self.site_position)
        job = '%s_%s' % (label, 'Clean')
        self.children[job] = self.atoms.copy()
        if self.site_symbol not in ['H', 'O']:
            for prefix, mol in self.mols.items():
                # print(prefix, mol)
                job = '%s_%s' % (label, prefix)
                self.log('    %s' % job)
                ads = mol.copy()
                natoms = atoms.copy()
                ads.translate(self.site_position - ads[0].position)
                natoms = natoms + ads
                self.children[job] = natoms
        # O
        ind = self.site
        if self.site_symbol == 'O':
            job = '%s_%s' % (label, 'O')
            self.log('    %s' % job)
            self.children[job] = self.atoms.copy()
            atoms = self.atoms.copy()
            del atoms[[ind]]
            job = '%s_%s' % (label, 'Clean')
            self.log('    %s' % job)
            self.children[job] = atoms.copy()
            for prefix, mol in self.mols.items():
                # print(prefix, mol)
                if prefix == 'O':
                    continue
                job = '%s_%s' % (label, prefix)
                self.log('    %s' % job)
                ads = mol.copy()
                natoms = atoms.copy()
                ads.translate(self.atoms[ind].position - ads[0].position)
                natoms = natoms + ads
                self.children[job] = natoms
        elif self.site_symbol == 'H':
            atoms = self.atoms.copy()
            dis = atoms.get_distance(ind, range(len(atoms)))
            indo = dis.index(min(dis))
            # O
            job = '%s_%s' % (label, 'O')
            self.log('    %s' % job)
            natoms = atoms.copy()
            del natoms[[ind]]
            self.children[job] = natoms
            # O2
            job = '%s_%s' % (label, 'Clean')
            self.log('    %s' % job)
            natoms = atoms.copy()
            del natoms[[ind, indo]]
            self.children[job] = natoms
            # OOH
            job = '%s_%s' % (label, 'OOH')
            self.log('    %s' % job)
            ads = self.mols[job].copy()
            natoms = atoms.copy()
            ads.translate(self.atoms[indo].position - ads[0].position)
            natoms = natoms + ads
            self.children[job] = natoms
        self.log('Total number of OER adsorbate: {0}\n'.format(
            len(self.children)))

    def get_free_energy(self, results=None, G_O=None, G_OH=None, G_OOH=None, E_H=None):
        """
        """
        G_h2o = self.molecule_energies['H2O']
        G_h2 = self.molecule_energies['H2']
        if not G_O:
            G_O = G_h2o - G_h2
        if not G_OH:
            G_OH = G_h2o - G_h2/2.0
        if not G_OOH:
            G_OOH = 2*G_h2o - 3*G_h2/2.0
        if not E_H:
            E_H = G_h2/2.0
        # fig, ax = plt.subplots()
        label = self.label.replace('/', '_')
        E0 = self.results['%s_%s' % (label, 'Clean')]['energy']
        # if self.atoms[self.site].symbol not in []:
        dG_OH = self.results['%s_%s' % (
            label, 'OH')]['energy'] - E0 - G_OH + self.zpes['%s_%s' % (label, 'OH')]
        dG_O = self.results['%s_%s' % (
            label, 'O')]['energy'] - E0 - G_O + self.zpes['%s_%s' % (label, 'O')]
        dG_OOH = self.results['%s_%s' % (
            label, 'OOH')]['energy'] - E0 - G_OOH + self.zpes['%s_%s' % (label, 'OOH')]
        dG_O2 = 4.92
        self.steps = ['OH', 'O', 'OOH', 'O2']
        self.free_energies = [dG_OH, dG_O, dG_OOH, dG_O2]
        # self.plot_free_energy([0, G_OH, G_O, G_OOH, G_O2], ['*', 'OH', 'O', 'OOH', 'O2'], ax = ax)
        # if self.atoms[self.site].symbol == 'O':
        #     dG_OOH = self.results['%s_%s'%(label, 'OOH')]['energy'] - E0 - G_OH + self.zpes['%s_%s'%(label, 'OOH')]
        #     dG_O2 = self.results['%s_%s'%(label, 'O2')]['energy'] - E0 + G_O + self.zpes['%s_%s'%(label, 'O2')]
        #     dG_OH = self.results['%s_%s'%(label, 'OH')]['energy'] - E0 - E_H + self.zpes['%s_%s'%(label, 'OH')]
        #     dG_O = 4.92
        #     self.steps = ['OOH', 'O2', 'OH', 'O']
        #     self.free_energies = [dG_OOH, dG_O2, dG_OH, dG_O]
        #     # self.plot_free_energy([0, G_OOH, G_O2, G_OH, G_O], ['O', 'OOH', 'O2', 'OH', 'O'], ax = ax)
        # if self.atoms[self.site].symbol == 'H':
        #     dG_O = self.results['%s_%s'%(label, 'OH')]['energy'] - E0 + E_H + self.zpes['%s_%s'%(label, 'OH')]
        #     dG_OOH = self.results['%s_%s'%(label, 'O')]['energy'] - E0 - G_O + self.zpes['%s_%s'%(label, 'O')]
        #     dG_O2 = self.results['%s_%s'%(label, 'O2')]['energy'] - E0 + G_OH + self.zpes['%s_%s'%(label, 'O2')]
        #     dG_OH = 4.92
        #     self.steps = ['O', 'OOH', '2', 'OH']
        #     self.free_energies = [dG_O, dG_OOH, dG_O2, dG_OH]
        #     # self.plot_free_energy([0, G_O, G_OOH, G_O2, G_OH], ['OH', 'O', 'OOH', 'O2', 'OH'], ax = ax)
        # fig.savefig(os.path.join(self.label, '%s-free-energy.png'%self.prefix))
        self.over_potential = max([self.free_energies[1] - self.free_energies[0],
                                   self.free_energies[2] -
                                   self.free_energies[1],
                                   self.free_energies[3] - self.free_energies[2]]) - 1.23

    def write_image(self, ):
        for job in self.results:
            atoms = self.results[job]['atoms']
            write(os.path.join(self.label, '%s-top.png' %
                               job), atoms * (2, 2, 1))
            write(os.path.join(self.label, '%s-side.png' %
                               job), atoms * (2, 2, 1), rotation='-90x')
            write(os.path.join(self.label, '%s.png' % job),
                  atoms * (2, 2, 1), rotation='10z,-80x')

    def plot_free_energy(self, E, labels, ax=None):
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
