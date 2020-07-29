"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"

Run pw.x jobs.
"""


import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent
from ase.calculators.espresso import Espresso
import os
import numpy as np
import pickle


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'
# 
    
class XEspresso(Espresso):
    """
    """
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms', 'time']
    command = 'pw.x -in PREFIX.pwi > PREFIX.pwo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='xespresso', atoms=None, 
                 queue = None,
                 **kwargs):
        """
        All options for pw.x are copied verbatim to the input file, and put
        into the correct section. Use ``input_data`` for parameters that are
        already in a dict, all other ``kwargs`` are passed as parameters.

        Accepts all the options for pw.x as given in the QE docs, plus some
        additional options:

        input_data: dict
            A flat or nested dictionary with input parameters for pw.x
        pseudopotentials: dict
            A filename for each atomic species, e.g.
            ``{'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}``.
            A dummy name will be used if none are given.
        kspacing: float
            Generate a grid of k-points with this as the minimum distance,
            in A^-1 between them in reciprocal space. If set to None, kpts
            will be used instead.
        kpts: (int, int, int), dict, or BandPath
            If kpts is a tuple (or list) of 3 integers, it is interpreted
            as the dimensions of a Monkhorst-Pack grid.
            If kpts is a dict, it will either be interpreted as a path
            in the Brillouin zone (*) if it contains the 'path' keyword,
            otherwise it is converted to a Monkhorst-Pack grid (**).
            (*) see ase.dft.kpoints.bandpath
            (**) see ase.calculators.calculator.kpts2sizeandoffsets
        koffset: (int, int, int)
            Offset of kpoints in each direction. Must be 0 (no offset) or
            1 (half grid offset). Setting to True is equivalent to (1, 1, 1).
        queue: dict
            A dictionary with parameters for job submission, e.g.
             queue = {'nodes': 4, 'nta]sks-per-node': 20, 
                      'account': 'xxx', 'partition': 'normal', 
                      'time': '23:59:00'}


        .. note::
           Set ``tprnfor=True`` and ``tstress=True`` to calculate forces and
           stresses.

        .. note::
           Band structure plots can be made as follows:


           1. Perform a regular self-consistent calculation,
              saving the wave functions at the end, as well as
              getting the Fermi energy:

              >>> input_data = {<your input data>}
              >>> calc = Espresso(input_data=input_data, ...)
              >>> atoms.set_calculator(calc)
              >>> atoms.get_potential_energy()
              >>> fermi_level = calc.get_fermi_level()

           2. Perform a non-self-consistent 'band structure' run
              after updating your input_data and kpts keywords:

              >>> calc.nscf(run_type = 'bands', kpts={<your Brillouin zone path>})
              >>> calc.nscf_calculate(atoms)

           3. Make the plot using the BandStructure functionality,
              after setting the Fermi level to that of the prior
              self-consistent calculation:

              >>> bs = calc.band_structure()
              >>> bs.reference = fermi_energy
              >>> bs.plot()

        """
        Espresso.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.calc = None
        if atoms:
            self.atoms = atoms
        self.prefix = label.split('/')[-1]
        if 'input_data' in kwargs:
            # if 'prefix' not in kwargs['input_data']:
            kwargs['input_data']['prefix'] = self.prefix
        self.set_queue(queue)
        self.scf_directory = None
        self.scf_parameters = None
        self.scf_results = None
    def set_queue(self, queue):
        ''' '''
        if queue:
            directory = './' + self.label[0:-len(self.prefix)]
            if not os.path.exists(directory):
                os.makedirs(directory)
            # Write the file
            command = os.environ.get('ASE_ESPRESSO_COMMAND')
            if 'PREFIX' in command:
                command = command.replace('PREFIX', self.prefix)
            with open('%s.job_file' % directory, 'w') as fh:
                fh.writelines("#!/bin/bash\n")
                fh.writelines("#SBATCH --job-name=%s \n" % self.prefix)
                fh.writelines("#SBATCH --output=%s.out\n" % self.prefix)
                fh.writelines("#SBATCH --error=%s.err\n" % self.prefix)
                fh.writelines("#SBATCH --wait\n")
                for key, value in queue.items():
                    if value:
                        fh.writelines("#SBATCH --%s=%s\n" %(key, value))
                fh.writelines("module load QuantumESPRESSO/6.5-iomkl-2018a \n")
                fh.writelines("%s \n" % command)
            self.command = "sbatch {0}".format('.job_file')
        else:
            self.command = os.environ.get('ASE_ESPRESSO_COMMAND')

    def read_results(self):
        output = io.read(self.label + '.pwo')
        self.calc = output.calc
        self.results = output.calc.results
        self.results['atoms'] = output
    def nscf(self, run_type = 'dos', queue = False, kpts = (10, 10, 10), **kwargs):
        import copy
        # save scf parameters
        if not self.scf_directory:
            self.scf_directory = self.directory
        if not self.scf_parameters:
            self.scf_parameters = copy.deepcopy(self.parameters)
        if not self.scf_results:
            self.scf_results = copy.deepcopy(self.results)
            self.efermi = self.get_fermi_level()
        # read new atoms from results and set magmoms
        self.atoms = self.results['atoms']
        if 'magmoms' in self.results:
            self.atoms.set_initial_magnetic_moments(self.results['magmoms'])
        self.parameters = copy.deepcopy(self.scf_parameters)
        self.parameters['kpts'] = kpts
        kwargs['verbosity'] = 'high'
        # nscf or bands
        if run_type.upper() == 'BANDS':
            kwargs['calculation'] = 'bands'
        elif run_type.upper() == 'DOS':
            kwargs['calculation'] = 'nscf'
            self.parameters['input_data']['occupations'] = 'tetrahedra'
            for key in ['smearing', 'degauss']:
                if key in self.parameters['input_data']:
                    del self.parameters['input_data'][key]
        for key, value in kwargs.items():
            self.parameters['input_data'][key] = value
        # create working directory
        self.directory = self.scf_directory + '/' + '%s/'%(kwargs['calculation'])
        self.label = self.directory + self.prefix
        self.set_queue(queue)
    def nscf_calculate(self, ):
        import shutil
        directory = self.directory + '/{0}.save'.format(self.prefix)
        if not os.path.exists(directory):
            os.makedirs(directory)
        files = ['charge-density.hdf5', 'charge-density.dat', 'data-file-schema.xml', 'paw.txt', 'occup.txt']
        for species, pseudo in self.parameters['pseudopotentials'].items():
            files.append(pseudo)
        for file in files:
            file = self.scf_directory + '/%s.save/%s' % (self.prefix, file)
            if os.path.exists(file):
                shutil.copy(file, directory)
        print('-'*30)
        print('\n {0} calculation in {1} ......'.format(
                    self.parameters['input_data']['calculation'], 
                    self.directory))
        self.calculate()
        # self.read_results()
    def calc_dos(self, Emin = -10, Emax = 10, DeltaE = 0.1, 
                 bz_sum = None, ngauss = None, sigma = None):
        # create input for dos.x
        with open(self.directory+'/%s-dos.inp' % self.prefix , 'w') as fpdos:
            fpdos.write('&DOS\n  prefix={0},\n  outdir=\'.\',\n'.format(self.prefix))
            fpdos.write('  Emin = %s,\n' %(Emin + self.efermi))
            fpdos.write('  Emax = %s,\n' %(Emax + self.efermi))
            fpdos.write('  DeltaE = %s, \n' % DeltaE)
            if bz_sum:
                fpdos.write('  lbz_sum =%s,\n' % bz_sum)
            if ngauss is not None:
                fpdos.write('  ngauss = %s,\n' %ngauss)
            if sigma is not None:
                fpdos.write('  degauss = %s,\n' %(sigma/Rydberg))
            fpdos.write('/')
        print('Running dos.x ......')
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system('dos.x -in {0}-dos.inp > {0}-dos.out'.format(self.prefix))
        os.chdir(cwd)
        self.read_dos()
    def read_dos(self):
        dos = np.loadtxt(self.directory+'/%s.dos' % self.prefix)
        self.nspins = self.get_number_of_spins()
        if self.nspins == 2:
            self.dos = [dos[:,1], dos[:, 2]]
        else:
            self.dos = [dos[:,1]]
        self.dos_energies = dos[:,0] - self.efermi
    def plot_dos(self, output = None, fill = True):
        import matplotlib.pyplot as plt
        self.nspins = self.get_number_of_spins()
        fig, ax = plt.subplots(figsize = (6, 3))
        label = 'total'
        lspins = ['up', 'dw']
        for i in range(self.nspins):
            if self.nspins == 2:    label = '%s' % lspins[i]
            plt.plot(self.dos_energies, (-1)**i*self.dos[i], linewidth=0.7, label = label)
            if fill:
                plt.fill_between(self.dos_energies, (-1)**i*self.dos[i], 0, alpha = 0.2)
        plt.legend()
        plt.grid(linestyle = '--')
        plt.xlabel('Energy (eV)')
        plt.ylabel('DOS (a.u.)')
        plt.title('%s' % self.prefix)
        if not output:
            output = '{0}-dos.png'.format(self.prefix) 
        plt.savefig('%s' %output)

    def calc_pdos(self, Emin = -10, Emax = 10, DeltaE = 0.1,
        ngauss = None, sigma = None, lwrite_overlaps = False,):
        # print('Fermi energy: ', self.efermi)
        # create input for projwfc.x
        with open(self.directory+'/%s-projwfc.inp' % self.prefix , 'w') as fpdos:
            fpdos.write('&PROJWFC\n  prefix={0},\n  outdir=\'.\',\n'.format(self.prefix))
            fpdos.write('  Emin = %s,\n' %(Emin + self.efermi))
            fpdos.write('  Emax = %s,\n' %(Emax + self.efermi))
            fpdos.write('  DeltaE = %s, \n' % DeltaE)
            fpdos.write('  filproj = %s, \n' % self.prefix)
            if ngauss is not None:
                fpdos.write('  ngauss = %s,\n' %ngauss)
            if sigma is not None:
                fpdos.write('  degauss = %s,\n' %(sigma/Rydberg))
            fpdos.write('/')
        # run projwfc.x
        print('Running projwfc.x ......')
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system('projwfc.x -in {0}-projwfc.inp > {0}-projwfc.out'.format(self.prefix))
        os.chdir(cwd)
        self.read_pdos()

    #
    def read_pdos_info(self, ):
        import re
        reg = re.compile('\w+')
        species = []
        pdos_info = {}
        orbitals = ['s', 'p', 'd', 'f']
        with open(self.directory+'/%s-projwfc.out' % self.prefix, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'state #' in line:
                    data = reg.findall(line)
                    istate = int(data[3])
                    iatom = int(data[3])
                    label = data[4]
                    ichannel = int(data[6])
                    l = int(data[8])
                    m = int(data[10])
                    # print(iatom, label, ichannel, l, m)
                    if label not in species:
                        species.append(label)
                        pdos_info[label] = {}
                        pdos_info[label]['iatom'] = []
                        pdos_info[label]['ichannel'] = {}
                        # pdos_info[label]['orbital'] = []
                    if iatom not in pdos_info[label]['iatom']:
                        pdos_info[label]['iatom'].append(iatom)
                    if ichannel not in pdos_info[label]['ichannel']:
                        pdos_info[label]['ichannel'][ichannel] = l
                        # pdos_info[label]['orbital'].append(orbitals[l])
        self.species = species
        self.pdos_info = pdos_info

    def read_pdos(self):
        self.read_pdos_info()
        # read in total density of states
        dos = np.loadtxt(self.directory+'/%s.pdos_tot' % self.prefix)
        self.nspins = self.get_number_of_spins()
        if self.nspins == 2:
            self.pdos_tot = [dos[:,1], dos[:, 2]]
        else:
            self.pdos_tot = [dos[:,1]]
        self.dos_energies = dos[:,0] - self.efermi
        npoints = len(self.dos_energies)
        orbitals = ['s', 'p', 'd', 'f']
        # read in projections onto atomic orbitals
        self.natoms = len(self.atoms)
        self.nspecies = len(self.species)
        self.pdos = []
        self.pdos_atom = []
        self.pdos_species = {}
        for species, info in self.pdos_info.items():
            pdos_species = {}
            for ichannel, l in info['ichannel'].items():
                ncomponents = (2*l+1) * self.nspins + 1
                # print(ncomponents)
                channel = '{0}{1}'.format(ichannel, orbitals[l])
                pdos_species[channel] = np.zeros((ncomponents, npoints), np.float)
            for iatom in info['iatom']:
                pdos_atom = {}
                for ichannel, l in info['ichannel'].items():
                    # print(species, ichannel, l)
                    filename = self.directory + '/{0}.pdos_atm#{1}({2})_wfc#{3}({4})'.format(self.prefix, iatom, species, ichannel, orbitals[l])
                    channel = '{0}{1}'.format(ichannel, orbitals[l])
                    pdosinp = np.genfromtxt(filename)
                    # print(inpfile, channel, jpos)
                    ncomponents = (2*l+1) * self.nspins + 1
                    pdos_atom[channel] = np.zeros((ncomponents, npoints), np.float)
                    # pdos_species[channel][0] += pdosinp[:, 0]
                    for j in range(ncomponents):
                        # print(j)
                        pdos_atom[channel][j] += pdosinp[:, j + 1]
                         # sum over species
                        pdos_species[channel][j] += pdosinp[:,j + 1]
                self.pdos_atom.append(pdos_atom)
            self.pdos_species[species] = pdos_species
        # sum over orbital
        # print(self.pdos_species)
        return self.dos_energies, self.pdos_tot, self.pdos_atom, self.pdos_species
    def read_proj(self, ):
        '''

        read projwfc output file *.up and *.dw
        
        '''
        #
        self.read_pdos_info()
        self.nspins = self.get_number_of_spins()
        proj_files = [self.directory + '/%s.projwfc_up'%self.prefix, self.directory + '/%s.projwfc_down'%self.prefix]
        projs = []
        projs_species = []
        channels = []
        for i in range(self.nspins):
            with open(proj_files[i]) as f:
                lines = f.readlines()
                for j in range(len(lines)):
                    if '       1       1' in lines[j]: break
                nstates = int(lines[j - 3].split()[0])
                nkpts = int(lines[j - 3].split()[1])
                nbands = int(lines[j - 3].split()[2])
                print(nstates, nkpts, nbands)
                proj = np.zeros((nstates, nkpts, nbands), np.float)
                proj_species = {}
                for istate in range(nstates):
                    iline = j + istate*(nkpts*nbands + 1) - 1
                    species = lines[iline].split()[2]
                    orbital = lines[iline].split()[3]
                    channel = species + orbital
                    if channel not in channels: 
                        channels.append(channel)
                        proj_species[channel] = np.zeros((nkpts, nbands), np.float)
                    # l = lines[iline].split()[2]
                    # print(species, ichannel)
                    for ikpt in range(nkpts):
                        for iband in range(nbands):
                            iline = j + istate*(nkpts*nbands + 1) + ikpt*nbands + iband
                            # print('line: ', iline)
                            proj[istate, ikpt, iband] = float(lines[iline].split()[2])
                            proj_species[channel][ikpt, iband] += float(lines[iline].split()[2])
                            # proj_atoms[channel, ikpt, iband] += 
            projs.append(proj)
            projs_species.append(proj_species)
        self.projs = projs
        self.projs_species = projs_species
    def plot_pdos_tot(self, fill = True, output = None):
        import matplotlib.pyplot as plt
        # print(self.pdos_tot)
        fig, ax = plt.subplots(figsize = (6, 3))
        self.nspins = self.get_number_of_spins()
        lspins = ['up', 'dw']
        label = 'total'
        for i in range(self.nspins):
            if self.nspins == 2:    label = '%s' % lspins[i]
            plt.plot(self.dos_energies, (-1)**i*self.pdos_tot[i], linewidth=0.7, label = label)
            if fill:
                plt.fill_between(self.dos_energies, (-1)**i*self.pdos_tot[i], 0, alpha = 0.2)
        plt.legend()
        plt.xlabel('Energy (eV)')
        plt.ylabel('PDOS (a.u.)')
        plt.title('%s' % self.prefix)
        if not output:
            output = '{0}-pdos-tot.png'.format(self.prefix) 
        plt.savefig('%s' %output)
    def plot_pdos(self, total = False, select = None, fill = True, output = None):
        import matplotlib.pyplot as plt
        self.nspins = self.get_number_of_spins()
        lspins = ['up', 'dw']
        fig, ax = plt.subplots(figsize = (6, 3))
        if total:
            label = 'total'
            for i in range(self.nspins):
                if self.nspins == 2:    label = 'total %s' % lspins[i]
                plt.plot(self.dos_energies, (-1)**i*self.pdos_tot[i], linewidth=0.7, label = label, color = 'grey')
                # if fill:
                    # plt.fill_between(self.dos_energies, (-1)**i*self.pdos_tot[i], 0, alpha = 0.2)
        for species, channels in self.pdos_species.items():
            if select and species not in select: continue
            for channel, pdos in channels.items():
                if select and channel[-1] not in select[species]: continue
                for i in range(self.nspins):
                    label = '{0}-{1}'.format(species, channel)
                    if self.nspins == 2:
                        label = '{0}-{1}-{2}'.format(species, channel[-1], lspins[i])
                    plt.plot(self.dos_energies, (-1)**i*pdos[i], linewidth=0.7, label = label)
                    if fill:
                        plt.fill_between(self.dos_energies, (-1)**i*pdos[i], 0, alpha = 0.2)
        plt.legend()
        plt.grid(linestyle = '--')
        plt.xlabel('Energy (eV)')
        plt.ylabel('PDOS (a.u.)')
        plt.title('%s' % self.prefix)
        if not output:
            output = '{0}-pdos.png'.format(self.prefix) 
        plt.savefig('%s' %output)
    def get_time(self, ):
        t = 0
        with open(self.directory + '/%s.pwo'%self.prefix) as f:
            lines = f.readlines()
            for line in lines[-200:]:
                if 'init_run     :' in line or 'electrons    :' in line:
                    t += float(line.split('CPU')[1].split('s WALL')[0])
        self.results['time'] = t
        return t
    def ppx(self, **kwargs):
        # create input for pp.x
        pp_parameters = {
            'INPUTPP': ['prefix', 'outdir', 'filplot', 'plot_num', 
                       'spin_component', 'spin_component', 'emin', 'emax', 
                       'delta_e', 'degauss_ldos', 'sample_bias', 'kpoint', 
                       'kband', 'lsign', 'spin_component', 'emin', 'emax', 
                       'spin_component', 'spin_component', 'spin_component', 
                       'spin_component'], 
            'PLOT': ['nfile', 'filepp', 'weight', 'iflag', 'output_format', 
                    'fileout', 'interpolation', 'e1', 'x0', 'nx', 'e1', 'e2', 
                    'x0', 'nx', 'ny', 'e1', 'e2', 'e3', 'x0', 'nx', 'ny', 'nz', 
                    'radius', 'nx', 'ny']
        }
        kwargs['prefix'] = self.prefix
        with open(self.directory+'/pp-%s-%s.inp' %(self.prefix, kwargs['plot_num']), 'w') as fpdos:
            for section, parameters in pp_parameters.items():
                fpdos.write('&%s\n '%section)
                for key, value in kwargs.items():
                    if key in parameters:
                        fpdos.write('  %s = %s, \n' %(key, value))
                fpdos.write('/ \n')
        # run pp.x
        print('Running pp.x ......')
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system('pp.x -in pp-{0}-{1}.inp > pp-{0}-{1}.out'.format(self.prefix, kwargs['plot_num']))
        os.chdir(cwd)
    def get_work_function(self, input = 'potential.cube', output = 'workfunction.png'):
        import matplotlib.pyplot as plt
        from ase.io.cube import read_cube_data
        from ase.units import create_units
        units = create_units('2006')
        #
        filename = self.directory + '/' + input
        data, atoms = read_cube_data(filename)
        data = data * units['Ry']
        ef = self.get_fermi_level()
        # x, y, z, lp = calc.get_local_potential()
        nx, ny, nz = data.shape
        axy = np.array([np.average(data[:, :, z]) for z in range(nz)])
        # setup the x-axis in realspace
        uc = atoms.get_cell()
        xaxis = np.linspace(0, uc[2][2], nz)
        plt.plot(xaxis, axy)
        plt.plot([min(xaxis), max(xaxis)], [ef, ef], 'k:')
        plt.xlabel('Position along z-axis')
        plt.ylabel('x-y averaged electrostatic potential')
        plt.savefig('%s'%output)
        # plt.show()
        atoms = self.results['atoms']
        pos = max(atoms.positions[:, 2] + 3)
        ind = (xaxis > pos) & (xaxis < pos + 3)
        wf = np.average(axy[ind]) - ef
        print('min: %s, max: %s'%(pos, pos + 3))
        print('The workfunction is {0:1.2f} eV'.format(wf))
    def read_convergence(self):
        converged = 'FAILED'
        output = self.label + '.pwo'
        if not os.path.exists(output):
            return 'NONE'
        with open(output, 'r') as f:
            lines = f.readlines()
        nlines = len(lines)
        n = min([100, nlines])
        for line in lines[-n:-1]:
            if line.rfind('Maximum CPU time exceeded') > -1:
                return 'RESTART'
            if line.rfind('JOB DONE.') > -1:
                converged = 'DONE'
        return converged


#====================================================
# tools
def summary(updates = [], prefix = 'datas'):
    datas = {}
    calc = Espresso()
    print('Reading.....')
    for update in updates:
        cwd = os.getcwd()
        for i,j,y in os.walk(update):
            output = is_espresso(i)
            if output:
                os.chdir(i)
                print('dire:', i)
                calc.directory = cwd + '/' + i
                calc.prefix = output[0:-4]
                try:
                    calc.results = {}
                    calc.read_results()
                    t = calc.get_time()
                    calc.results['time'] = t
                    datas[i] = calc.results
                except Exception as e:
                    print('='*30, '\n', i, e)
            os.chdir(cwd)
    with open('%s.pickle' % prefix, 'wb') as f:
        pickle.dump(datas, f)
    print('Finished')

def is_espresso(path):
    '''
    check espresso 
    '''
    dirs = os.listdir(path)
    # print(dirs)
    # flag = True
    for qefile in ['.pwo']:
        flag = False
        for file in dirs:
            if qefile in file:
                return file
        if not flag:
            return False
    # return flag