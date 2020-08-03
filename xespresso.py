"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_COMMAND="/path/to/pw.x -in PREFIX.pwi > PREFIX.pwo"

Run pw.x jobs.
"""


import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent, CalculationFailed
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
        self.directory = './' + label[0:-len(self.prefix)]
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.set_queue('pw', queue)
        self.scf_directory = None
        self.scf_parameters = None
        self.scf_results = None
    def set_queue(self, package = 'pw', queue = None):
        ''' '''
        command = os.environ.get('ASE_ESPRESSO_COMMAND')
        if 'PACKAGE' in command:
            command = command.replace('PACKAGE', package)
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)
        print('Command: {0}'.format(command))
        if queue:
            # Write the job file
            with open('%s/.job_file' % self.directory, 'w') as fh:
                fh.writelines("#!/bin/bash\n")
                fh.writelines("#SBATCH --job-name=%s \n" % self.prefix)
                fh.writelines("#SBATCH --output=%s.out\n" % self.prefix)
                fh.writelines("#SBATCH --error=%s.err\n" % self.prefix)
                fh.writelines("#SBATCH --wait\n")
                for key, value in queue.items():
                    if value:
                        fh.writelines("#SBATCH --%s=%s\n" %(key, value))
                fh.writelines("ulimit -s unlimited \n")
                fh.writelines("module load QuantumESPRESSO/6.5-iomkl-2018a \n")
                fh.writelines("%s \n" % command)
            self.command = "sbatch {0}".format('.job_file')
        else:
            self.command = command
        print('Command: {0}'.format(self.command))
    def read_results(self):
        output = io.read(self.label + '.pwo')
        self.calc = output.calc
        self.results = output.calc.results
        self.results['atoms'] = output
        self.efermi = self.get_fermi_level()
        self.nspins = self.get_number_of_spins()
        self.atoms = output
    def read_convergence(self):
        '''
        slurmstepd: error: *** JOB 31580279 ON anode150 CANCELLED AT 
        2020-07-29T18:45:16 DUE TO NODE FAILURE, 
        SEE SLURMCTLD LOG FOR DETAILS 
        '''
        '''
        ORTE has lost communication with a remote daemon.       

          HNP daemon   : [[16919,0],0] on node anode031
          Remote daemon: [[16919,0],3] on node anode127     

        This is usually due to either a failure of the TCP network
        connection to the node, or possibly an internal failure of
        the daemon itself. We cannot recover from this failure, and
        therefore will terminate the job.
        '''
        savedir = self.label + '.save'
        # print(savedir)
        fromfile = False
        if os.path.exists(savedir):
            fromfile = True
        converged = 'DONE'
        errfile = self.label + '.err'
        if not os.path.exists(errfile):
            return converged, fromfile, 'NO OUTPUT'
        errs = ['RESTART', 'pw.x', 'out-of-memory', 'NODE FAILURE', 'TIME LIMIT', 'COMMUNICATION']
        with open(errfile, 'r') as f:
            lines = f.readlines()
            # print(line)
            for line in lines:
                for err in errs:
                    if err in line:
                        converged = 'RESTART'
                        return converged, fromfile, line
            if 'CANCELLED' in lines[0]:
                    converged = 'CANCELLED'
                    return converged, fromfile, lines[0]
        # 
        output = self.label + '.pwo'
        if not os.path.exists(output):
            return converged, fromfile, 'NO OUTPUT'
        with open(output, 'r') as f:
            lines = f.readlines()
            stime = lines[1].split('starts on')[1]
            nlines = len(lines)
            n = min([100, nlines])
            for line in lines[-n:-1]:
                if line.rfind('Maximum CPU time exceeded') > -1:
                    fromfile = True
                    return 'RESTART', fromfile, '%sMaximum CPU time exceeded'%stime
        return converged, fromfile, 'JOB DONE'
    def run(self, atoms = None, restart = False):
        import datetime
        if not restart:
            try:
                atoms.get_potential_energy()
            except:
                print('Not converge.')
        converged, fromfile, meg0 = self.read_convergence()
        print(converged, fromfile, meg0)
        i = 0
        while converged == 'RESTART':
            print(datetime.datetime.now())
            print('RESTART %s'%i)
            if fromfile:
                self.parameters['input_data']['restart_mode'] = 'restart'
                self.parameters['input_data']['startingwfc'] = 'file'
            else:                
                self.parameters['input_data']['restart_mode'] = 'from_scratch'
            try:
                atoms.get_potential_energy()
            except:
                print('Not converge.')
            converged, fromfile, meg = self.read_convergence()
            print(converged, fromfile, meg)
            if meg == meg0: break
            meg0 = meg
            i += 1
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
        self.set_queue('pw', queue)
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

    def get_time(self, ):
        t = 0
        with open(self.directory + '/%s.pwo'%self.prefix) as f:
            lines = f.readlines()
            for line in lines[-200:]:
                if 'init_run     :' in line or 'electrons    :' in line:
                    t += float(line.split('CPU')[1].split('s WALL')[0])
        self.results['time'] = t
        return t
    def post(self, queue = False, package = 'dos', **kwargs):
        '''
        '''
        import subprocess
        post_parameters = {
        'dos':  {'DOS': ['prefix', 'outdir', 'bz_sum', 'ngauss', 'degauss', 
                        'Emin', 'Emax', 'DeltaE', 'fildo', ]},
        'pp':   {'INPUTPP': ['prefix', 'outdir', 'filplot', 'plot_num', 
                       'spin_component', 'spin_component', 'emin', 'emax', 
                       'delta_e', 'degauss_ldos', 'sample_bias', 'kpoint', 
                       'kband', 'lsign', 'spin_component', 'emin', 'emax', 
                       'spin_component', 'spin_component', 'spin_component', 
                       'spin_component'], 
                'PLOT': ['nfile', 'filepp', 'weight', 'iflag', 'output_format', 
                    'fileout', 'interpolation', 'e1', 'x0', 'nx', 'e1', 'e2', 
                    'x0', 'nx', 'ny', 'e1', 'e2', 'e3', 'x0', 'nx', 'ny', 'nz', 
                    'radius', 'nx', 'ny']},
        'projwfc': {'PROJWFC': ['prefix', 'outdir', 'ngauss', 'degauss', 
                    'Emin', 'Emax', 'DeltaE', 'lsym', 'pawproj', 'filpdos', 'filproj', 
                    'lwrite_overlaps', 'lbinary_data', 'kresolveddos', 'tdosinboxes', 
                    'n_proj_boxes', 'irmin(3,n_proj_boxes)', 'irmax(3,n_proj_boxes)', 
                    'plotboxes', ]},
        }
        kwargs['prefix'] = self.prefix
        package_parameters = post_parameters[package]
        with open(self.directory+'/%s.%si' %(self.prefix, package), 'w') as f:
            for section, parameters in package_parameters.items():
                f.write('&%s\n '%section)
                for key, value in kwargs.items():
                    if key in parameters:
                        f.write('  %s = %s, \n' %(key, value))
                f.write('/ \n')
        self.set_queue(package, queue)
        # self.post_calculate()
        # self.read_dos()
        # if not command:
        command = self.command
        print('running %s'%package)
        print(command)
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format(self.name, command,
                                                  path, errorcode))
            raise CalculationFailed(msg)
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
