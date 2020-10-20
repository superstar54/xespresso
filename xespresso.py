"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_COMMAND="/path/to/package.x -in PREFIX.pwi > PREFIX.pwo"

Run pw.x jobs.
"""


import warnings
from ase import io
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent, CalculationFailed
import ase.calculators.espresso
from xespresso.xio import write_espresso_in
from xespresso.xespressorc import xespressorc
import os
import shutil
import numpy as np
import pickle
from datetime import datetime


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'
# 
    
class Espresso(ase.calculators.espresso.Espresso):
    """
    """
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms', 'time']
    command = 'pw.x -in PREFIX.pwi > PREFIX.pwo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='xespresso', atoms=None, package = 'pw', parallel = '',
                 queue = None,
                 xc = 'pbe', pseudo = 'kjpaw',
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
        package: str
            Choose the quantum espresso pacakge: pw, dos, band, ph, ..
            neb is not implemented yet.
        parallel: str
            A str which control the parallelization parameters: -nimage, -npools, 
            -nband, -ntg, -ndiag or -northo (shorthands, respectively: 
            -ni, -nk, -nb, -nt, -nd).


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


        """
        
        
        
        self.directory = os.path.split(label)[0]
        self.prefix = os.path.split(label)[1]
        if 'input_data' not in kwargs:
            kwargs['input_data'] = {}
        kwargs['input_data']['prefix'] = self.prefix
        ase.calculators.espresso.Espresso.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.calc = None
        if atoms:
            self.atoms = atoms                                
        self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
        self.set_queue(package = package, parallel = '', queue = queue)
        self.package = package
        self.parallel = parallel
        self.scf_directory = None
        self.scf_parameters = None
        self.scf_results = None
    def set_queue(self, package = 'pw', parallel = '', queue = None):
        '''
        If queue, change command to "sbatch .job_file".
        The queue information are written into '.job_file'
        '''
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        command = os.environ.get('ASE_ESPRESSO_COMMAND')
        if 'PACKAGE' in command:
            if 'pw' in package:
                command = command.replace('PACKAGE', package, 1)
                command = command.replace('PACKAGE', 'pw', 2)
            else:
                command = command.replace('PACKAGE', package)
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)
        if 'PARALLEL' in command:
            command = command.replace('PARALLEL', parallel)
        if queue:
            # Write the job file
            xespressorc['queue'].update(queue)
            jobname = self.prefix
            # jobname = '%s-%s-%s' %(self.prefix, package, 
                # self.parameters['input_data']['calculation'])
            with open('%s/.job_file' % self.directory, 'w') as fh:
                fh.writelines("#!/bin/bash\n")
                fh.writelines("#SBATCH --job-name=%s \n" % jobname)
                fh.writelines("#SBATCH --output=%s.out\n" % self.prefix)
                fh.writelines("#SBATCH --error=%s.err\n" % self.prefix)
                fh.writelines("#SBATCH --wait\n")
                for key, value in xespressorc['queue'].items():
                    if value:
                        fh.writelines("#SBATCH --%s=%s\n" %(key, value))
                fh.writelines("%s \n"% xespressorc['script'])
                fh.writelines("%s \n" % command)
            self.command = "sbatch {0}".format('.job_file')
        else:
            self.command = command
        # print('Command: {0}'.format(self.command))
    def read_results(self):
        pwo = self.prefix + '.pwo'
        pwos = [file for file in os.listdir(self.directory) if pwo in file]
        output = None
        for pwo in pwos:
            pwo = os.path.join(self.directory, pwo)
            # print('read results from: %s'%pwo)
            try:
                output = io.read(pwo)
            except:
                pass
                # print('\nread %s failed\n'%pwo)
        if output:
            self.calc = output.calc
            self.results = output.calc.results
            self.results['atoms'] = output
            self.efermi = self.get_fermi_level()
            self.nspins = self.get_number_of_spins()
            self.atoms = output
        else:
            print('\nread result failed\n')
    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write_espresso_in(self.label + '.pwi', atoms, **self.parameters)

    def read_xml_file(self):
        '''
        If data-file-schema.xml is readable, then fromfile = True
        '''
        self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
        if not os.path.exists(self.save_directory): return False
        xml0 = os.path.join(self.save_directory, 'data-file-schema.xml')
        xmls = [file for file in os.listdir(self.save_directory) if 'data-file-schema.xml' in file]
        goodxml = []
        for xml in xmls:
            xml = os.path.join(self.save_directory, xml)
            # print('read results from: %s'%xml)
            with open(xml, 'r') as f:
                lines = f.readlines()
                if len(lines) < 1: continue
                if '</qes:espresso>' in lines[-1]:
                    goodxml.append(xml)
        fromfile = False
        nxml = len(goodxml)
        # print(nxml, goodxml, xml0)
        if nxml > 0:
            if xml0 not in goodxml:
                shutil.copy(goodxml[-1], xml0)
            fromfile = True
        return fromfile
    def read_convergence(self):
        '''
        Read the status of the calculation.
        {
        '0': 'Done',
        '-1': 'Cancelled',
        '1': 'Need restart'
        '2': 'Pending or not submit',
        '3': 'Not start',
        '4': 'Not finish',
        }
        '''
        fromfile = self.read_xml_file()
        errfile = self.label + '.err'
        if not os.path.exists(errfile):
            # print('%s not exists'%errfile)
            return 1, fromfile, 'Pending or not submit'
        errs = ['RESTART', 'pw.x', 'out-of-memory', 
                'NODE FAILURE', 'TIME LIMIT', 'COMMUNICATION', 
                'segmentation fault', 'PARSE_ERR', 'mpirun']
        cancelled = False
        with open(errfile, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'CANCELLED' in line:
                    cancelled = True
                for err in errs:
                    if err in line:
                        return 1, fromfile, line
            # cancelled by owner
            if len(lines) > 0:
                if cancelled:
                    return -1, fromfile, lines[0]
        # 
        output = self.label + '.pwo'
        if not os.path.exists(output):
            # print('%s not exists'%output)
            return 3, fromfile, 'No pwo output file'
        with open(output, 'r') as f:
            lines = f.readlines()
            if len(lines) == 0: return 1, fromfile, 'pwo file has nothing'
            stime = lines[1].split('starts on')[1]
            nlines = len(lines)
            n = min([100, nlines])
            for line in lines[-n:-1]:
                if line.rfind('too many bands are not converged') > -1:
                    return 1, fromfile, '%s, %s'%(stime, line)
                if line.rfind('convergence NOT achieved after') > -1:
                    fromfile = True
                    return 1, fromfile, '%s, %s'%(stime, line)
                if line.rfind('Maximum CPU time exceeded') > -1:
                    fromfile = True
                    return 1, fromfile, '%s, %s'%(stime, line)
                if line.rfind('JOB DONE.') > -1:
                    fromfile = True
                    return 0, fromfile, line
        return 4, fromfile, line
    def run(self, atoms = None, restart = False):
        '''
        run and restart
        '''
        if not restart:
            try:
                self.calculate()
            except:
                print('Not converge.')
        converged, fromfile, meg0 = self.read_convergence()
        print(converged, fromfile, meg0)
        i = 0
        # restart 0, 1, 2
        while converged > 0 or restart == 2:
            restart = 1
            # backup old files, *pwo, *xml
            self.backup_file('%s.pwo'%self.prefix, directory = self.directory)
            self.backup_file('data-file-schema.xml', directory = self.save_directory)
            print('RESTART %s'%i)
            if fromfile:
                self.parameters['input_data']['restart_mode'] = 'restart'
                # self.parameters['input_data']['startingwfc'] = 'file'
            else:                
                self.parameters['input_data']['restart_mode'] = 'from_scratch'
            try:
                atoms.get_potential_energy()
            except:
                print('Not converge.')
            converged, fromfile, meg = self.read_convergence()
            print(converged, fromfile, meg)
            if meg == meg0:
                print('Same message! \n %s'%meg)
                break
            meg0 = meg
            i += 1
    def backup_file(self, src, directory = '.'):
        '''
        compare files
        backup
        '''
        import filecmp
        tnow = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        new_src = os.path.join(directory, '%s-%s'%(tnow, src))
        src = os.path.join(directory, src)
        if not os.path.exists(src): return 0
        flist = []
        ftype = src.split('.')[-1]
        flag = True
        for dst in os.listdir(directory):
            dst = os.path.join(directory, dst)
            if dst == src: continue
            if dst.endswith(ftype):
                if filecmp.cmp(dst, src):
                    flag = False
        # print('backup files: ', flag, src, new_src)
        if flag:
            shutil.copy(src, new_src)
        return 0
    def nscf(self, calculation = 'nscf', package = None, queue = False, parallel = '', kpts = (10, 10, 10), **kwargs):
        import copy
        # save scf parameters
        if not self.scf_directory:
            print(self.directory)
            self.scf_directory = self.directory
        if not self.scf_parameters:
            self.scf_parameters = copy.deepcopy(self.parameters)
        if not self.scf_results:
            self.scf_results = copy.deepcopy(self.results)
            self.efermi = self.get_fermi_level()
        # read new atoms from results and set magmoms
        self.atoms = self.results['atoms']
        self.parameters = copy.deepcopy(self.scf_parameters)
        self.parameters['kpts'] = kpts
        self.parameters['calculation'] = calculation
        self.parameters['verbosity'] = 'high'
        self.parameters['restart_mode'] = 'restart'
        for key, value in kwargs.items():
            self.parameters[key] = value
        # create working directory
        self.save_directory_old = os.path.join(self.scf_directory, '%s.save'%self.prefix)
        self.directory = os.path.join(self.scf_directory, '%s/'%calculation)
        self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
        self.label = os.path.join(self.directory, self.prefix)
        if not package:
            package = self.package
        self.set_queue(package = package, parallel = parallel, queue = queue)
    def nscf_calculate(self, ):
        import shutil
        if not os.path.exists(self.save_directory):
            os.makedirs(self.save_directory)
        files = ['charge-density.hdf5', 'charge-density.dat', 'data-file-schema.xml', 'paw.txt', 'occup.txt']
        for species, pseudo in self.parameters['pseudopotentials'].items():
            files.append(pseudo)
        for file in files:
            file = os.path.join(self.save_directory_old, '%s'%file)
            if os.path.exists(file):
                # print(file)
                shutil.copy(file, self.save_directory)
        print('-'*30)
        print('\n {0} calculation in {1} ......'.format(
                    self.parameters['calculation'], 
                    self.directory))
        self.calculate()
        # self.read_results()

    def get_time(self, ):
        t = 0
        filename = os.path.join(self.directory, '%s.pwo'%self.prefix)
        with open(filename) as f:
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
        'hp': {'INPUTHP': ['prefix','outdir','iverbosity','max_seconds','nq1','nq2',
                           'nq3','skip_equivalence_q','determine_num_pert_only',
                           'find_atpert','docc_thr','skip_type','equiv_type',
                           'perturb_only_atom','start_q','last_q','sum_pertq',
                           'compute_hp','conv_thr_chi','thresh_init','ethr_nscf',
                           'niter_max','alpha_mix(i)','nmix','num_neigh','lmin','rmax',]},
        'bands': {'BANDS': ['prefix', 'outdir', 'filband', 'spin_component', 'lsigma', 
                            'lp', 'filp', 'lsym', 'no_overlap', 'plot_2d', 'firstk', 'lastk']},
        'ph': {'INPUTHP': ['amass', 'outdir', 'prefix', 'niter_ph', 'tr2_ph', 
                           'alpha_mix(niter)', 'nmix_ph', 'verbosity', 'reduce_io', 
                           'max_seconds', 'fildyn', 'fildrho', 'fildvscf', 'epsil', 
                           'lrpa', 'lnoloc', 'trans', 'lraman', 'eth_rps', 'eth_ns', 
                           'dek', 'recover', 'low_directory_check', 'only_init', 'qplot', 
                           'q2d', 'q_in_band_form', 'electron_phonon', 'el_ph_nsigma', 
                           'el_ph_sigma', 'ahc_dir', 'ahc_nbnd', 'ahc_nbndskip', 
                           'skip_upperfan', 'lshift_q', 'zeu', 'zue', 'elop', 'fpol', 
                           'ldisp', 'nogg', 'asr', 'ldiag', 'lqdir', 'search_sym', 
                           'nq1', 'nq2', 'nq3', 'nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3', 
                           'diagonalization', 'read_dns_bare', 'ldvscf_interpolate', 
                           'wpot_dir', 'do_long_range', 'do_charge_neutral', 'start_irr', 
                           'last_irr', 'nat_todo', 'modenum', 'start_q', 'last_q', 
                           'dvscf_star', 'drho_star', ]},
        
        }
        kwargs['prefix'] = self.prefix
        package_parameters = post_parameters[package]
        filename = os.path.join(self.directory, '%s.%si' %(self.prefix, package))
        with open(filename, 'w') as f:
            for section, parameters in package_parameters.items():
                f.write('&%s\n '%section)
                for key, value in kwargs.items():
                    if key in parameters:
                        f.write('  %s = %s, \n' %(key, value))
                f.write('/ \n')
        self.set_queue(package = package, parallel = '', queue = queue)
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
    def get_work_function(self, ax = None, inpfile = 'potential.cube', output = None, shift = False):
        import matplotlib.pyplot as plt
        from ase.io.cube import read_cube_data
        from ase.units import create_units
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        units = create_units('2006')
        #
        filename = os.path.join(self.directory, inpfile)
        data, atoms = read_cube_data(filename)
        data = data * units['Ry']
        ef = self.get_fermi_level()
        # x, y, z, lp = calc.get_local_potential()
        nx, ny, nz = data.shape
        axy = np.array([np.average(data[:, :, z]) for z in range(nz)])
        # setup the x-axis in realspace
        uc = atoms.get_cell()
        xaxis = np.linspace(0, uc[2][2], nz)
        if shift:
            ax.plot(xaxis, axy - ef)
            ef = 0
        else:
            ax.plot(xaxis, axy)
            ax.plot([min(xaxis), max(xaxis)], [ef, ef], 'k:')
        ax.set_xlabel('Position along z-axis ($\AA$)')
        ax.set_ylabel('Potential (eV)')
        if output:
            plt.savefig('%s'%output)
        # plt.show()
        atoms = self.results['atoms']
        pos = max(atoms.positions[:, 2] + 3)
        ind = (xaxis > pos) & (xaxis < pos + 3)
        wf = np.average(axy[ind]) - ef
        print('min: %s, max: %s'%(pos, pos + 3))
        print('The workfunction is {0:1.2f} eV'.format(wf))
    def get_bader_charge(self, inpfile = None):
        '''
        '''
        if not inpfile:
            inpfile = '%s.cube' %self.prefix
        command = 'bader %s'%inpfile
        print(command)
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err
        acf = os.path.join(self.directory, 'ACF.dat')


    def clean(self):
        '''
        remove wfc, hub files
        '''
        keys = ['.wfc', '.hub']
        files = os.listdir(self.directory)
        for file in files:
            for key in keys:
                if key in file:
                    os.remove(os.path.join(self.directory, file))
