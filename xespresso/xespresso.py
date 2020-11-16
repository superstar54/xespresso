"""Quantum ESPRESSO Calculator

export PYTHONPATH="Your-Location":$PYTHONPATH
export ASE_ESPRESSO_COMMAND="/path/to/PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo"
export ESPRESSO_PSEUDO="/path/to/pseudo"

Run PACKAGE.x jobs.
"""


from ase import io
from ase.calculators.calculator import FileIOCalculator, CalculationFailed, equal, compare_atoms
import ase.calculators.espresso
from xespresso.xio import write_espresso_in, read_espresso_input, read_espresso_asei, write_espresso_asei, get_atomic_species
import os
import shutil
import numpy as np
from datetime import datetime
import copy
config_files = [os.path.join(os.environ['HOME'], '.xespressorc'),
            '.xespressorc']



class Espresso(ase.calculators.espresso.Espresso):
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms', 'time']
    command = 'PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='xespresso', prefix = None, atoms=None, package = 'pw', parallel = '',
                 queue = None, **kwargs):
        """
        Accepts all the options for pw.x as given in the QE docs, plus some
        additional options:

        input_data, pseudopotentials, kspacing, kpts, koffset
            Please have a look at Espresso module in ASE
        queue: dict
            A dictionary with parameters for job submission, e.g.
             queue = {'nodes': 4, 'ntasks-per-node': 20, 
                      'account': 'xxx', 'partition': 'normal', 
                      'time': '23:59:00'}
        package: str
            Choose the quantum espresso pacakge: pw, dos, projwfc, band, pp, ph, ..
            For NEB calculation, please use neb.NEBEspresso module.
            Calculaiton use phonon is not implemented yet.
        parallel: str
            A str which control the parallelization parameters: -nimage, -npools, 
            -nband, -ntg, -ndiag or -northo (shorthands, respectively: 
            -ni, -nk, -nb, -nt, -nd).

        General working procedure is as follows:
            1. Perform a regular self-consistent calculation:
              >>> calc = Espresso(input_data=input_data, ...)
              >>> atoms.set_calculator(calc)
              >>> atoms.get_potential_energy()
            2. post calculation can be made as follows:
              >>> calc.post('pacakge' = 'dos', **kwargs)
              >>> calc.post('pacakge' = 'pp', **kwargs)
            3. non-self-consistent calculation can be made as follows:
              >>> calc.nscf(calculation = 'nscf', kpts=(4, 4, 4))
              >>> calc.nscf_calculate()
        """
        self.set_label(restart, label, prefix)
        kwargs = self.set_kwargs(kwargs)
        FileIOCalculator.__init__(self, restart = self.directory, ignore_bad_restart_file = True,
                                  label = self.label, atoms = atoms, **kwargs)
        if atoms:
            self.atoms = atoms
        self.set_queue(package = package, parallel = parallel, queue = queue)
        self.parallel = parallel
        
    def set_label(self, restart, label, prefix):
        '''
        '''
        if restart:
            self.directory = os.path.split(label)[0]
            self.prefix = os.path.split(label)[1]
        else:
            self.directory = label
            if not prefix:
                self.prefix = os.path.split(label)[1]
            else:
                self.prefix = prefix
        self.label = os.path.join(self.directory, self.prefix)
        self.pwi = os.path.join(self.directory, '%s.pwi'%self.prefix)
        self.pwo = os.path.join(self.directory, '%s.pwo'%self.prefix)
        self.asei = os.path.join(self.directory, '%s.asei'%self.prefix)
        self.asei_temp = os.path.join(self.directory, '.%s.asei_temp'%self.prefix)
        self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
        self.scf_directory = None
        self.scf_parameters = None
        self.scf_results = None
    def set_kwargs(self, kwargs):
        
        if 'input_data' not in kwargs:
            kwargs['input_data'] = {}
        ase_parameters = copy.deepcopy(kwargs)
        for key, value in kwargs.items():
            if key not in ['pseudopotentials', 'kpts', 'kspacing', 'koffset', 'input_data', 'climbing_images', 'path_data']:
                ase_parameters['input_data'][key] = value
                del ase_parameters[key]
        ase_parameters['input_data']['prefix'] = self.prefix
        ase_parameters['input_data']['verbosity'] = 'high'
        self.ase_parameters = ase_parameters
        return ase_parameters
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
            script = ''
            for cf in config_files:
                if os.path.exists(cf):
                    file = open(cf, 'r')
                    script = file.read()
            jobname = self.prefix
            with open('%s/.job_file' % self.directory, 'w') as fh:
                fh.writelines("#!/bin/bash\n")
                fh.writelines("#SBATCH --job-name=%s \n" % jobname)
                fh.writelines("#SBATCH --output=%s.out\n" % self.prefix)
                fh.writelines("#SBATCH --error=%s.err\n" % self.prefix)
                fh.writelines("#SBATCH --wait\n")
                for key, value in queue.items():
                    if value:
                        fh.writelines("#SBATCH --%s=%s\n" %(key, value))
                fh.writelines("%s \n"%script)
                fh.writelines("%s \n" % command)
            self.command = "sbatch {0}".format('.job_file')
        else:
            self.command = command
    def read(self, label):
        """Read the files in a calculation if they exist.
        This function reads self.parameters and atoms.
        """
        # Else read a regular calculation. we start with reading stuff
        # that is independent of the calculation state.
        self.directory = label
        try:
            atoms, parameters = read_espresso_asei(self.asei)
            self.restart_atoms = atoms
            self.restart_parameters = parameters
            self.read_results()
        except:
            pass
    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        
    def check_state(self, atoms, tol=1e-4):
        """Check for any system changes since last calculation."""
        if not os.path.exists(self.asei):
            return True
        if isinstance(self.restart_atoms, list):
            if len(self.restart_atoms) != len(atoms):
                return True
            for i in range(len(atoms)):
                if compare_atoms(self.restart_atoms[i], atoms[i], tol=tol,
                             excluded_properties=set(self.ignored_changes)):
                    return True
        else:
            if compare_atoms(self.restart_atoms, atoms, tol=tol,
                             excluded_properties=set(self.ignored_changes)):
                return True
        if self.restart_parameters != self.ase_parameters:
            return True
        #print('Same geometry and parameters, try to check the results are available or not.')
        return False

    def read_results(self):
        '''
        get atomic species
        '''
        pwo = self.prefix + '.pwo'
        pwos = [file for file in os.listdir(self.directory) if pwo in file]
        output = None
        for pwo in pwos:
            atomic_species = None
            pwo = os.path.join(self.directory, pwo)
            atomic_species = get_atomic_species(pwo)
            try:
                output = io.read(pwo)
                atomic_species = get_atomic_species(pwo)
                if atomic_species: output.info['species'] = atomic_species
            except Exception as e:
                print(pwo, e)
        if output:
            self.calc = output.calc
            self.results = output.calc.results
            self.results['atoms'] = output
            self.efermi = self.get_fermi_level()
            self.nspins = self.get_number_of_spins()
        else:
            print('\nRead result failed\n')
    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write_espresso_in(self.label + '.pwi', atoms, **self.parameters)
        write_espresso_asei(self.label + '.asei', atoms, self.parameters)

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
            except Exception as e:
                print('Not converge: %s'%e)
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
            except Exception as e:
                print('Not converge: %s'%e)
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
    def nscf(self, calculation = 'nscf', package = 'pw', queue = False, parallel = '', kpts = (10, 10, 10), **kwargs):
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
    def post(self, queue = False, package = 'dos', parallel = '', **kwargs):
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
        self.set_queue(package = package, parallel = parallel, queue = queue)
        command = self.command
        print('Running %s'%package)
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
        print('Done: %s' % package)
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
        files = os.listdir(self.save_directory)
        keys = ['wfc']
        for file in files:
            for key in keys:
                if key in file:
                    os.remove(os.path.join(self.save_directory, file))
