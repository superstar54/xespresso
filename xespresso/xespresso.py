"""Quantum ESPRESSO Calculator

export PYTHONPATH="Your-Location"/xespresso:$PYTHONPATH
export ASE_ESPRESSO_COMMAND="/path/to/PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo"
export ESPRESSO_PSEUDO="/path/to/pseudo"

Run PACKAGE.x jobs.
"""


from ase import io
from ase.calculators.calculator import FileIOCalculator, CalculationFailed, equal, compare_atoms, PropertyNotPresent
from xespresso.xio import write_espresso_in, read_espresso_input, read_espresso_asei, write_espresso_asei, get_atomic_species, get_atomic_constraints
from ase.io.espresso import read_espresso_in
import os
import shutil
import numpy as np
from datetime import datetime
import copy
import logging
import sys
import subprocess
import warnings
from pprint import pprint

logging.basicConfig(stream=sys.stdout,
                    format=('%(levelname)-8s '
                            '[%(funcName)-20s]: %(message)s'),
                    level=logging.INFO)

logger = logging.getLogger(__name__)


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'


class Espresso(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms', 'time']
    command = 'PACKAGE.x  PARALLEL  -in  PREFIX.PACKAGEi  >  PREFIX.PACKAGEo'
    discard_results_on_any_change = False

    def __init__(self, label='xespresso', prefix=None, atoms=None, package='pw', parallel='',
                 queue=None, debug=False, **kwargs):
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
            Choose the quantum espresso package: pw, dos, projwfc, band, pp, ph, ..
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
              >>> calc.post('package' = 'dos', **kwargs)
              >>> calc.post('package' = 'pp', **kwargs)
            3. non-self-consistent calculation can be made as follows:
              >>> calc.nscf(calculation = 'nscf', kpts=(4, 4, 4))
              >>> calc.nscf_calculate()
        """
        print("{0:=^60}".format(package))
        if debug:
            logger.setLevel(debug)
        self.set_label(label, prefix)
        self.scf_directory = None
        self.scf_parameters = None
        self.scf_results = None
        kwargs = self.check_input(kwargs, prefix=self.prefix)
        self.ase_parameters = kwargs
        FileIOCalculator.__init__(self, restart=self.directory,
                                  label=self.label, atoms=atoms, **kwargs)
        if atoms:
            self.atoms = atoms
        self.queue = queue
        self.parallel = parallel
        self.package = package
        self.parallel = parallel
        self.debug = debug
        # self.discard_results_on_any_change = False

    def set_label(self, label, prefix):
        '''
        set directory and prefix from label
        '''
        self.directory = label
        if not prefix:
            self.prefix = os.path.split(label)[1]
        else:
            self.prefix = prefix
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        self.label = os.path.join(self.directory, self.prefix)
        self.pwi = os.path.join(self.directory, '%s.pwi' % self.prefix)
        self.pwo = os.path.join(self.directory, '%s.pwo' % self.prefix)
        self.asei = os.path.join(self.directory, '%s.asei' % self.prefix)
        self.asei_temp = os.path.join(
            self.directory, '.%s.asei_temp' % self.prefix)
        self.save_directory = os.path.join(
            self.directory, '%s.save' % self.prefix)
        logger.debug('Directory: %s' % (self.directory))
        logger.debug('Prefix: %s' % (self.prefix))

    def check_pseudopotentials(self, pseudopotentials):
        """
        """
        # pseudopotentials
        if 'species' not in self.atoms.arrays:
            all_species = set(self.atoms.get_chemical_symbols())
        else:
            all_species = set(self.atoms.arrays['species'])
        new_pseudopotentials = {}
        for species in all_species:
            assert species in pseudopotentials, '\n  Species: %s, does not have a pseudopotentials!' % (
                species)
            new_pseudopotentials[species] = pseudopotentials[species]

    @classmethod
    def check_input(cls, kwargs, prefix='pw', package='PW'):
        """
        Chekc input paramter type.
        set parameters for pw.x
        all other parameters are stored in 'input_data'
        """
        from xespresso.xio import sort_qe_input, check_qe_input
        sorted_parameters, unuse_parameters = sort_qe_input(kwargs)
        check_qe_input(sorted_parameters['input_data'])
        if unuse_parameters:
            logger.debug(
                'Warnning: parameter %s not used, please check.' % unuse_parameters)
        sorted_parameters['input_data']['CONTROL']['prefix'] = prefix
        sorted_parameters['input_data']['CONTROL']['verbosity'] = 'high'
        # pprint(sorted_parameters)
        return sorted_parameters

    def read(self, label):
        """Read the files in a calculation if they exist.
        This function reads self.parameters and atoms.
        """
        # Else read a regular calculation. we start with reading stuff
        # that is independent of the calculation state.
        self.directory = label
        self.restart_atoms = None
        self.restart_parameters = {}
        if os.path.exists(label) and os.path.exists(self.asei):
            logger.debug('Try to read information from folder: %s.' % (label))
            try:
                # atoms, parameters = self.read_xml_file()
                atoms, parameters = read_espresso_asei(self.asei)
                self.restart_atoms = atoms
                self.restart_parameters = copy.deepcopy(parameters)
                # self.parameters = copy.deepcopy(parameters)
                self.read_results()
            except Exception as e:
                logger.debug(
                    'Directory: %s exist, faild to read. %s' % (label, e))
        else:
            logger.debug(
                'Directory %s is not a espresso folder, start a new calculation:' % (label))

    def set(self, **kwargs):
        """
        """
        from xespresso.input_parameters import restart_ignore
        self.changed_parameters, self.igonre_parameters = self.compare_parameters(
            self.restart_parameters, self.ase_parameters, ignore=restart_ignore['PW'])
        self.parameters = kwargs
        if self.discard_results_on_any_change and self.changed_parameters:
            self.reset()
        return self.changed_parameters

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
                    logger.debug('Atoms changed, run a new calculation')
                    return True
        else:
            if compare_atoms(self.restart_atoms, atoms, tol=tol,
                             excluded_properties=set(self.ignored_changes)):
                logger.debug('Atoms changed, run a new calculation')
                return True
        # if self.igonre_parameters:
            # logger.debug('Parameter: %s are ignored, results are not affected.'%self.igonre_parameters)
        if self.changed_parameters:
            logger.debug('Parameter: %s changed, results are affected' %
                         self.changed_parameters)
            return True
        logger.debug('Same geometry and parameters, use previous results.')
        converged, meg = self.read_convergence()
        if converged > 0:
            logger.debug('Not converged: %s' % meg)
            return True
        return False

    def compare_parameters(self, para1, para2, ignore=[]):
        """
        """
        from xespresso.input_parameters import qe_namespace, default_parameters, restart_ignore

        changed_parameters = []
        igonre_parameters = []
        if not para1:
            changed_parameters = ['all']
            return changed_parameters, igonre_parameters
        default_parameters = default_parameters['PW']
        # pseudopotentials
        key = 'pseudopotentials'
        try:
            for species, value in para1[key].items():
                if value != para2[key][species]:
                    changed_parameters.append(key)
                    continue
        except:
            changed_parameters.append(key)
        # kpts
        key = 'kpts'
        try:
            if para1[key] != para2[key]:
                changed_parameters.append(key)
        except:
            changed_parameters.append(key)
        # input_data
        for section, paras in para1['input_data'].items():
            if section == 'INPUT_NTYP':
                changed_parameters1, igonre_parameters1 = self.compare_dict(
                    para1['input_data'][section], para2['input_data'][section], restart_ignore['PW'])
            else:
                changed_parameters1, igonre_parameters1 = self.compare_dict(
                    para1['input_data'][section], para2['input_data'][section], restart_ignore['PW'], default=default_parameters[section])
            changed_parameters.extend(changed_parameters1)
            igonre_parameters.extend(igonre_parameters1)
        return changed_parameters, igonre_parameters

    def compare_dict(self, dict1, dict2, ignore=[], default=None):
        igonre_parameters = []
        changed_parameters = []
        keys = set(list(dict1.keys()) + list(dict2.keys()))
        print
        for key in keys:
            if key in ignore:
                igonre_parameters.append(key)
            elif key not in dict1:
                changed_parameters.append(key)
            elif key not in dict2:
                if default and self.compare_value(dict1[key], default[key]):
                    continue
                changed_parameters.append(key)
            elif not self.compare_value(dict1[key], dict2[key]):
                changed_parameters.append(key)
        return changed_parameters, igonre_parameters

    def compare_value(self, v1, v2, tol=1e-5):
        """
        """
        if isinstance(v1, str):
            if v1.upper() != v2.upper():
                return False
        elif isinstance(v1, dict):
            changed_parameters, igonre_parameters = self.compare_dict(v1, v2)
            if changed_parameters:
                return False
            for key, value in v1.items():
                if not self.compare_value(value, v2[key]):
                    return False
        elif isinstance(v1, bool):
            if v1 != v2:
                return False
        else:
            if abs(v1 - v2) > tol:
                return False
        return True

    def read_results(self):
        '''
        get atomic species
        '''
        pwo = os.path.join(self.directory, '%s.pwo' % self.prefix)
        pwi = os.path.join(self.directory, '%s.pwi' % self.prefix)
        convergence, meg = self.read_convergence()
        if convergence != 0:
            logger.debug('Not converged. %s' % (meg))
        try:
            atoms = read_espresso_in(pwi)
            self.results['atoms'] = atoms
        except Exception as e:
            logger.debug('Read input: %s, failed! %s' % (pwi, e))
        try:
            output = io.read(pwo)
            atomic_species = get_atomic_species(pwo)
            constraints = get_atomic_constraints(pwo, len(output))
            output.set_constraint(None)
            output.set_constraint(constraints)
            if atomic_species:
                output.new_array('species', np.array(
                    atomic_species, dtype='U20'))
            self.calc = output.calc
            self.results = output.calc.results
            self.results['atoms'] = output
            self.efermi = self.get_fermi_level()
            # self.nspins = self.get_number_of_spins()
            logger.debug('Read result successfully!')
        except Exception as e:
            logger.debug('Read output: %s, failed! %s' % (pwo, e))
        self.results['convergence'] = convergence
        # logger.debug('Read result failed!')
        # pwos = [file for file in os.listdir(self.directory) if pwo in file]
        # output = None
        # for pwo in pwos:
        # atomic_species = None
        # pwo = os.path.join(self.directory, pwo)
        # atomic_species = get_atomic_species(pwo)

    def write_input(self, atoms, properties=None, system_changes=None):
        from xespresso.scheduler import set_queue
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        set_queue(self)
        write_espresso_in(self.label + '.pwi', atoms, **self.parameters)
        write_espresso_asei(self.label + '.asei', atoms, self.parameters)
        logger.debug('Write input successfully')

    def check_xml_file(self):
        '''
        If data-file-schema.xml is readable, then xml_end = True
        '''
        self.save_directory = os.path.join(
            self.directory, '%s.save' % self.prefix)
        if not os.path.exists(self.save_directory):
            return False
        xml0 = os.path.join(self.save_directory, 'data-file-schema.xml')
        xmls = [file for file in os.listdir(
            self.save_directory) if 'data-file-schema.xml' in file]
        goodxml = []
        for xml in xmls:
            xml = os.path.join(self.save_directory, xml)
            with open(xml, 'r') as f:
                lines = f.readlines()
                if len(lines) < 1:
                    continue
                if '</qes:espresso>' in lines[-1]:
                    goodxml.append(xml)
        xml_end = False
        nxml = len(goodxml)
        # print(nxml, goodxml, xml0)
        if nxml > 0:
            if xml0 not in goodxml:
                shutil.copy(goodxml[-1], xml0)
            xml_end = True
        return xml_end

    def read_xml_file(self):
        from xespresso.utils.xml_parser import xml_parser
        if not self.check_xml_file():
            logger.debug('xml file not finished.')
        #
        xmlfile = os.path.join(self.save_directory, 'data-file-schema.xml')
        if not os.path.exists(xmlfile):
            return None, {}
        try:
            atoms, xml_input = xml_parser(xmlfile)
            logger.debug(
                'Read geometry and parameters from xml file successfully.')
        except:
            logger.debug('Read xml file failed.')
            return None, {}
        return atoms, xml_input

    def read_convergence(self):
        '''
        Read the status of the calculation.
        {
        '0': 'Done',
        '-1': ' manual cancelled',
        '1': 'Maximum CPU time exceeded', convergence NOT achieved after n iterations
             Then, change the parameter: 1) mixing_beta or 2) degauss
        '2': 'manual restart', 'Pending or not submit',
        '3': 'other errors',
        '4': 'unknow error',
        }
        '''
        # read the error file from queue
        # if self.queue:
        #     errfile = self.label + '.err'
        #     if not os.path.exists(errfile):
        #         logger.debug('%s not exists'%errfile)
        #         # return 2, 'Pending or not submit'
        #     else:
        #         errs = ['RESTART', 'pw.x', 'out-of-memory',
        #                 'NODE FAILURE', 'TIME LIMIT', 'COMMUNICATION',
        #                 'segmentation fault', 'PARSE_ERR', 'mpirun']
        #         cancelled = False
        #         with open(errfile, 'r') as f:
        #             lines = f.readlines()
        #             for line in lines:
        #                 if 'CANCELLED' in line:
        #                     cancelled = True
        #                 for err in errs:
        #                     if err in line:
        #                         logger.debug('Need restart')
        #                         return 3, line
        #             # cancelled by owner
        #             if len(lines) > 0:
        #                 if cancelled:
        #                     logger.debug('Cancelled')
        #                     return -1, lines[0]
        # no error from queue, or job not run in the queue
        output = self.label + '.pwo'
        if not os.path.exists(output):
            # print('%s not exists'%output)
            return 3, 'No pwo output file'
        with open(output, 'r') as f:
            lines = f.readlines()
            if len(lines) == 0:
                return 1, 'pwo file has nothing'
            stime = lines[1].split('starts on')[1]
            nlines = len(lines)
            n = min([200, nlines])
            lastlines = lines[-n:-1]
            for line in lastlines:
                if line.rfind('too many bands are not converged') > -1:
                    logger.debug('Need restart')
                    return 1, 'Reason: %s' % (line)
                if line.rfind('convergence NOT achieved after') > -1:
                    logger.debug('Need restart: %s' % line)
                    return 1, 'Reason: %s' % (line)
                if line.rfind('Maximum CPU time exceeded') > -1:
                    logger.debug('Need restart')
                    return 2, 'Reason: %s' % (line)
                if line.rfind('JOB DONE.') > -1:
                    logger.debug('JOB DONE')
                    return 0, line
            logger.debug('Not converged, %s' % lastlines[-1])
        return 4, line

    def run(self, atoms=None, restart=False):
        '''
        run and restart
        '''
        if atoms is not None:
            self.atoms = atoms
        icount = 1
        converged = 4
        meg0 = 'None'
        self.atoms.calc = self
        if not restart:
            try:
                logger.debug('Run step: %s' % icount)
                self.atoms.get_potential_energy()
                converged, meg0 = self.read_convergence()
            except Exception as e:
                logger.debug('Not converge: %s' % e)
                converged, meg0 = self.read_convergence()
                icount += 1
                if converged == 0:
                    converged = 4
        # print(converged, xml_end, meg0)
        # restart 0, 1, 2
        exitfile = os.path.join(self.directory, 'EXIT')
        while converged > 0 or restart == 2:
            if os.path.exists(exitfile):
                logger.debug('Exit file exists!')
                os.remove(exitfile)
                sys.exit()
            restart = 1
            self.backup_file('%s.pwo' % self.prefix, directory=self.directory)
            self.backup_file('data-file-schema.xml',
                             directory=self.save_directory)
            logger.debug('Run step: %s' % icount)
            xml_end = self.check_xml_file()
            if xml_end:
                self.parameters['input_data']['CONTROL']['restart_mode'] = 'restart'
            else:
                self.parameters['input_data']['CONTROL']['restart_mode'] = 'from_scratch'
            if converged == 1:
                if 'electron_maxstep' in self.parameters['input_data']['ELECTRONS']:
                    electron_maxstep = int(
                        self.parameters['input_data']['ELECTRONS']['electron_maxstep'])
                else:
                    electron_maxstep = 100
                if 'mixing_beta' in self.parameters['input_data']['ELECTRONS']:
                    mixing_beta = float(
                        self.parameters['input_data']['ELECTRONS']['mixing_beta'])
                else:
                    mixing_beta = 0.7
                if electron_maxstep < 100:
                    self.parameters['input_data']['ELECTRONS']['electron_maxstep'] = electron_maxstep + 50
                    logger.debug('electron_maxstep change from %s to %s' % (
                        electron_maxstep, electron_maxstep + 50))
                self.parameters['input_data']['ELECTRONS']['mixing_beta'] = mixing_beta/2.0
                logger.debug('mixing_beta change from %s to %s' %
                             (mixing_beta, mixing_beta/2.0))
            try:
                icount += 1
                atoms.get_potential_energy()
                converged, meg = self.read_convergence()
            except Exception as e:
                logger.debug('Not converge: %s' % e)
                converged, meg = self.read_convergence()
                meg = str(e)
                restart = 2
                # if converged == 0:
                converged = 4
            logger.debug('%s, %s' % (converged, meg))
            if meg == meg0 and converged != 2:
                logger.debug(
                    'Exit! Maybe deadblock! Same error message! %s' % meg)
                sys.exit()
            meg0 = meg

    def backup_file(self, src, directory='.'):
        '''
        compare files
        backup
        '''
        import filecmp
        tnow = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        new_src = os.path.join(directory, '%s-%s' % (tnow, src))
        src = os.path.join(directory, src)
        if not os.path.exists(src):
            return 0
        ftype = src.split('.')[-1]
        flag = True
        for dst in os.listdir(directory):
            dst = os.path.join(directory, dst)
            if dst == src:
                continue
            if dst.endswith(ftype):
                if filecmp.cmp(dst, src):
                    flag = False
        # print('backup files: ', flag, src, new_src)
        if flag:
            shutil.copy(src, new_src)
        return 0

    def plot_phdos(self, fldos=None, ax=None, output=None):
        '''
        '''
        import matplotlib.pyplot as plt
        if fldos is None:
            fldos = self.prefix
        phdos = np.loadtxt(self.directory+'/%s.phdos' % fldos)
        self.phdos = phdos
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 3))
            # ax = plt.gca()
        ax.plot(self.phdos[:, 0], self.phdos[:, 1], linewidth=0.7)
        ax.set_xlabel('Frequency (cm^-1)')
        ax.set_ylabel('DOS (a.u.)')
        if output is not None:
            plt.savefig('%s' % output)
        return ax

    def get_work_function(self, ax=None, inpfile='potential.cube', output=None, shift=False):
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
            plt.savefig('%s' % output)
        # plt.show()
        atoms = self.results['atoms']
        pos = max(atoms.positions[:, 2] + 3)
        ind = (xaxis > pos) & (xaxis < pos + 3)
        wf = np.average(axy[ind]) - ef
        print('min: %s, max: %s' % (pos, pos + 3))
        print('The workfunction is {0:1.2f} eV'.format(wf))

    def get_bader_charge(self, inpfile=None):
        '''
        '''
        from ase.io.bader import attach_charges
        if not inpfile:
            inpfile = '%s.cube' % self.prefix
        command = 'bader %s' % inpfile
        print(command)
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err
        acf = os.path.join(self.directory, 'ACF.dat')
        attach_charges(self.results['atoms'], acf)

    def read_time(self, ):
        '''
        '''
        pwo = os.path.join(self.directory, self.prefix + '.pwo')
        with open(pwo, 'r') as f:
            ts = 0
            lines = f.readlines()
            for line in lines[::-1]:
                if 'PWSCF' in line and 'WALL' in line:
                    t = line.split('CPU')[-1].split('WALL')[0]
                    if 's' in t:
                        ts = float(t.split('m')[-1].split('s')[0])
                    if 'm' in t:
                        tm = float(t.split('h')[-1].split('m')[0])
                        ts += tm*60
                    if 'h' in t:
                        th = float(t.split('h')[0])
                        ts += th*3600
                    break
                self.results['time'] = ts
        return ts

    def get_time(self, ):
        t = 0
        filename = os.path.join(self.directory, '%s.pwo' % self.prefix)
        with open(filename) as f:
            lines = f.readlines()
            for line in lines[-200:]:
                if 'init_run     :' in line or 'electrons    :' in line:
                    t += float(line.split('CPU')[1].split('s WALL')[0])
        self.results['time'] = t
        return t

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
        # keys = ['wfc']
        # for file in files:
        #     for key in keys:
        #         if key in file:
        #             os.remove(os.path.join(self.save_directory, file))

    def get_fermi_level(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Fermi level')
        return self.calc.get_fermi_level()

    def get_ibz_k_points(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'IBZ k-points')
        ibzkpts = self.calc.get_ibz_k_points()
        if ibzkpts is None:
            warnings.warn(warn_template % 'IBZ k-points')
        return ibzkpts

    def get_k_point_weights(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'K-point weights')
        k_point_weights = self.calc.get_k_point_weights()
        if k_point_weights is None:
            warnings.warn(warn_template % 'K-point weights')
        return k_point_weights

    def get_eigenvalues(self, **kwargs):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Eigenvalues')
        eigenvalues = self.calc.get_eigenvalues(**kwargs)
        if eigenvalues is None:
            warnings.warn(warn_template % 'Eigenvalues')
        return eigenvalues

    def get_number_of_spins(self):
        if self.calc is None:
            raise PropertyNotPresent(error_template % 'Number of spins')
        nspins = self.calc.get_number_of_spins()
        if nspins is None:
            warnings.warn(warn_template % 'Number of spins')
        return nspins
