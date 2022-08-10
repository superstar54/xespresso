
from ase.calculators.calculator import FileIOCalculator, CalculationFailed, equal, compare_atoms, PropertyNotPresent
from xespresso.xio import read_espresso_asei
from xespresso.xespresso import Espresso
import os
import logging
logger = logging.getLogger(__name__)


class EspressoNscf(FileIOCalculator):

    package='pw'

    def __init__(self, scf_directory, prefix, parallel='',
            queue=None, debug=False, 
            kpts = (10, 10, 10), 
            **kwargs):
        print("{0:=^60}".format('nscf'))
        self.scf_directory = scf_directory
        self.prefix = prefix
        self.parallel = parallel
        self.queue = queue
        if debug:
            logger.setLevel(debug)
        self.load_scf()
        self.parameters['kpts'] = kpts
        self.parameters.update(kwargs)
        self.parameters = Espresso.check_input(self.parameters, self.prefix)
        #
        self.directory = os.path.join(self.scf_directory, 'nscf')
        self.label = os.path.join(self.directory, self.prefix)
        self.asei = os.path.join(self.directory, '%s.nscf_asei'%self.prefix)
        # create working directory
        self.set_label(self.directory, self.prefix)
    
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
        self.pwi = os.path.join(self.directory, '%s.pwi'%self.prefix)
        self.pwo = os.path.join(self.directory, '%s.pwo'%self.prefix)
        self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
        logger.debug('Directory: %s'%(self.directory))
        logger.debug('Prefix: %s'%(self.prefix))

    def load_scf(self):
        """load parameters and results from privous scf calculation
        """
        import copy
        asei = os.path.join(self.scf_directory, '%s.asei'%self.prefix)
        self.atoms, self.scf_parameters = read_espresso_asei(asei, 'PW')
        self.parameters = copy.deepcopy(self.scf_parameters)
        self.parameters['input_data']['CONTROL']['calculation'] = 'nscf'
        self.parameters['input_data']['CONTROL']['verbosity'] = 'high'
        self.parameters['input_data']['CONTROL']['outdir'] = '../'
        self.parameters['input_data']['SYSTEM']['occupations'] = 'tetrahedra'
        self.parameters['input_data']['SYSTEM'].pop('degauss', None)

    def write_input(self, atoms, properties=None, system_changes=None):
        from xespresso.xio import write_espresso_asei, write_espresso_in
        from xespresso.scheduler import set_queue
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        set_queue(self, package = self.package, parallel = self.parallel,
                queue = self.queue)
        write_espresso_in(self.label + '.pwi', atoms, **self.parameters)
        logger.debug("write asei file: {}".format(self.asei))
        logger.debug("  with parameters: {}".format(self.parameters))
        write_espresso_asei(self.asei, self.state_info, self.parameters)
        logger.debug('Write input successfully')

    @property
    def state_info(self):
        from xespresso.utils import get_hash
        chargeFile = os.path.join(self.scf_directory, '%s.save/charge-density.dat'%self.prefix)
        state_info = get_hash(chargeFile)
        return state_info
    
    def check_state(self, ):
        # read information of the charge-density file
        self.state_parameters = self.parameters
        output, meg = self.read_convergence_post('pw')
        logger.debug("Check state: {}".format(meg))
        if output:
            if os.path.isfile(self.asei):
                system_changes =  self.check_state_post(self.asei, package='PW')
                if not system_changes:
                    logger.debug('Use previous results!')
                    return False
            else:
                logger.debug('No asei output. Start a new calculation.')

        else:
            logger.debug('No pw output. Start a new calculation.')

        return True
    
    def check_state_post(self, asei, package):
        logger.debug("asei file: {}".format(asei))
        old_state_info, old_state_parameters = read_espresso_asei(asei, package)
        if not self.state_info == old_state_info:
            logger.debug('File in save changed')
            return True
        elif not self.state_parameters == old_state_parameters:
            logger.debug(self.state_parameters)
            logger.debug(old_state_parameters)
            logger.debug('Parameters changed')
            return True
        else:
            return False

    def read_convergence_post(self, package = 'pw'):
        '''
        Read the status of the calculation.
        {
        '0': 'Done',
        }
        '''
        output = self.label + '.%so'%package
        logger.debug("read output: {}".format(output))
        if not os.path.exists(output):
            # print('%s not exists'%output)
            return False, 'No pwo output file'
        with open(output, 'r') as f:
            lines = f.readlines()
            if len(lines) == 0: 
                return False, 'pwo file has nothing'
            nlines = len(lines)
            n = min([100, nlines])
            for line in lines[-n:-1]:
                if line.rfind('JOB DONE.') > -1:
                    logger.debug('JOB DONE.')
                    return True, line
        return False, line
    
    def run(self):
        if self.check_state():
            self.calculate()