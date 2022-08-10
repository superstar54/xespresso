import os
from ase.calculators.calculator import FileIOCalculator, CalculationFailed
from xespresso.xio import read_espresso_asei, write_espresso_asei
import copy
import logging
logger = logging.getLogger(__name__)


class PostCalculation:
    
    package = 'dos'
    package_parameters = {}

    def __init__(self, parent_directory, prefix, queue = False, parallel = '', **kwargs) -> None:
        self.parent_directory = parent_directory
        self.prefix = prefix
        self.queue = queue
        self.parallel = parallel
        self.parameters = kwargs

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
        self.asei = os.path.join(self.directory, '%s.asei'%self.prefix)
        self.asei_temp = os.path.join(self.directory, '.%s.asei_temp'%self.prefix)
        self.post_asei = os.path.join(self.directory, '%s.post_asei'%self.prefix)
        self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
        logger.debug('Directory: %s'%(self.directory))
        logger.debug('Prefix: %s'%(self.prefix))

    def run(self):
        '''
        todo: 
        '''
        from xespresso.utils import get_hash
        from xespresso.scheduler import set_queue
        
        print('{0:=^60}'.format(self.package))
        self.directory = os.path.join(self.parent_directory, '%s/'%self.package)
        self.set_label(self.directory, self.prefix)
        self.parameters['prefix'] = self.prefix
        self.parameters['outdir'] = '../'
        self.state_parameters = copy.deepcopy(self.parameters)
        self.state_info = None
        for wfc in ['wfc1', 'wfcdw1']:
            wfcFile = os.path.join(self.parent_directory, '%s.save/%s.dat'%(self.prefix, wfc))
            if os.path.isfile(wfcFile):
                self.state_info = get_hash(wfcFile)
        output, meg = self.read_convergence_post(self.package)
        if output:
            logger.debug('Previous calculation done.')
            if os.path.isfile(self.post_asei):
                system_changes =  self.check_state_post(self.post_asei, self.package)
                if not system_changes:
                    logger.debug('File and Parameters did not change. Use previous results!')
                    return 0
        self.post_write_input(self.package, **self.state_parameters)
        write_espresso_asei(self.post_asei, self.state_info, self.state_parameters)
        set_queue(self, package = self.package, parallel = self.parallel, queue = self.queue)
        self.post_calculate()
        self.post_read_results()

    def check_state_post(self, asei, package):
        old_state_info, old_state_parameters = read_espresso_asei(asei, package)
        if not self.state_info == old_state_info:
            logger.debug('File in save changed')
            return True
        elif not self.state_parameters == old_state_parameters:
            logger.debug('Parameters changed')
            return True
        else:
            return False
    

    def post_write_input(self, package, **kwargs):
        filename = os.path.join(self.directory, '%s.%si' %(self.prefix, package))
        with open(filename, 'w') as f:
            for section, parameters in self.package_parameters.items():
                # logger.debug(section)
                if section != 'LINE':
                    f.write('&%s\n'%section)
                    for key, value in kwargs.items():
                        if key in parameters:
                            if isinstance(value, dict):
                                for subkey, subvalue in value.items():
                                    if isinstance(subvalue, str):
                                        f.write('  %s(%s) = "%s", \n' %(key, subkey, subvalue))
                                    else:
                                        f.write('  %s(%s) = %s, \n' %(key, subkey, subvalue))
                            else:
                                if isinstance(value, str):
                                    f.write('  {0:10s} =  "{1}" \n'.format(key, value))
                                else:
                                    f.write('  {0:10s} =  {1} \n'.format(key, value))
                    f.write('/ \n')
                else:
                    for key, value in kwargs.items():
                        if key in parameters:
                            f.write('  %s \n' %(value))


    def post_calculate(self):
        import subprocess
        
        command = self.command
        print('Running %s'%self.package)
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
        print('Done: %s' %self.package)
    
    def post_read_results(self):
        '''
        '''
        pass

    def read_convergence_post(self, package = 'pw'):
        '''
        Read the status of the calculation.
        {
        '0': 'Done',
        }
        '''
        
        output = self.label + '.%so'%package
        if not os.path.exists(output):
            # print('%s not exists'%output)
            return False, 'No pwo output file'
        with open(output, 'r') as f:
            lines = f.readlines()
            if len(lines) == 0: return False, 'pwo file has nothing'
            nlines = len(lines)
            n = min([100, nlines])
            for line in lines[-n:-1]:
                if line.rfind('JOB DONE.') > -1:
                    logger.debug('JOB DONE.')
                    return True, line
        return False, line
    