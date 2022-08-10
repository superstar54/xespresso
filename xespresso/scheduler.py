import os
import logging
logger = logging.getLogger(__name__)


config_files = [os.path.join(os.environ['HOME'], '.xespressorc'),
            '.xespressorc']
            
def set_queue(calc, package = None, parallel = None, queue = None, command = None):
    '''
    If queue, change command to "sbatch .job_file".
    The queue information are written into '.job_file'
    '''
    if queue is None:
        queue = calc.queue
    else:
        calc.queue = queue
    if package is None:
        package = calc.package
    else:
        calc.package = package
    if parallel is None:
        parallel = calc.parallel
    else:
        calc.parallel = parallel
    if command is None:
        command = os.environ.get('ASE_ESPRESSO_COMMAND')
    if 'PACKAGE' in command:
        if 'pw' in package:
            command = command.replace('PACKAGE', package, 1)
            command = command.replace('PACKAGE', 'pw', 2)
        else:
            command = command.replace('PACKAGE', package)
    if 'PREFIX' in command:
        command = command.replace('PREFIX', calc.prefix)
    if 'PARALLEL' in command:
        command = command.replace('PARALLEL', parallel)
    logger.debug('Espresso command: %s'%(command))
    if queue:
        script = ''
        if 'config' in queue:
            cf = os.path.join(os.environ['HOME'], queue['config'])
            file = open(cf, 'r')
            script = file.read()
            # del queue['config']
        else:
            for cf in config_files:
                if os.path.exists(cf):
                    file = open(cf, 'r')
                    script = file.read()
        jobname = calc.prefix
        with open('%s/.job_file' % calc.directory, 'w') as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s \n" % jobname)
            fh.writelines("#SBATCH --output=%s.out\n" % calc.prefix)
            fh.writelines("#SBATCH --error=%s.err\n" % calc.prefix)
            fh.writelines("#SBATCH --wait\n")
            for key, value in queue.items():
                if key == 'config': continue
                if value:
                    fh.writelines("#SBATCH --%s=%s\n" %(key, value))
            fh.writelines("%s \n"%script)
            fh.writelines("%s \n" % command)
        calc.command = "sbatch {0}".format('.job_file')
    else:
        calc.command = command
    logger.debug('Queue command: %s'%(calc.command))