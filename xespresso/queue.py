import os


class Queue():
    def __init__(self, package = None, parallel = None, queue = None, command = None, type = 'slurm') -> None:
        if command is None:
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
        pass

    def set_queue(self):
        '''
        If queue, change command to "sbatch .job_file".
        The queue information are written into '.job_file'
        '''
        
        if self.queue:
            script = ''
            if 'config' in self.queue:
                cf = os.path.join(os.environ['HOME'], self.queue['config'])
                file = open(cf, 'r')
                script = file.read()
                # del self.queue['config']
            else:
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
                for key, value in self.queue.items():
                    if key == 'config': continue
                    if value:
                        fh.writelines("#SBATCH --%s=%s\n" %(key, value))
                fh.writelines("%s \n"%script)
                fh.writelines("%s \n" % command)
            self.command = "sbatch {0}".format('.job_file')
        else:
            self.command = command
        self.logger.debug('Command: %s'%(self.command))