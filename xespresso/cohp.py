import numpy as np
from xespresso import Espresso
import matplotlib.pyplot as plt
import os
import shutil
import subprocess

class COHP:
    def __init__(self, directory = '.', prefix = 'cohp', indexs = [[1, 2]], queue = False, command = 'lobster', **kwargs):
        """
        Energetic window [COHPstartEnergy, COHPendEnergy]
        """
        self.directory = directory
        self.prefix = prefix
        self.parameters = kwargs
        self.indexs = indexs
        self.command = command
    def run(self, cpu = 1):
        '''
        '''
        shutil.copy(os.path.join(self.directory, '%s.pwi'%self.prefix), os.path.join(self.directory, '%s.scf.in'%self.prefix))
        self.write_input()
        command = self.command
        os.system('export OMP_NUM_THREADS=%s'%cpu)
        print('Running %s on %s cpus' %(command, cpu))
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
            print(msg)
            exit()
        print('Done: %s' %command)
        #
    def write_input(self, ):
        filename = os.path.join(self.directory, 'lobsterin')
        with open(filename, 'w') as f:
            for key, value in self.parameters.items():
                if isinstance(value, list):
                    for subvalue in value:
                        f.write('  %s  %s \n' %(key, subvalue))
                else:
                    f.write('  %s  %s \n' %(key, value))
            for ind in self.indexs:
                f.write('cohpbetween atom %s atom %s \n' %(ind[0], ind[1]))
    def read_cohp(self, ):
        from ase.calculators.vasp import VaspDos
        dos = VaspDos(os.path.join(self.directory, 'DOSCAR.lobster'))
        self.dos = dos._total_dos[1]
        self.dos_energies = dos._total_dos[0]
        #
        datas = np.genfromtxt(os.path.join(self.directory, 'COHPCAR.lobster'), skip_header = 3 + len(self.indexs))
        self.cohp = datas[:,1]
        self.cohp_energies = datas[:,0]
        #
        datas = np.genfromtxt(os.path.join(self.directory, 'COOPCAR.lobster'), skip_header = 3 + len(self.indexs))
        self.coop = datas[:,1]
        self.coop_energies = datas[:,0]
    def read_icohp(self, ):
        from ase.calculators.vasp import VaspDos
        datas = np.genfromtxt(os.path.join(self.directory, 'ICOOPLIST.lobster'), skip_header = 1 + len(self.indexs))
        self.icohp = datas
        # self.icohp_energies = datas[:, 7]
        #
        datas = np.genfromtxt(os.path.join(self.directory, 'ICOOPLIST.lobster'), skip_header = 1 + len(self.indexs))
        self.icoop = datas
        # self.icoop_energies = datas[:, 7]
    def plot_cohp(self, ax = None, output = None):
        '''
        '''
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.dos, self.dos_energies)
        ax[0].set_ylabel('Energy (eV)')
        ax[1].plot(-self.cohp, self.cohp_energies)
        ax[1].axvline(x = 0, color = 'k')
        ax[0].set_xlabel('DOS')
        ax[1].set_xlabel('-COHP')
        if output:
            plt.savefig(output)
        return ax
    def plot_coop(self, ax = None, output = None):
        '''
        '''
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.dos, self.dos_energies)
        ax[0].set_ylabel('Energy (eV)')
        ax[1].plot(self.coop, self.coop_energies)
        ax[1].axvline(x = 0, color = 'k')
        ax[0].set_xlabel('DOS')
        ax[1].set_xlabel('COOP')
        if output:
            plt.savefig(output)
        return ax
    
    