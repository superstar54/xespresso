import numpy as np
from xespresso import Espresso
import matplotlib.pyplot as plt
import os
import shutil
import subprocess


class COHP:
    def __init__(self,
                 directory='.',
                 prefix='cohp',
                 indexs=[[1, 2]],
                 queue=False,
                 command='lobster',
                 **kwargs):
        """
        Energetic window [COHPstartEnergy, COHPendEnergy]

        Args:
            command (str, optional): Path to lobster executable file.
        """
        self.directory = directory
        self.prefix = prefix
        self.parameters = kwargs
        self.indexs = indexs
        self.command = command

    def run(self, cpu=1):
        """Run COHP analysis using lobster

        Args:
            cpu (int, optional): CPU threads used. Defaults to 1.

        Raises:
            EnvironmentError: [description]
        """
        shutil.copy(os.path.join(self.directory, '%s.pwi' % self.prefix),
                    os.path.join(self.directory, '%s.scf.in' % self.prefix))
        self.write_input()
        command = self.command
        os.system('export OMP_NUM_THREADS=%s' % cpu)
        print('Running %s on %s cpus' % (command, cpu))
        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Command "{}" failed in '
                   '{} with error code {}'.format(command, path, errorcode))
            print(msg)
            exit()
        print('Done: %s' % command)
        #

    def write_input(self, ):
        """make input file for lobster
        """
        filename = os.path.join(self.directory, 'lobsterin')
        with open(filename, 'w') as f:
            for key, value in self.parameters.items():
                if isinstance(value, list):
                    for subvalue in value:
                        f.write('  %s  %s \n' % (key, subvalue))
                else:
                    f.write('  %s  %s \n' % (key, value))
            for ind in self.indexs:
                f.write('cohpbetween atom %s atom %s \n' % (ind[0], ind[1]))

    def read_cohp(self, ):
        """read COHP anlysis results
        """
        from ase.calculators.vasp import VaspDos
        dos = VaspDos(os.path.join(self.directory, 'DOSCAR.lobster'))
        self.dos = dos._total_dos[1]
        self.dos_energies = dos._total_dos[0]
        # read COHP results
        datas = np.genfromtxt(os.path.join(self.directory, 'COHPCAR.lobster'),
                              skip_header=3 + len(self.indexs))
        self.cohp = datas[:, 1]
        self.cohp_energies = datas[:, 0]
        # read COOP results
        datas = np.genfromtxt(os.path.join(self.directory, 'COOPCAR.lobster'),
                              skip_header=3 + len(self.indexs))
        self.coop = datas[:, 1]
        self.coop_energies = datas[:, 0]
        # read COBI results
        datas = np.genfromtxt(os.path.join(self.directory, 'COBICAR.lobster'),
                              skip_header=3 + len(self.indexs))
        self.cobi = datas[:, 1]
        self.cobi_energies = datas[:, 0]

    def read_icohp(self, ):
        """read icohp data for bond strength evaluation
        """
        from ase.calculators.vasp import VaspDos
        import pandas as pd
        # read icohp
        datas = np.genfromtxt(os.path.join(self.directory,
                                           'ICOHPLIST.lobster'),
                                           dtype = [int, 'S5','S5',float,float], 
                                           usecols=(0,1,2,3,7),
                                           names=['interaction','atom1','atom2','distance','ICOHP'],
                                           skip_header=1)
        self.icohp = pd.DataFrame(datas)
        # self.icohp_energies = datas[:, 7]
        # read icoop
        datas = np.genfromtxt(os.path.join(self.directory,
                                           'ICOOPLIST.lobster'),
                                           dtype = [int, 'S5','S5',float,float], 
                                           usecols=(0,1,2,3,7),
                                           names=['interaction','atom1','atom2','distance','ICOOP'],
                                           skip_header=1)
        self.icoop = pd.DataFrame(datas)
        # self.icoop_energies = datas[:, 7]
        # read icobi
        datas = np.genfromtxt(os.path.join(self.directory,
                                           'ICOHPLIST.lobster'),
                                           dtype = [int,'S5','S5',float,float], 
                                           usecols=(0,1,2,3,7),
                                           names=['interaction','atom1','atom2','distance','ICOBI'],
                                           skip_header=1)
        self.icobi = pd.DataFrame(datas)
        # self.icobi_energies = datas[:, 7]

    def plot_cohp(self, ax=None, output=None):
        """create a DOS and COHP plot

        Args:
            ax (COHP object, optional):  Defaults to None.
            output (strng, optional): filename if want to save plot. Defaults to None.

        Returns:
            plot: a DOS with COHP plot
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.dos, self.dos_energies)
        ax[0].set_ylabel('E-Ef (eV)')
        ax[0].axhline(y=0, color='k', linewidth = 0.5, linestyle = '--')
        ax[1].plot(-self.cohp, self.cohp_energies)
        ax[1].axvline(x=0, color='k', linewidth = 0.2)
        ax[1].axhline(y=0, color='k', linewidth = 0.2, linestyle = '--')
        ax[0].set_xlabel('DOS')
        ax[1].set_xlabel('-pCOHP')
        if output:
            plt.savefig(output)
        return ax

    def plot_coop(self, ax=None, output=None):
        """create a DOS with COOP plot

        Args:
            ax (COHP object, optional): Defaults to None.
            output (string, optional): filename if want to save plot. Defaults to None.

        Returns:
            plot: a DOS with COOP plot
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.dos, self.dos_energies)
        ax[0].set_ylabel('E-Ef (eV)')
        ax[0].axhline(y=0, color='k', linewidth = 0.5, linestyle = '--')
        ax[1].plot(self.coop, self.coop_energies)
        ax[1].axvline(x=0, color='k', linewidth = 0.5)
        ax[1].axhline(y=0, color='k', linewidth = 0.5, linestyle = '--')
        ax[0].set_xlabel('DOS')
        ax[1].set_xlabel('pCOOP')
        if output:
            plt.savefig(output)
        return ax

    def plot_cobi(self, ax=None, output=None):
        """create a DOS with COBI plot

        Args:
            ax (COHP object, optional): Defaults to None.
            output (string, optional): filename if want to save plot. Defaults to None.

        Returns:
            plot: a DOS with COBI plot
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(1, 2)
        ax[0].plot(self.dos, self.dos_energies)
        ax[0].set_ylabel('E-Ef (eV)')
        ax[0].axhline(y=0, color='k', linewidth = 0.5, linestyle = '--')
        ax[1].plot(self.cobi, self.cobi_energies)
        ax[1].axvline(x=0, color='k', linewidth = 0.5)
        ax[1].axhline(y=0, color='k', linewidth = 0.5, linestyle = '--')
        ax[0].set_xlabel('DOS')
        ax[1].set_xlabel('COBI')
        if output:
            plt.savefig(output)
        return ax
