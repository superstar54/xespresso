"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_COMMAND="/path/to/package.x -in PREFIX.pwi > PREFIX.pwo"

Run pw.x jobs.
"""

import numpy as np
import warnings
from ase import io
from xespresso import Espresso
from xespresso.xio import write_neb_in


error_template = 'Property "%s" not available. Please try running Quantum\n' \
                 'Espresso first by calling Atoms.get_potential_energy().'

warn_template = 'Property "%s" is None. Typically, this is because the ' \
                'required information has not been printed by Quantum ' \
                'Espresso at a "low" verbosity level (the default). ' \
                'Please try running Quantum Espresso with "high" verbosity.'
# 
    
class NEBEspresso(Espresso):
    """
    """
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms', 'time']
    command = 'pw.x -in PREFIX.pwi > PREFIX.pwo'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='xespresso', images=None, climbing_images = False, path_data = None, package = 'neb', parallel = '',
                 queue = None,
                 **kwargs):
        """

        """
        atoms = images[0]
        Espresso.__init__(self, restart, ignore_bad_restart_file,
                 label, atoms, package, parallel,
                 queue, **kwargs)
        self.images = images
        self.climbing_images = climbing_images
        self.path_data = path_data
    def write_input(self, images, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, self.images, properties, system_changes)
        write_neb_in(self.label + '.nebi', self.images, self.climbing_images, self.path_data, **self.parameters)
    def read_results(self):
        try:
            paths, energies, error = self.read_dat()
            self.paths = paths
            self.energies = energies
            self.error = error
            interpolation_paths, interpolation_energies = self.read_int()
            self.interpolation_paths = interpolation_paths
            self.interpolation_energies = interpolation_energies
        except Exception as e:
            print('\nread paths and energies failed\n')
            print(e)
        try:
            images = self.read_xyz()
            self.images = images
        except Exception as e:
            print('\nread xyz failed\n')
            print(e)
    def read_dat(self):
        data = np.loadtxt(self.directory+'/%s.dat' % self.prefix)
        paths, energies, error = data[:, 0], data[:, 1], data[:, 2]
        return paths, energies, error
    def read_int(self):
        data = np.loadtxt(self.directory+'/%s.int' % self.prefix)
        paths, energies = data[:, 0], data[:, 1]
        return paths, energies
    def read_xyz(self):
        images = io.read(self.directory+'/%s.xyz' % self.prefix)
        return images
    def read_path(self):
        pass
    def plot(self):
        import matplotlib.pyplot as plt
        plt.plot(self.paths, self.energies, 'o')
        plt.plot(self.interpolation_paths, self.interpolation_energies, '--')
        plt.xlabel('Path')
        plt.ylabel('Energy (eV)')
        plt.show()