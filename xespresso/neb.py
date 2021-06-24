"""
Quantum ESPRESSO Calculator for NEB
"""
import os
import numpy as np
from ase import Atoms
from ase import io
from xespresso import Espresso
from xespresso.xio import write_neb_in, write_espresso_asei
from ase.calculators.calculator import FileIOCalculator, compare_atoms
from ase.utils.forcecurve import fit_raw
from ase.units import create_units
units = create_units('2006')


    
class NEBEspresso(Espresso):
    """
    """
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms', 'time', 'neb']
    command = 'neb.x  PARALLEL  -in  PREFIX.nebi  >  PREFIX.nebo'

    def __init__(self, label='xespresso', prefix = None, images=[Atoms('')], package = 'neb', parallel = '',
                 queue = None, debug = False,
                 **kwargs):
        """

        """
        atoms = images[0]
        Espresso.__init__(self, label = label, prefix = prefix, atoms = atoms, package = package, parallel = parallel,
                 queue = queue, debug = debug, **kwargs)
        self.images = images
    def write_input(self, images, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, self.images, properties, system_changes)
        write_neb_in(self.label + '.nebi', self.images, **self.parameters)
        write_espresso_asei(self.label + '.asei', self.images, self.parameters)
    def get_neb_path_energy(self, images = None):
        if images is None:
            images = self.images
        paths, energies = self.get_property('neb', images)
        return paths, energies
    def read_results(self):
        paths, energies, error = self.read_dat()
        try:
            paths, energies, error = self.read_dat()
            self.paths = paths
            self.energies = energies
            self.error = error
            self.results['neb'] = [self.paths, self.energies]
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
        try:
            positions, gradients = self.read_path()
            self.positions = positions
            self.gradients = gradients
        except Exception as e:
            print('\nread path failed\n')
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
        images = io.read(self.directory+'/%s.xyz' % self.prefix, index = ':')
        return images
    def read_path(self):
        filename = os.path.join(self.directory, '%s.path'%self.prefix)
        nimages = len(self.images)
        natoms = len(self.images[0])
        image_gradients = []
        image_positions = []
        with open(filename, 'r') as f:
            pwo_lines = f.readlines()
            for i in range(nimages):
                index = i*(natoms + 2) + 9
                # position
                positions = [
                    [float(x) for x in line.split()[0:3]] for line
                    in pwo_lines[index:index + natoms]]
                positions = np.array(positions) * units['Bohr']
                image_positions.append(positions)
                # gradients
                gradients = [
                    [float(x) for x in line.split()[3:6]] for line
                    in pwo_lines[index:index + natoms]]
                gradients = np.array(gradients) * units['Ry'] / units['Bohr']
                image_gradients.append(-gradients)
        return image_positions, image_gradients
    def get_fit(self):
        """Returns the parameters for fitting images to band."""
        images = self.images
        R = [atoms.positions for atoms in images]
        R = self.positions
        E = self.energies
        F = self.gradients
        # print(F)
        A = images[0].cell
        pbc = images[0].pbc
        s, E, Sfit, Efit, lines = fit_raw(E, F, R, A, pbc)
        return s, E, Sfit, Efit, lines
    def plot(self):
        import matplotlib.pyplot as plt
        plt.plot(self.paths, self.energies, 'o')
        plt.plot(self.interpolation_paths, self.interpolation_energies, '--')
        plt.xlabel('Path')
        plt.ylabel('Energy (eV)')
        Eb = max(self.energies) - self.energies[0]
        Ef = max(self.energies) - self.energies[-1]
        plt.title('Eb = %1.2f(eV), Ef = %1.2f'%(Eb, Ef))
        plt.show()
    def plot_fit(self, ax=None):
        """Plots the NEB band on matplotlib axes object 'ax'. If ax=None
        returns a new figure object."""
        if not ax:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = None
        s, E, Sfit, Efit, lines = self.get_fit()
        ax.plot(s, E, 'o')
        for x, y in lines:
            ax.plot(x, y, '-g')
        ax.plot(Sfit, Efit, 'k-')
        ax.set_xlabel('Reaction path [$\AA$]')
        ax.set_ylabel('Energy profile [eV]')
        Ef = max(Efit) - E[0]
        Er = max(Efit) - E[-1]
        self.Ef = Ef
        self.Er = Er
        dE = E[-1] - E[0]
        ax.set_title('$E_\mathrm{f} \\approx$ %.3f eV; '
                     '$E_\mathrm{r} \\approx$ %.3f eV; '
                     '$\\Delta E$ = %.3f eV'
                     % (Ef, Er, dE))
        return fig


def interpolate(images, n = 10):
    '''
    Interpolate linearly the potisions of the middle temp:
    '''
    from ase.neb import NEB
    if len(images) == 2:
        new_images = [images[0]]
        for i in range(n):
            new_images.append(images[0].copy())
        new_images.append(images[1])
        neb = NEB(new_images)
        neb.interpolate()
        return new_images
    elif len(images) == 3:
        new_images1 = [images[0]]
        for i in range(n):
            new_images.append(images[0].copy())
        new_images.append(images[1])
        neb = NEB(new_images1)
        neb.interpolate()
        new_images = neb.images.copy()
        #
        new_images2 = [images[1]]
        for i in range(n):
            new_images.append(images[1].copy())
        new_images.append(images[2])
        neb = NEB(new_images2)
        neb.interpolate()
        new_images.extend(neb.images)
        return new_images

