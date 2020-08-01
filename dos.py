from math import pi, sqrt
import numpy as np
from ase.dft.kpoints import get_monkhorst_pack_size_and_offset
from ase.parallel import world
from ase.utils.cext import cextension
import matplotlib.pyplot as plt

lspins = ['up', 'dw']
orbitals = ['s', 'p', 'd', 'f']

class DOS:
    def __init__(self, calc, dos = True, pdos = False, width=0.1, 
                 Emin = -10, Emax = 5, npts=401):
        """Electronic Density Of States object.

        calc: calculator object
            Any ASE compliant calculator object.
        width: float
            Width of guassian smearing.  Use width=0.0 for linear tetrahedron
            interpolation.
        window: tuple of two float
            Use ``window=(emin, emax)``.  If not specified, a window
            big enough to hold all the eigenvalues will be used.
        npts: int
            Number of points.

        """
        self.directory = calc.directory
        self.label = calc.label
        self.prefix = calc.prefix
        self.efermi = calc.get_fermi_level()
        self.atoms = calc.atoms
        self.npts = npts
        self.width = width
        self.w_k = calc.get_k_point_weights()
        self.nspins = calc.get_number_of_spins()
        if self.nspins == 0:  self.nspins = 1
        if dos:
            self.read_dos()
        if pdos:
            self.read_pdos()
    def read_dos(self):
        dos = np.loadtxt(self.directory+'/%s.dos' % self.prefix)
        if self.nspins == 2:
            self.dos = [dos[:,1], dos[:, 2]]
        else:
            self.dos = [dos[:,1]]
        self.dos_energies = dos[:,0] - self.efermi
    def read_pdos_info(self, ):
        import re
        reg = re.compile('\w+')
        species = []
        pdos_info = {}
        with open(self.directory+'/%s-projwfc.out' % self.prefix, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'state #' in line:
                    data = reg.findall(line)
                    istate = int(data[3])
                    iatom = int(data[3])
                    label = data[4]
                    ichannel = int(data[6])
                    l = int(data[8])
                    m = int(data[10])
                    # print(iatom, label, ichannel, l, m)
                    if label not in species:
                        species.append(label)
                        pdos_info[label] = {}
                        pdos_info[label]['iatom'] = []
                        pdos_info[label]['ichannel'] = {}
                        # pdos_info[label]['orbital'] = []
                    if iatom not in pdos_info[label]['iatom']:
                        pdos_info[label]['iatom'].append(iatom)
                    if ichannel not in pdos_info[label]['ichannel']:
                        pdos_info[label]['ichannel'][ichannel] = l
                        # pdos_info[label]['orbital'].append(orbitals[l])
        self.species = species
        self.pdos_info = pdos_info

    def read_pdos(self):
        self.read_pdos_info()
        # read in total density of states
        dos = np.loadtxt(self.directory+'/%s.pdos_tot' % self.prefix)
        if self.nspins == 2:
            self.pdos_tot = [dos[:,1], dos[:, 2]]
        else:
            self.pdos_tot = [dos[:,1]]
        self.pdos_energies = dos[:,0] - self.efermi
        npoints = len(self.pdos_energies)
        # read in projections onto atomic orbitals
        self.natoms = len(self.atoms)
        self.nspecies = len(self.species)
        self.pdos = []
        self.pdos_atom = {}
        self.pdos_species = {}
        for species, info in self.pdos_info.items():
            pdos_species = {}
            for ichannel, l in info['ichannel'].items():
                ncomponents = (2*l+1) * self.nspins + 1
                # print(ncomponents)
                channel = '{0}{1}'.format(ichannel, orbitals[l])
                pdos_species[channel] = np.zeros((ncomponents, npoints), np.float)
            for iatom in info['iatom']:
                pdos_atom = {}
                for ichannel, l in info['ichannel'].items():
                    # print(species, ichannel, l)
                    filename = self.directory + '/{0}.pdos_atm#{1}({2})_wfc#{3}({4})'.format(self.prefix, iatom, species, ichannel, orbitals[l])
                    channel = '{0}{1}'.format(ichannel, orbitals[l])
                    pdosinp = np.genfromtxt(filename)
                    # print(inpfile, channel, jpos)
                    ncomponents = (2*l+1) * self.nspins + 1
                    pdos_atom[channel] = np.zeros((ncomponents, npoints), np.float)
                    # pdos_species[channel][0] += pdosinp[:, 0]
                    for j in range(ncomponents):
                        # print(j)
                        pdos_atom[channel][j] += pdosinp[:, j + 1]
                         # sum over species
                        pdos_species[channel][j] += pdosinp[:,j + 1]
                self.pdos_atom[iatom] = pdos_atom
            self.pdos_species[species] = pdos_species
        # sum over orbital
        # print(self.pdos_species)
        return self.pdos_energies, self.pdos_tot, self.pdos_atom, self.pdos_species
    def read_proj(self, ):
        '''
        read projwfc output file *.up and *.dw
        '''
        #
        self.read_pdos_info()
        proj_files = [self.directory + '/%s.projwfc_up'%self.prefix, self.directory + '/%s.projwfc_down'%self.prefix]
        projs = []
        projs_species = []
        channels = []
        for i in range(self.nspins):
            with open(proj_files[i]) as f:
                lines = f.readlines()
                for j in range(len(lines)):
                    if '       1       1' in lines[j]: break
                nstates = int(lines[j - 3].split()[0])
                nkpts = int(lines[j - 3].split()[1])
                nbands = int(lines[j - 3].split()[2])
                print(nstates, nkpts, nbands)
                proj = np.zeros((nstates, nkpts, nbands), np.float)
                proj_species = {}
                for istate in range(nstates):
                    iline = j + istate*(nkpts*nbands + 1) - 1
                    species = lines[iline].split()[2]
                    orbital = lines[iline].split()[3]
                    channel = species + orbital
                    if channel not in channels: 
                        channels.append(channel)
                        proj_species[channel] = np.zeros((nkpts, nbands), np.float)
                    # l = lines[iline].split()[2]
                    # print(species, ichannel)
                    for ikpt in range(nkpts):
                        for iband in range(nbands):
                            iline = j + istate*(nkpts*nbands + 1) + ikpt*nbands + iband
                            # print('line: ', iline)
                            proj[istate, ikpt, iband] = float(lines[iline].split()[2])
                            proj_species[channel][ikpt, iband] += float(lines[iline].split()[2])
                            # proj_atoms[channel, ikpt, iband] += 
            projs.append(proj)
            projs_species.append(proj_species)
        self.projs = projs
        self.projs_species = projs_species
    def plot_data(self, energies, dos, label, ax = None, fill = True):
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        for i in range(self.nspins):
            if self.nspins == 2:
                newlabel = '%s-%s' % (label, lspins[i])
            else:
                newlabel = label
            ax.plot(energies, (-1)**i*dos[i], linewidth=0.7, label = newlabel)
            if fill:
                ax.fill_between(energies, (-1)**i*dos[i], 0, alpha = 0.2)
        return ax
    def plot_dos(self, ax = None, fill = True, output = None):
        '''
        '''
        fig, ax = plt.subplots(figsize = (6, 3))
        ax = self.plot_data(self.dos_energies, self.dos, label = 'dos', ax = ax, fill = fill)
        ax.legend()
        # plt.grid(linestyle = '--')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('DOS (a.u.)')
        ax.set_title('%s' % self.prefix)
        if not output:
            output = '{0}-dos.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return ax
    def plot_pdos_tot(self, fill = True, output = None):
        # print(self.pdos_tot)
        fig, ax = plt.subplots(figsize = (6, 3))
        ax = self.plot_data(self.pdos_energies, self.pdos_tot, label = 'pdos', ax = ax, fill = fill)
        ax.legend()
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('PDOS (a.u.)')
        ax.set_title('%s' % self.prefix)
        if not output:
            output = '{0}-pdos-tot.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return ax
    def plot_pdos(self, total = False, select = None, fill = True, output = None):
        '''
        '''
        fig, ax = plt.subplots(figsize = (6, 3))
        if total:
            ax = self.plot_data(self.pdos_energies, self.pdos_tot, label = 'pdos', ax = ax, fill = fill)
        for species, channels in self.pdos_species.items():
            if select and species not in select: continue
            for channel, pdos in channels.items():
                if select and channel[-1] not in select[species]: continue
                for i in range(self.nspins):
                    label = '{0}-{1}'.format(species, channel)
                    ax = self.plot_data(self.pdos_energies, pdos, label = label, ax = ax, fill = fill)
        ax.legend()
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('PDOS (a.u.)')
        ax.set_title('%s' % self.prefix)
        if not output:
            output = '{0}-pdos.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return ax
    def plot_pdos_layer(self, atoms = None, 
                        miller = (0, 0, 1), 
                        tolerance = 0.5, 
                        total = False, 
                        select = None, 
                        fill = True, 
                        output = None):
        
        from ase.geometry import get_layers
        atoms = self.atoms
        layers = get_layers(atoms, miller, tolerance)[0]
        nlayers = max(layers)
        fig, axs = plt.subplots(nlayers, 1, figsize = (6, 3), sharex = True)
        for i in range(nlayers):
            self.plot_pdos()
            if total:
                ax = self.plot_data(self.pdos_energies[layers==i], pdos_tot[layers==i], label = 'pdos', ax = ax[i], fill = fill)
            for species, channels in self.pdos_species.items():
                if select and species not in select: continue
                for channel, pdos in channels.items():
                    if select and channel[-1] not in select[species]: continue
                    for i in range(self.nspins):
                        label = '{0}-{1}'.format(species, channel)
                        # ax = self.plot_data(self.pdos_energies[layers==i], pdos[layers==i], label = label, ax = ax[i], fill = fill)
                        pass
        # ax.legend()
        # plt.set_xlabel('Energy (eV)')
        # plt.set_ylabel('PDOS (a.u.)')
        # ax.set_title('%s' % self.prefix)
        if not output:
            output = '{0}-pdos.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return axs