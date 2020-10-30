from math import pi, sqrt
import numpy as np
from ase.dft.kpoints import get_monkhorst_pack_size_and_offset
from ase.parallel import world
from ase.utils.cext import cextension
import matplotlib.pyplot as plt
from xespresso import Espresso
from ase.geometry import get_layers
from ase.visualize import view
import copy
import matplotlib.pyplot as plt

lspins = ['up', 'dw']
orbitals = ['s', 'p', 'd', 'f']



class DOS:
    def __init__(self, calc = None, label = None,
                 ):
        """Electronic Density Of States object.

        calc: quantum espresso calculator object
        label: the directory and prefix
        """
        if not calc and not label:
            raise ValueError('Please give one of them: calc or label')
        if label:
            calc = Espresso(label = label)
        calc.read_results()
        self.directory = calc.directory
        self.label = calc.label
        self.prefix = calc.prefix
        self.efermi = calc.get_fermi_level()
        self.atoms = calc.atoms
        self.nspins = calc.get_number_of_spins()
        if self.nspins == 0:  self.nspins = 1
    def read_dos(self):
        dos = np.loadtxt(self.directory+'/%s.dos' % self.prefix)
        if self.nspins == 2:
            self.dos = [dos[:,1], dos[:, 2]]
        else:
            self.dos = [dos[:,1]]
        self.dos_energies = dos[:,0] - self.efermi
    def read_pdos_info(self, ):
        '''
        read pdos information from output of projwfc
        state #   1: atom   1 (Al ), wfc  1 (l=0 m= 1)
        state #   2: atom   1 (Al ), wfc  2 (l=1 m= 1)
        state #   3: atom   1 (Al ), wfc  2 (l=1 m= 2)
        state #   4: atom   1 (Al ), wfc  2 (l=1 m= 3)
        state #   5: atom   2 (Al ), wfc  1 (l=0 m= 1)
        state #   6: atom   2 (Al ), wfc  2 (l=1 m= 1)
        state #   7: atom   2 (Al ), wfc  2 (l=1 m= 2)
        ...

        reg
        ['state', '1', 'atom', '1', 'Al', 'wfc', '1', 'l', '0', 'm', '1']
        Output: 
        post_info: dict,
            e.g. {'Al': {'iatom': [1, 2, 3, 4], 'istate': {1: 0, 2: 1}}}
        '''
        import re
        reg = re.compile('\w+')
        # 
        kinds = []
        pdos_info = {}
        with open(self.directory+'/%s.projwfco' % self.prefix, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'state #' in line:
                    data = reg.findall(line)
                    iatom = int(data[3])
                    kind = data[4]
                    istate = int(data[6])
                    l = int(data[8])
                    if kind not in kinds:
                        kinds.append(kind)
                        pdos_info[kind] = {}
                        pdos_info[kind]['iatom'] = []
                        pdos_info[kind]['istate'] = {}
                    if iatom not in pdos_info[kind]['iatom']:
                        pdos_info[kind]['iatom'].append(iatom)
                    if istate not in pdos_info[kind]['istate']:
                        pdos_info[kind]['istate'][istate] = l
        self.kinds = kinds
        self.pdos_info = pdos_info
    def read_pdos(self):
        '''
        projwfc
        '''
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
        self.pdos_atoms = {}
        self.pdos_kinds = {}
        for kind, info in self.pdos_info.items():
            pdos_kind = {}
            for istate, l in info['istate'].items():
                ncomponents = (2*l+1) * self.nspins + 1
                channel = '{0}{1}'.format(istate, orbitals[l])
                pdos_kind[channel] = np.zeros((ncomponents, npoints), np.float)
            for iatom in info['iatom']:
                pdos_atom = {}
                for istate, l in info['istate'].items():
                    filename = self.directory + '/{0}.pdos_atm#{1}({2})_wfc#{3}({4})'.format(self.prefix, iatom, kind, istate, orbitals[l])
                    channel = '{0}{1}'.format(istate, orbitals[l])
                    pdosinp = np.genfromtxt(filename)
                    ncomponents = (2*l+1) * self.nspins + 1
                    pdos_atom[channel] = np.zeros((ncomponents, npoints), np.float)
                    for j in range(ncomponents):
                        pdos_atom[channel][j] += pdosinp[:, j + 1]
                         # sum over kind
                self.pdos_atoms[iatom] = pdos_atom
            self.pdos_kinds[kind] = pdos_kind
        self.pdos_kinds = self.merge_kinds(index = range(len(self.atoms)))
        # sum over orbital
        return self.pdos_energies, self.pdos_tot, self.pdos_atoms, self.pdos_kinds
    def merge_kinds(self, index = []):
        '''
        merge atoms in the index
        '''
        pdos_kinds = {}
        npoints = len(self.pdos_energies)
        for kind, info in self.pdos_info.items():
            # print(kind, info)
            pdos_kind = self.pdos_kinds[kind]
            for channel, pdos in pdos_kind.items():
                # print(channel)
                pdos[:, :] = 0
            for iatom in info['iatom']:
                if iatom - 1 not in index: continue
                # print(iatom)
                for channel, pdos_atom in self.pdos_atoms[iatom].items():
                    pdos_kind[channel] += pdos_atom
            pdos_kinds[kind] = pdos_kind
        return pdos_kinds
    def read_proj(self, ):
        '''
        read projwfc output file *.up and *.dw
        '''
        #
        self.read_pdos_info()
        proj_files = [self.directory + '/%s.projwfc_up'%self.prefix, self.directory + '/%s.projwfc_down'%self.prefix]
        projs = []
        projs_kinds = []
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
                proj_kinds = {}
                for istate in range(nstates):
                    iline = j + istate*(nkpts*nbands + 1) - 1
                    kinds = lines[iline].split()[2]
                    orbital = lines[iline].split()[3]
                    channel = kinds + orbital
                    if channel not in channels: 
                        channels.append(channel)
                        proj_kinds[channel] = np.zeros((nkpts, nbands), np.float)
                    # l = lines[iline].split()[2]
                    # print(kinds, istate)
                    for ikpt in range(nkpts):
                        for iband in range(nbands):
                            iline = j + istate*(nkpts*nbands + 1) + ikpt*nbands + iband
                            # print('line: ', iline)
                            proj[istate, ikpt, iband] = float(lines[iline].split()[2])
                            proj_kinds[channel][ikpt, iband] += float(lines[iline].split()[2])
                            # proj_atoms[channel, ikpt, iband] += 
            projs.append(proj)
            projs_kinds.append(proj_kinds)
        self.projs = projs
        self.projs_kinds = projs_kinds
    def plot_data(self, energies, dos, label, xindex, ax = None, fill = True, color = None):
        if ax is None:
            ax = plt.gca()
        for i in range(self.nspins):
            if self.nspins == 2:
                newlabel = '%s-%s' % (label, lspins[i])
            else:
                newlabel = label
            if color:
                ax.plot(energies[xindex], (-1)**i*dos[i][xindex], linewidth=0.7, label = newlabel, color = color)
            else:
                ax.plot(energies[xindex], (-1)**i*dos[i][xindex], linewidth=0.7, label = newlabel)
            if fill:
                if color:
                    ax.fill_between(energies[xindex], (-1)**i*dos[i][xindex], 0, alpha = 0.2, color = color)
                else:
                    ax.fill_between(energies[xindex], (-1)**i*dos[i][xindex], 0, alpha = 0.2)
        return ax
    def plot_dos(self, energies = None, dos = None, 
        Emin = -5, Emax = 5, color = None,
        ax = None, fill = True, output = None):
        '''
        '''
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(figsize = (6, 3))
            # ax = plt.gca()
        if dos is None: dos = self.dos
        if energies is None: energies = self.dos_energies
        xindex = (energies>Emin) & (energies<Emax)
        ax = self.plot_data(energies, dos, label = 'dos', ax = ax, xindex = xindex, fill = fill, color = color)
        # ax.legend()
        # plt.grid(linestyle = '--')
        # ax.set_xlabel('Energy (eV)')
        # ax.set_ylabel('DOS (a.u.)')
        # ax.set_title('%s' % self.prefix)
        if output is None:
            output = '{0}-dos.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return ax
    def plot_pdos_tot(self, energies = None, dos = None, 
        Emin = -5, Emax = 5,
        ax = None, fill = True, output = None):
        # print(self.pdos_tot)
        if dos is None: dos = self.pdos_tot
        if energies is None: energies = self.pdos_energies
        xindex = (energies>Emin) & (energies<Emax)
        fig, ax = plt.subplots(figsize = (6, 3))
        ax = self.plot_data(energies, dos, label = 'pdos', ax = ax, xindex = xindex, fill = fill)
        ax.legend()
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('PDOS (a.u.)')
        ax.set_title('%s' % self.prefix)
        if output is None:
            output = '{0}-pdos-tot.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return ax
    def plot_pdos(self, energies = None, pdos_kinds = None, 
                  Emin = -5, Emax = 5,
                  ax = None, total = False, select = None, fill = True, 
                  output = None, legend = False, xylabel = True):
        '''
        '''
        if energies is None: energies = self.pdos_energies
        if pdos_kinds is None: pdos_kinds = self.pdos_kinds
        xindex = (energies>Emin) & (energies<Emax)
        # print(xindex)
        if ax is None:
            fig, ax = plt.subplots(figsize = (6, 3))
        # if total:
            # ax = self.plot_data(self.pdos_energies, self.pdos_tot, label = 'pdos', ax = ax, fill = fill)
        for kind, channels in pdos_kinds.items():
            if select and kind not in select: continue
            for channel, pdos in channels.items():
                if select and channel[-1] not in select[kind]: continue
                for i in range(self.nspins):
                    label = '{0}-{1}'.format(kind, channel)
                    ax = self.plot_data(energies, pdos, label = label, ax = ax, xindex = xindex, fill = fill)
        if legend: ax.legend()
        if xylabel:
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('PDOS (a.u.)')
        if output is None:
            output = '{0}-pdos.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        return ax
    def plot_pdos_layer(self, atoms = None, 
                        layers = None,
                        miller = (0, 0, 1), 
                        dz = 0.3, 
                        Emin = -5, 
                        Emax = 3, 
                        total = False, 
                        select = None, 
                        fill = True, 
                        output = None):
        '''
        
        '''
        atoms = self.atoms
        if not layers:
            layers = get_layers(atoms, miller, dz)[0]
        nlayers = max(layers) + 1
        fig, axs = plt.subplots(nlayers, 1, figsize = (10, nlayers*2), sharex = True)
        indexs = range(len(atoms))
        images = []
        xindex = (self.pdos_energies>Emin) & (self.pdos_energies<Emax)
        for ilayer in range(nlayers):
            iax = nlayers - ilayer - 1
            # axs[iax].set_yticks([])
            index = [j for j in indexs if layers[j] ==  ilayer]
            images.append(atoms[index])
            pdos_kinds = self.merge_kinds(index)
            self.plot_pdos(energies = self.pdos_energies, pdos_kinds = pdos_kinds, select = select, 
                           Emin = Emin, Emax = Emax, ax = axs[iax], 
                           legend = False, xylabel = False)
            axs[iax].axvline(0, color = 'b')
            # print(iax, ilayer)
            # print(index)
            # axs[iax].set_xlim([Emin, Emax])
        plt.xlabel('E - E$_{Fermi}$ (eV)', size='16')
        plt.subplots_adjust(hspace=0)
        fig.text(0.05, 0.5, 'PDOS (a.u.)', size='16', va='center', rotation='vertical')
        if output is None:
            output = '{0}-pdos.png'.format(self.prefix) 
            plt.savefig('%s' %output)
        # view(images)
        return axs, images
    def get_pdos(self, species, orbital):
        dos = self.pdos_kinds[species][orbital]
        return self.pdos_energies, dos
    def get_d_band_center(sef, species = None, orbital = 'd'):
        pdos_kinds
        energies, dos = self.get_pdos(species, orbital)
        Nstates = np.trapz(dos, energies)
        occupied = energies <= 0.0
        N_occupied_states = np.trapz(dos[occupied], energies[occupied])
        # first moment
        ed = np.trapz(energies * dos, energies) / Nstates
        # second moment
        wd2 = np.trapz(energies**2 * dos, energies) / Nstates
        return d_center, d_width
    
def compare_pdos(dos1, index1, dos2, index2):
    '''
    '''
    if len(index1) != len(index2):
        assert ('number of atoms wrong!')
    natoms = len(index1)
    for i in range(natoms):
        if dos1.atoms[i].symbol != dos2.atoms[i].symbol:
            assert('atoms symbol not match!')
    ddos = copy.deepcopy(dos1)
    # print(dos1.atoms.get_chemical_symbols())
    # print(dos2.atoms.get_chemical_symbols())
    atoms = ddos.atoms
    for i in range(natoms):
        iatom1 = index1[i]
        iatom2 = index2[i]
        for channel in ddos.pdos_atoms[iatom1 + 1]:
            # print(i, atoms[i].symbol, channel)
            npoints = len(ddos.pdos_atoms[iatom1 + 1][channel][0])
            ddos.pdos_atoms[iatom1 + 1][channel] = dos1.pdos_atoms[iatom1 + 1][channel][:, 0:npoints] - dos2.pdos_atoms[iatom2 + 1][channel][:, 0:npoints]
    ddos.pdos_kinds = ddos.merge_kinds(index = range(len(ddos.atoms)))
    return ddos