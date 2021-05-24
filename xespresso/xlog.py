from ase import __version__ as ase_version
from ase.utils import convert_string_to_fd
from ase.utils import search_current_git_hash
from ase.visualize import view
import numpy as np
import os
import time
import sys
import xespresso
import ase

class XLogger():
    """Class for handling all text output."""
    def __init__(self):
        self.verbose = False
        self._fd = None
        self.oldfd = None

    @property
    def fd(self):
        return self._fd

    @fd.setter
    def fd(self, fd):
        """Set the stream for text output.
        """
        if fd == self.oldfd:
            return
        self.oldfd = fd
        self._fd = convert_string_to_fd(fd)
        self.logo()
        self.header()

    def __call__(self, *args, **kwargs):
        flush = kwargs.pop('flush', False)
        print(*args, file=self._fd, **kwargs)
        if flush:
            self._fd.flush()

    def flush(self):
        self._fd.flush()


    def logo(self):
        self(' xespresso  ')

    def header(self):
        self()
        nodename, machine = os.uname()[1::3]
        self('User:  ', os.getenv('USER', '???') + '@' + nodename)
        self('Date:  ', time.asctime())
        self('Arch:  ', machine)
        self('Pid:   ', os.getpid())
        self('Python: {0}.{1}.{2}'.format(*sys.version_info[:3]))
        # Espresso
        line = os.path.dirname(xespresso.__file__)
        githash = search_current_git_hash(xespresso)
        if githash is not None:
            line += ' ({:.10})'.format(githash)
        self('xespresso:  ', line)

        # ASE
        line = '%s (version %s' % (os.path.dirname(ase.__file__), ase_version)
        githash = search_current_git_hash(ase)
        if githash is not None:
            line += '-{:.10}'.format(githash)
        line += ')'
        self('ase:   ', line)
        self('\n'*2)

    def print_atoms(self, atoms):
        self()
        self('{0:15s}: {1}'.format('    Formula', atoms.get_chemical_formula()))
        self('{0:15s}: {1}'.format('    Cell', np.round(np.asarray(atoms.cell.cellpar()), 3)))
        self('{0:15s}: {1}'.format('    PBC', atoms.pbc))
        self('{0:15s}:'.format('    Info'))
        self.print_dict(atoms.info, sep = '        ')

    def print_calculator(self, calc):
        self()
        self('Calculator:')
        self('Package: ')
        self('  k-point: ', )
        self('  Pseudopotentials: \n')
        
    def print_dict(self, dct, sep='  '):
        options = np.get_printoptions()
        try:
            np.set_printoptions(threshold=4, linewidth=50)
            for key, value in sorted(dct.items()):
                if isinstance(value, dict):
                    sep = ',\n     ' + ' ' * len(key)
                    keys = sorted(value, key=lambda k: (str(type(k)), k))
                    s = sep.join('{0}: {1}'.format(k, value[k]) for k in keys)
                    self('  {0}: {{{1}}}'.format(key, s))
                elif hasattr(value, '__len__'):
                    value = np.asarray(value)
                    sep = ',\n    ' + ' ' * len(key)
                    s = sep.join(str(value).splitlines())
                    self('  {0}: {1}'.format(key, s))
                else:
                    self('  {0}: {1}'.format(key, value))
        finally:
            np.set_printoptions(**options)