from ase.build import fcc111
atoms = fcc111('Al', size=(1, 1, 2), vacuum = 4.0)
# view(atoms)
from gpaw import GPAW
calc = GPAW(mode='pw',
            kpts=(8, 8, 1),
            symmetry='off',
            txt='al111.txt')
atoms.calc = calc
energy = atoms.get_potential_energy()
calc.write('al111.gpw', 'all')

# Creates: 2d.png, 2d_I.png, line.png, dIdV.png
from ase.dft.stm import STM
from gpaw import GPAW
calc = GPAW('al111.gpw')
atoms = calc.get_atoms()
stm = STM(atoms)
z = 8.0
bias = 1.0
c = stm.get_averaged_current(bias, z)
x, y, h = stm.scan(bias, c, repeat=(3, 5))

import matplotlib.pyplot as plt
plt.gca(aspect='equal')
plt.contourf(x, y, h, 40)
plt.colorbar()
plt.savefig('images/2d_H-gpaw.png')

plt.figure()
plt.gca(aspect='equal')
z = 7.8
x, y, I = stm.scan2(bias, z, repeat=(3, 5))
plt.contourf(x, y, I, 40)
plt.colorbar()
plt.savefig('images/2d_I-gpaw.png')