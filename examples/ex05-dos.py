from ase.build import bulk
from ase.visualize import view
from xespresso import Espresso
from xespresso.dos import DOS
from xespresso.tools import get_nbnd
import matplotlib.pyplot as plt

atoms = bulk('Fe')
pseudopotentials = {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
calc = Espresso(pseudopotentials = pseudopotentials, 
                label  = 'scf/fe',  # 'scf/fe' is the directory, 'fe' is the prefix
                ecutwfc = 40,
                occupations = 'smearing',
                degauss = 0.03,
                kpts=(6, 6, 6))
atoms.calc = calc
# e = atoms.get_potential_energy()
calc.read_results()
fermi = calc.get_fermi_level()
#===============================================================
# nscf calculation
calc.nscf(kpts=(8, 8, 8))
calc.nscf_calculate()
# post calculation
calc.post(package='dos', Emin = fermi - 20, Emax = fermi + 10, DeltaE = 0.1)
calc.post(package='projwfc', Emin = fermi - 20, Emax = fermi + 10, DeltaE = 0.1)
# # DOS analysis
dos = DOS(calc)
dos.read_dos()
dos.plot_dos()
plt.show()
#

'''
Energy: -3368.435
'''
