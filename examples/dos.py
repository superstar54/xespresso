from ase.build import bulk
from ase.visualize import view
from xespresso import Espresso
from xespresso.dos import DOS
import matplotlib.pyplot as plt

# build Fe
atoms = bulk('Fe')
# input parameters for pw
input_data = {
'verbosity': 'high', 
'ecutwfc': 40.0,
'ecutrho': 302.0,
'occupations': 'smearing',
'smearing': 'gaussian',
'degauss': 0.03,
#
'mixing_beta': 0.3,
'conv_thr': 1.0e-8,
}
pseudopotentials = {
'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
}
calc = Espresso(pseudopotentials = pseudopotentials, 
                label  = 'scf/fe/fe',      # 'scf/fe' is the directory, 'fe' is the prefix
                input_data = input_data, 
                kpts=(6, 6, 6))
atoms.calc = calc
# e = atoms.get_potential_energy()
# print('Energy: {0:1.3f}'.format(e))
# read previous scf results

#===============================================================
# start nscf calculation, and dos, projwfc
calc.read_results()
fermi = calc.get_fermi_level()
# nscf calculation
calc.nscf(kpts=(10, 10, 10))
calc.nscf_calculate()
calc.read_results()
# post calculation
calc.post(package='dos', Emin = fermi - 20, Emax = fermi + 10, DeltaE = 0.1)
calc.post(package='projwfc', Emin = fermi - 20, Emax = fermi + 10, DeltaE = 0.1)
# DOS analysis
dos = DOS(calc, dos = True, pdos = True)
dos.read_pdos()
dos.plot_dos()
plt.show()
#

'''
Energy: -3368.435
'''
