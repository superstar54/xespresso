from ase.build import bulk
from ase.visualize import view
from xespresso import Espresso

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
e = atoms.get_potential_energy()
print('Energy: {0:1.3f}'.format(e))

'''
Energy: -3368.435
'''
