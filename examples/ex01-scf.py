from ase.build import bulk
from xespresso import Espresso

atoms = bulk('Fe')
pseudopotentials = {'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'}
input_data = {'smearing': 'gaussian'}
calc = Espresso(pseudopotentials = pseudopotentials, 
                label  = 'scf/fe',  # 'scf/fe' is the directory, 'fe' is the prefix
                ecutwfc = 40,
                occupations = 'smearing',
                degauss = 0.03,
                input_data = input_data,
                kpts=(6, 6, 6))
atoms.calc = calc
e = atoms.get_potential_energy()
print('Energy: {0:1.3f}'.format(e))

'''
Energy: -3368.488
'''
