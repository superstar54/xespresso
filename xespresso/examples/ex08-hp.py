from ase.build import bulk
from ase.visualize import view
from xespresso import Espresso
import numpy as np

atoms = bulk('Fe', cubic = True)
# set new species for afm state
atoms.new_array('species', np.array(atoms.get_chemical_symbols(), dtype = 'U20'))
atoms.arrays['species'][1] = 'Fe1'
input_ntyp = {'starting_magnetization': {'Fe': 1.0, 'Fe1': -1.0, }}

input_data = {
'verbosity': 'high', 
'ecutwfc': 40.0,
'ecutrho': 320.0,
'occupations': 'smearing',
'smearing': 'gaussian',
'degauss': 0.03,
'nspin': 2,
'input_ntyp': input_ntyp,
#
'mixing_beta': 0.3,
'conv_thr': 1.0e-8,
}

pseudopotentials = {
'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
'Fe1': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF',
}
calc = Espresso(pseudopotentials = pseudopotentials, 
                label  = 'scf/fe-afm',
                input_data = input_data, kpts=(6, 6, 6))
atoms.calc = calc
e = atoms.get_potential_energy()
calc.read_results()
e = calc.results['energy']
print('Energy  {0:1.3f}'.format(e))


parameters_hp = {
'nq1': 1,
'nq2': 1,
'nq3': 1,
'equiv_type': 0,
'conv_thr_chi': 1e-5,
'perturb_only_atom': {'1':True}   #  key(i)
}
calc.post(package='hp', **parameters_hp)


'''
Energy  -6737.122
'''
