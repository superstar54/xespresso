#!/usr/bin/env python3
from ase.build import bulk
from xespresso import XEspresso
atoms = bulk('Al', cubic = 'True')
input_data = {
'calculate': 'scf',
'verbosity': 'high', 
'ecutwfc': 40.0,
'ecutrho': 320.0,
'occupations': 'smearing',
'smearing': 'gaussian',
'degauss': 0.01,
}
queue = {'nodes': 1, 
         'ntasks-per-node': 4, 
		 'account': 'dcb', 
		 'partition': 'all', 
		 'time': '0:10:00'}
pseudopotentials = {'Al': 'Al.pbe-n-rrkjus_psl.1.0.0.UPF'}
calc = XEspresso(pseudopotentials = pseudopotentials, 
				 # queue = queue,
				 label  = 'scf/al',
				 input_data = input_data, kpts=(8, 8, 8))
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
print('Energy: ', e)
#
calc.read_results()
calc.nscf(kpts = (12, 12, 12))
calc.nscf_calculate()
#
calc.read_results()
calc.calc_pdos()
calc.plot_pdos(total = True)


'''
energy:  -468.75500017962105
'''
