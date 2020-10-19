'''
H2 + H --> H + H2
'''
from ase.atoms import Atoms
from ase.constraints import FixAtoms, FixCartesian
from ase.visualize import view
from xespresso import Espresso
from xespresso.neb import NEBEspresso

#=============================================================
# first optmize the inital and final structure
atoms = Atoms('H3', positions = [[0, 0, 0], [0.8, 0, 0.0], [3.0, 0, 0]], 
                    cell = [12.0, 5, 5],
                    pbc = [True, True, True])
constraints = FixAtoms(indices = [0, 2])
atoms.set_constraint(constraints)
# view(atoms)
input_data = {
'calculation': 'relax',
'verbosity': 'high', 
'ecutwfc': 20.0,
'ecutrho': 100.0,
'occupations': 'smearing',
'smearing': 'gaussian',
'degauss': 0.03,
'nspin': 2,
'starting_magnetization': 0.5E0,
#
'mixing_beta': 0.3,
'conv_thr': 1.0e-8,
}
queue = {'nodes': 1, 'ntasks-per-node': 4, 'account': 'dcb', 
         'partition': 'all', 'time': '0:59:00'}
pseudopotentials = {
'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
}
calc = Espresso(pseudopotentials = pseudopotentials, 
                #  queue = queue,
                 label  = 'relax/initial/initial',
                 input_data = input_data, kpts=(1, 1, 1))
atoms.calc = calc
# atoms.get_potential_energy()
calc.read_results()
print('Initial energy: %s'%calc.results['energy'])


#=============================================================
# start NEB 

initial = calc.results['atoms']
final = initial.copy()
final.positions[1] = [3.0 - final[1].x, 0, 0]
images = [initial, final]
# view(images)
#
path_data = {
# path
'restart_mode': 'from_scratch',
'string_method' : 'neb',
'nstep_path': 20,
'ds': 2.E0,
'opt_scheme': "broyden",
'num_of_images' : 7,
'k_max' : 0.3E0,
'k_min' : 0.2E0,
'CI_scheme' : "auto",
'path_thr': 0.1E0,
}
# print(queue)
calc = NEBEspresso(pseudopotentials = pseudopotentials, 
                 package = 'neb',
                #  queue = queue,
                 images = images,
                 climbing_images = [5],
                 label  = 'neb/h3/h3',
                 input_data = input_data, 
                 path_data = path_data,
                 kpts=(1, 1, 1))
#
# calc.calculate()
calc.read_results()
calc.plot()

