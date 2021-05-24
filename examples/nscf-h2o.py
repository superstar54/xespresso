from ase.build import molecule
from xespresso import Espresso
from pseudo import mypseudo

mol = molecule('H2O')
mol.cell = [10, 10, 10]       # Unit cell
mol.pbc = [True, True, True]  # We always use the periodic boundary condition for our DFT calculation


queue = {'nodes': 1, 'ntasks-per-node': 10, 'partition': 'epyc2', 'mem-per-cpu': '5G', 'time': '23:59:00', 'config': '.xespresso-intel-2020b'}
pseudopotentials = {'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
                    'H': 'H.pbe-rrkjus_psl.1.0.0.UPF'}


calc = Espresso(label = 'relax/h2o',
                pseudopotentials = pseudopotentials, 
                calculation = 'relax',   # allow atoms to move
                ecutwfc = 60,
                ecutrho = 480.0,
                occupations = 'smearing',
                degauss = 0.005,
                kpts = (1, 1, 1),
                queue = queue,
		debug = True)
calc.run(atoms = mol)
calc.read_results()
E = calc.results['energy']
print('Energy {0:3.3f}'.format(E) + ' eV')

fe = calc.get_fermi_level()
calc.nscf(queue = queue, occupations = 'tetrahedra', kpts = [2, 2, 2])
calc.nscf_calculate()
calc.post(queue = queue, package = 'dos', Emin = fe - 30, Emax = fe + 30, DeltaE = 0.01)
calc.post(queue = queue, package = 'projwfc', Emin = fe - 30, Emax = fe + 30, DeltaE = 0.01)

