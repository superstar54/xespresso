from ase.build import molecule
from ase.visualize import view
from xespresso import Espresso
from xespresso.dos import DOS
import matplotlib.pyplot as plt

atoms = molecule('CO')
atoms.center(5)
atoms.pbc = [True, True, True]
pseudopotentials = {'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
                    'C': 'C.pbe-n-rrkjus_psl.1.0.0.UPF'}
calc = Espresso(label = 'scf/co/co',
                pseudopotentials = pseudopotentials,
                ecutwfc = 30,
                kpts = (1, 1, 1),
                )
atoms.calc = calc
atoms.get_potential_energy()
calc.read_results()
print('Energy = {0:1.3f} eV'.format(calc.results['energy']))
#===============================================================
# start nscf calculation, and dos, projwfc
calc.read_results()
fermi = calc.get_fermi_level()
calc.nscf(kpts = (1, 1, 1))
calc.nscf_calculate()
calc.post(package = 'dos', Emin = fermi - 20, Emax = fermi + 10, DeltaE = 0.1)
# DOS analysis
dos = DOS(calc)
dos.read_dos()
dos.plot_dos()
plt.savefig('images/co-dos.png')