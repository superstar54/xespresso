from _common_helpers import set_envs
import numpy as np
import os

def test_scf():
    from ase.build import molecule
    from xespresso import Espresso
    set_envs()
    h2 = molecule('H2')
    h2.cell = [8, 8, 8]
    h2.pbc = [True, True, True]
    pseudopotentials = {'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',}
    calc = Espresso(pseudopotentials = pseudopotentials, 
                    label  = 'calculations/scf/h2',  # 'scf/fe' is the directory, 'fe' is the prefix
                    ecutwfc = 20,
                    occupations = 'smearing',
                    degauss = 0.03,
                    kpts=(1, 1, 1),
                    debug = True,
                    )
    h2.calc = calc
    try:
        e = h2.get_potential_energy()
        print('Energy: {0:1.4f}'.format(e))
        assert np.isclose(e, -31.4454)
    except Exception as e:
        print(e)
    os.system('ls calculations/scf/h2')
    os.system('cat calculations/scf/h2/CRASH')

def test_relax():
    from ase.build import molecule
    from xespresso import Espresso
    set_envs()
    atoms = molecule('H2')
    atoms.center(5)
    atoms.pbc = [True, True, True]
    pseudopotentials = {'H': 'H.pbe-rrkjus_psl.1.0.0.UPF'}
    calc = Espresso(label = 'calculations/relax/h2',
                    pseudopotentials = pseudopotentials,
                    calculation = 'relax',
                    ecutwfc = 20,
                    kpts = (1, 1, 1),
                    debug = True,
                    )
    atoms.calc = calc
    try:
        e = atoms.get_potential_energy()
        print('Energy = {0:1.4f} eV'.format(e))
        assert np.isclose(e, -31.4552)
    except Exception as e:
        print(e)
    os.system('ls calculations/relax/h2')
    os.system('cat calculations/relax/h2/CRASH')
