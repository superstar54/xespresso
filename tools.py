import os
from ase.geometry import get_layers
from ase.constraints import FixAtoms
from xespresso import Espresso
from ase.dft.bandgap import bandgap
import pickle
import multiprocessing
import numpy as np
#====================================================
# tools
def merge_slab(slab1, slab2, index = 2):
    '''
    '''
    slab2.cell[index] = slab1.cell[index]
    slab2.set_cell(slab1.cell, scale_atoms = True)
    slab1 = slab1 + slab2
    return slab1
def build_oer(atoms):
    '''
    '''
    from ase.atoms import Atoms
    #----------------------------------
    ooh = Atoms('O2H', positions = [[0, 0, 0], [1.4, 0, 0], [1.4, 0, 1.0]])
    ooh.rotate('z', np.pi/4)
    mols = {
    'o' : Atoms('O'),
    'oh' : Atoms('OH', positions = [[0, 0, 0], [0, 0, 1.0]]),
    'ooh' : ooh,
    }
    #
    jobs = {}
    maxz = max(atoms.positions[:, 2])
    indm = [atom.index for atom in atoms if atom.z > maxz - 1.0 and atom.symbol != 'O'][0]
    for job, mol in mols.items():
        # print(job, mol)
        ads = mol.copy()
        natoms = atoms.copy()
        ads.translate(atoms[indm].position - ads[0].position + [0, 0, 1.9])
        natoms = natoms + ads
        jobs[job] = natoms
    return jobs

def fix_layers(atoms, miller, tol = 1.0, n = [0, 4]):
    '''
    '''
    layers = get_layers(atoms, miller, tol)[0]
    index = [j for j in range(len(atoms)) if layers[j] in range(n[0], n[1])]
    constraint = FixAtoms(indices=index)
    atoms.set_constraint(constraint)
    return atoms

def mypool(jobs, func):
    '''
    '''
    pool = multiprocessing.Pool(len(jobs))
    results = []
    images = []
    for job, atoms in jobs.items():
        print(job, atoms)
        results.append(pool.apply_async(func, (job, atoms)))
    for r in results:
        r.get()
    pool.close()
    pool.join()

def dwubelix(updates = []):
    file = 'pw err out dos pdos projwfc'
    print('Downloading.....')
    cwd = os.getcwd()
    for update in updates:
        os.chdir(update)
        os.system('dwubelix.py %s'%file)
        os.chdir(cwd)
    print('Finished')
def ana(dire, calc):
    atoms = calc.results['atoms']
    results = [dire, atoms, atoms.cell, atoms.positions]
    for prop in ['energy', 'forces', 'stress', 'magmoms']:
        if prop in calc.results:
            prop = calc.results[prop]
        else:
            prop = None
        results.append(prop)
    return results
# tools
def summary(updates = [], prefix = 'datas'):
    import pandas as pd
    columns = ['label', 'atoms', 'cell', 'positions', 'energy', 'forces', 'stress']
    file = '%s.pickle' % prefix
    db = '%s.db' % prefix
    if os.path.exists(file):
        with open(file, 'rb') as f:
            datas, df = pickle.load(f)
    else:
        datas = {}
        df = pd.DataFrame(columns=columns)
    calc = Espresso()
    print('Reading.....')
    for update in updates:
        cwd = os.getcwd()
        for i,j,y in os.walk(update):
            output = is_espresso(i)
            if output:
                os.chdir(i)
                print('dire:', i)
                calc.directory = cwd + '/' + i
                calc.prefix = output[0:-4]
                try:
                    calc.results = {}
                    calc.read_results()
                    gap, p1, p2 = bandgap(calc)
                    calc.results['gap'] = gap
                    t = calc.get_time()
                    calc.results['time'] = t
                    datas[i] = calc.results
                    results = ana(i, calc)
                    df.loc[len(df)] = results
                except Exception as e:
                    print('='*30, '\n', i, e)
            os.chdir(cwd)
    with open(file, 'wb') as f:
        pickle.dump([datas, df], f)
    print('Finished')

def is_espresso(path):
    '''
    check espresso 
    '''
    dirs = os.listdir(path)
    # print(dirs)
    # flag = True
    for qefile in ['.pwo']:
        flag = False
        for file in dirs:
            if qefile in file:
                return file
        if not flag:
            return False
    # return flag
def grep_valence_configuration(pseudopotential):
    """
    Given a UPF pseudopotential file, find the valence configuration.
    
    Valence configuration:
    nl pn  l   occ       Rcut    Rcut US       E pseu
    3S  1  0  2.00      0.700      1.200    -6.910117

    """
    orbitals = ['S', 'P', 'D', 'F']
    valence = {}
    with open(pseudopotential) as psfile:
        lines = psfile.readlines()
        for i in range(len(lines)):
            if 'valence configuration:' in lines[i].lower():
                j = i + 2
                ob = lines[j].split()[0]
                while ob[1] in orbitals:
                    valence[ob] = lines[j].split()[3]
                    j += 1
                    ob = lines[j].split()[0]
                return valence
    if not valence:
        raise ValueError('Valence configuration missing in {}'.format(pseudopotential))

