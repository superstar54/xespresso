## xespresso
Quantum ESPRESSO Calculator for Atomic Simulation Environment (ASE).

For the introduction of ASE , please visit https://wiki.fysik.dtu.dk/ase/index.html


### Functions:
* Automatic submit job
* Automatic set up "nscf" calculation
* Plot pdos

### Author
* Xing Wang  <xingwang1991@gmail.com>

### Dependencies

* Python
* ASE

### Installation

Clone this repo. Add it to your PYTHONPATH and PATH. On windows, you can edit the system environment variables.

``` sh
export PYTHONPATH="Your-Location":$PYTHONPATH
export ASE_ESPRESSO_COMMAND="/path/to/PACKAGE.x PARALLEL -in PREFIX.PACKAGEi > PREFIX.PACKAGEo"
```


### Examples

#### Automatic submit job

A example of setting parameters for the queue.

``` python
#!/usr/bin/env python3
from ase.build import bulk
from xespresso import XEspresso
atoms = bulk('Fe', cubic = 'True')
input_data = {
'calculate': 'scf',
'verbosity': 'high', 
'ecutwfc': 40.0,
'ecutrho': 320.0,
'occupations': 'smearing',
'smearing': 'gaussian',
'degauss': 0.01,
'nspin': 2,
}
queue = {'nodes': 1, 
         'ntasks-per-node': 4, 
		 'account': 'dcb', 
		 'partition': 'all', 
		 'time': '0:10:00'}
pseudopotentials = {'Fe': 'Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'}
calc = XEspresso(pseudopotentials = pseudopotentials, 
				 queue = queue,
				 label  = 'scf/fe',
				 input_data = input_data, kpts=(8, 8, 8))
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
print('Energy: ', e)
```

#### Add new species
Some atoms are special:
+ atoms with different starting_magnetization
+ atoms with different U values
+ atoms with special basis set

For example, Fe with spin state AFM:

``` python
atoms = bulk('Fe', cubic=True)
atoms.arrays['species'] = atoms.get_chemical_symbols()
atoms.arrays['species'][0] = 'Fe1'
atoms.arrays['species'][1] = 'Fe2'
```

#### Setting parameters with "(i), i=1,ntyp"
Hubbard, starting_magnetization, starting_charge...

``` python
input_ntyp = {
'starting_magnetization': {'Fe1': 1.0, 'Fe2': -1.0},
'Hubbard_U': {'Fe1': 3.0, 'Fe2': 3.0},
}
```
then add input_ntyp into input_data.
``` python
input_data['input_ntyp'] = input_ntyp,
```

#### Control parallelization levels
To control the number of processors in each group: -ni,
-nk, -nb, -nt, -nd) are used.

``` python
calc = Espresso(pseudopotentials = pseudopotentials, 
                 package = 'pw',
                 parallel = '-nk 2 -nt 4 -nd 144',  # parallel parameters
				 }
```

#### nscf calculation

A example of nscf calculation following the above one.

``` python
calc.read_results()
calc.nscf(queue = queue, kpts = (12, 12, 12))
calc.nscf_calculate()
```

#### calculate dos and pdos

A example of calculating and plotting the pdos from the nscf calculation.

``` python
calc.read_results()
calc.post(queue = queue, package = 'dos', Emin = fe - 30, Emax = fe + 30, DeltaE = 0.1)
calc.post(queue = queue, package = 'projwfc', Emin = fe - 30, Emax = fe + 30, DeltaE = 0.1)
```
<!-- <img src="examples/figs/al-pdos.png" width="500"/> -->

#### calculate work function
``` python
calc.post(queue = queue, package = 'pp', plot_num = 11, fileout = 'potential.cube', iflag = 3, output_format=6)
calc.get_work_function()
```

#### restar from previous calculation
``` python
calc.read_results()
atoms = calc.results['atoms']       
calc.run(atoms = atoms, restart = 1)
```

