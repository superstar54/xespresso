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
```


### Examples

#### Automatic submit job

A example of setting parameters for the queue.

``` python
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
				 queue = queue,
				 label  = 'scf/al',
				 input_data = input_data, kpts=(8, 8, 8))
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
print('Energy: ', e)
````

#### nscf calculation

A example of nscf calculation following the above one.

``` python
calc.read_results()
calc.nscf(kpts = (12, 12, 12))
calc.nscf_calculate()
````

#### plot pdos

A example of calculating and plotting the pdos from the nscf calculation.

``` python
calc.read_results()
calc.calc_pdos()
calc.plot_pdos(total = True)
````
<img src="examples/figs/al-pdos.png" width="500"/>