import os
from os import path
import json
import pickle
import operator as op
import warnings
from ase import constraints
import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
from ase.constraints import FixAtoms, FixCartesian
from ase.data import chemical_symbols, atomic_numbers
from ase.io.espresso import read_espresso_in, construct_namelist, grep_valence, SSSP_VALENCE, read_fortran_namelist
from ase.units import create_units
from pprint import pprint
from xespresso.utils import check_type


def write_espresso_in(filename, atoms, input_data={}, pseudopotentials=None,
                      kspacing=None, kpts=None, koffset=(0, 0, 0),
                      crystal_coordinates=False, **kwargs):
    """
    Create an input file for pw.x.

    Modified the read_espresso_in function from ase.io.espresso
    Parameters:
        Please have a look at the ase.io.espresso module in ASE
    
    Implemented:

    - support parameters with ntyp
    - support atomic species
    """

    # Convert to a namelist to make working with parameters much easier
    # Note that the name ``input_data`` is chosen to prevent clash with
    # ``parameters`` in Calculator objects
    input_parameters = construct_namelist(input_data, **kwargs)

    # Construct input file into this
    pwi = []
    # atomic_species
    atomic_species_str, species_info, total_valence = build_atomic_species_str(atoms, input_parameters, pseudopotentials)
    # sections
    section_str, input_parameters = build_section_str(atoms, species_info, input_data, input_parameters)
    pwi.extend(section_str)
    # Pseudopotentials
    pwi.extend(atomic_species_str)
    # KPOINTS - add a MP grid as required
    pwi.extend(build_kpts_str(atoms, kspacing, kpts, koffset))
    # CELL block, if required
    pwi.extend(build_cell_str(atoms, input_parameters))
    # Positions - already constructed, but must appear after namelist
    engine_str = pwi.copy()
    pwi.extend(build_atomic_positions_str(atoms, crystal_coordinates))
    # write file or not
    if filename:
        with open(filename, 'w') as fd:
            fd.write(''.join(pwi))
    # return section_str, atomic_species_str, kpts_str, cell_str, atomic_positions_str
    return engine_str

def build_section_str(atoms, species_info, input_data, input_parameters):
    '''
    '''
    # Add computed parameters
    # different magnetisms means different types
    input_parameters['system']['ntyp'] = len(species_info)
    input_parameters['system']['nat'] = len(atoms)

    #
    if 'INPUT_NTYP' in input_data:
        for key, value in input_data['INPUT_NTYP'].items():
            for species in value:
                if species in species_info:
                    mag_str = '{0}({1})'.format(key, species_info[species]['index'])
                    input_parameters['system'][mag_str] = value[species]
    if 'hubbard_v' in input_data:
        for key, value in input_data['hubbard_v'].items():
                    mag_str = 'Hubbard_V{0}'.format(key)
                    input_parameters['system'][mag_str] = value
    # Use cell as given or fit to a specific ibrav
    if 'ibrav' in input_parameters['system']:
        ibrav = input_parameters['system']['ibrav']
        if ibrav != 0:
            celldm = cell_to_ibrav(atoms.cell, ibrav)
            regen_cell = ibrav_to_cell(celldm)[1]
            if not np.allclose(atoms.cell, regen_cell):
                warnings.warn('Input cell does not match requested ibrav'
                              '{} != {}'.format(regen_cell, atoms.cell))
            input_parameters['system'].update(celldm)
    else:
        # Just use standard cell block
        input_parameters['system']['ibrav'] = 0
    # Assume sections are ordered (taken care of in namelist construction)
    # and that repr converts to a QE readable representation (except bools)
    section_str = []
    for section in input_parameters:
        section_str.append('&{0}\n'.format(section.upper()))
        for key, value in input_parameters[section].items():
            if value is True:
                section_str.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                section_str.append('   {0:16} = .false.\n'.format(key))
            else:
                # repr format to get quotes around strings
                section_str.append('   {0:16} = {1!r:}\n'.format(key, value))
        section_str.append('/\n')  # terminate section
    section_str.append('\n')
    return section_str, input_parameters

def build_atomic_species_str(atoms, input_parameters, pseudopotentials):
    '''
    '''
    # Deal with pseudopotentials
    # Look in all possible locations for the pseudos and try to figure
    # out the number of valence electrons
    pseudo_dirs = []
    if 'pseudo_dir' in input_parameters['control']:
        pseudo_dirs.append(input_parameters['control']['pseudo_dir'])
    if 'ESPRESSO_PSEUDO' in os.environ:
        pseudo_dirs.append(os.environ['ESPRESSO_PSEUDO'])
    pseudo_dirs.append(path.expanduser('~/espresso/pseudo/'))

    # Species info holds the information on the pseudopotential and
    # associated for each element
    if 'species' not in atoms.info:
        # print('no species')
        atoms.info['species'] = atoms.get_chemical_symbols()
    if pseudopotentials is None:
        pseudopotentials = {}
    species_info = {}
    atomic_species_str = ['ATOMIC_SPECIES\n']
    ntyp = 0
    # Convert atoms into species.
    for i in range(len(atoms)):
        species = atoms.info['species'][i]
        if species not in species_info:
            ntyp += 1
            species_info[species] = {}
            species_info[species]['index'] = ntyp
            species_info[species]['mass'] = atoms[i].mass
            species_info[species]['element'] = atoms[i].symbol
            species_info[species]['count'] = 1
        else:
            species_info[species]['count'] += 1
    total_valence = 0
    for species in species_info:
        pseudo = pseudopotentials.get(species, '{}_dummy.UPF'.format(species))
        for pseudo_dir in pseudo_dirs:
            if path.exists(path.join(pseudo_dir, pseudo)):
                valence = grep_valence(path.join(pseudo_dir, pseudo))
                break
        else:  # not found in a file
            valence = SSSP_VALENCE[atomic_numbers[species_info[species]['element']]] 
        species_info[species]['pseudo'] = pseudo,
        species_info[species]['valence'] = valence
        total_valence += valence*species_info[species]['count']
        atomic_species_str.append(
            '{species} {mass} {pseudo}\n'.format(
                species=species, mass=species_info[species]['mass'],
                pseudo=pseudo))
    atomic_species_str.append('\n')
    return atomic_species_str, species_info, total_valence
              

def build_cell_str(atoms, input_parameters):
    '''
    '''
    # CELL block, if required
    if input_parameters['SYSTEM']['ibrav'] == 0:
        cell_str = ['CELL_PARAMETERS angstrom\n']
        cell_str.append('{cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
                   '{cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
                   '{cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}\n'
                   ''.format(cell=atoms.cell))
        cell_str.append('\n')
    return cell_str
def build_kpts_str(atoms, kspacing, kpts, koffset):
    '''
    '''
    # KPOINTS - add a MP grid as required
    if kspacing is not None:
        kgrid = kspacing_to_grid(atoms, kspacing)
    elif kpts is not None:
        if isinstance(kpts, dict) and 'path' not in kpts:
            kgrid, shift = kpts2sizeandoffsets(atoms=atoms, **kpts)
            koffset = []
            for i, x in enumerate(shift):
                assert x == 0 or abs(x * kgrid[i] - 0.5) < 1e-14
                koffset.append(0 if x == 0 else 1)
        else:
            kgrid = kpts
    else:
        kgrid = "gamma"

    # True and False work here and will get converted by ':d' format
    if isinstance(koffset, int):
        koffset = (koffset, ) * 3

    # BandPath object or bandpath-as-dictionary:
    if isinstance(kgrid, dict) or hasattr(kgrid, 'kpts'):
        kpts_str = ['K_POINTS crystal_b\n']
        assert hasattr(kgrid, 'path') or 'path' in kgrid
        kgrid = kpts2ndarray(kgrid, atoms=atoms)
        kpts_str.append('%s\n' % len(kgrid))
        for k in kgrid:
            kpts_str.append('{k[0]:.14f} {k[1]:.14f} {k[2]:.14f} 0\n'.format(k=k))
        kpts_str.append('\n')
    elif isinstance(kgrid, str) and (kgrid == "gamma"):
        kpts_str = ['K_POINTS gamma\n']
        kpts_str.append('\n')
    else:
        kpts_str = ['K_POINTS automatic\n']
        kpts_str.append('{0[0]} {0[1]} {0[2]}  {1[0]:d} {1[1]:d} {1[2]:d}\n'
                   ''.format(kgrid, koffset))
        kpts_str.append('\n')
    return kpts_str

def build_atomic_positions_str(atoms, crystal_coordinates):
    '''

    '''
    import sys
    if 'species' not in atoms.info:
        atoms.info['species'] = atoms.get_chemical_symbols()
    else:
        if len(atoms.info['species']) != len(atoms):
            sys.exit('Species is wrong, please check!')
    # Convert ase constraints to QE constraints
    # Nx3 array of force multipliers matches what QE uses
    # Do this early so it is available when constructing the atoms card
    constraint_mask = np.ones((len(atoms), 3), dtype='int')
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            constraint_mask[constraint.index] = 0
        elif isinstance(constraint, FixCartesian):
            constraint_mask[constraint.a] = constraint.mask
        else:
            warnings.warn('Ignored unknown constraint {}'.format(constraint))
    # Positions - already constructed, but must appear after namelist
    if crystal_coordinates:
        atomic_positions_str = ['ATOMIC_POSITIONS crystal\n']
    else:
        atomic_positions_str = ['ATOMIC_POSITIONS angstrom\n']
    for i in range(len(atoms)):
        atom = atoms[i]
        species = atoms.info['species'][i]
        # only inclued mask if something is fixed
        if not all(constraint_mask[atom.index]):
            mask = ' {mask[0]} {mask[1]} {mask[2]}'.format(
                mask=constraint_mask[atom.index])
        else:
            mask = ''
        if crystal_coordinates:
            coords = [atom.a, atom.b, atom.c]
        else:
            coords = atom.position
        atomic_positions_str.append(
            '{species} '
            '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f} '
            '{mask}\n'.format(species=species, coords=coords, mask=mask))
    return atomic_positions_str

def write_neb_in(filename, images, climbing_images = [], path_data = {}, 
                 input_data={}, pseudopotentials=None,
                 kspacing=None, kpts=None, koffset=(0, 0, 0),
                 crystal_coordinates=False, **kwargs):
    """
    Create an input file for neb.x.

    Parameters
    ----------
    filename: str
        A file to write the input file to.
    images: list
        A list of atomistic configuration including the first, last 
        and intermedinate positions.
    climbing_images: list
        Indices of the images which the Climbing-Image procedure apply
    path_data: dict
        A flat or nested dictionary with path parameters for neb.x
    """
    nebi = ['BEGIN\n']
    package_parameters = {
            'PATH': {'string_method', 'restart_mode', 'nstep_path', 'num_of_images', 
            'opt_scheme', 'CI_scheme', 'first_last_opt', 'minimum_image', 'temp_req', 
            'ds', 'k_max', 'k_min', 'path_thr', 'use_masses', 'use_freezing', 
            'lfcpopt', 'fcp_mu', 'fcp_tot_charge_first', 'fcp_tot_charge_last', },
        }
    path_str = ['\nBEGIN_PATH_INPUT\n']
    for section, parameters in package_parameters.items():
        path_str.append('&%s\n '%section)
        for key, value in path_data.items():
            if key in parameters:
                if value is True:
                    path_str.append('   {0:16} = .true.\n'.format(key))
                elif value is False:
                    path_str.append('   {0:16} = .false.\n'.format(key))
                else:
                    path_str.append('   {0:16} = {1!r:}\n'.format(key, value))
        path_str.append('/ \n')
    if climbing_images:
        path_str.append('CLIMBING_IMAGES\n')
        climbing_images_str = '   '
        for x in climbing_images: climbing_images_str += ' %s'%x
        path_str.append(climbing_images_str)
        path_str.append('\n')
    path_str.append('END_PATH_INPUT\n')
    nebi.extend(path_str)
    
    # BEGIN_ENGINE_INPUT
    engine_str = ['\nBEGIN_ENGINE_INPUT\n']
    pwi = write_espresso_in(None, images[0], input_data, pseudopotentials, 
                            kspacing, kpts, koffset, crystal_coordinates, **kwargs)
    engine_str.extend(pwi)
    
    # POSITION
    atomic_positions_str = ['\nBEGIN_POSITIONS\n']
    for i in range(len(images)):
        atoms = images[i]
        if i == 0:
            atomic_positions_str.append('\nFIRST_IMAGE\n')
        elif i == len(images) - 1:
            atomic_positions_str.append('\nLAST_IMAGE\n')
        else:
            atomic_positions_str.append('\nINTERMEDIATE_IMAGE\n')
        atomic_positions_str.extend(build_atomic_positions_str(atoms, crystal_coordinates))
    atomic_positions_str.append('END_POSITIONS\n')
    #
    engine_str.extend(atomic_positions_str)
    engine_str.append('END_ENGINE_INPUT\n')
    nebi.extend(engine_str)
    nebi.append('END\n')
    with open(filename, 'w') as fd:
        fd.write(''.join(nebi))


def get_atomic_species(pwo):
    '''
    '''
    atomic_species = []
    with open(pwo) as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            if 'number of atoms/cell' in line:
                nat = int(line.split()[-1])
            elif 'positions (alat units)' in line:
                atomic_species = [at_line.split()[1]
                    for at_line in lines[idx + 1:idx + 1 + nat]]
                break
    return atomic_species

def read_espresso_input(fileobj, neb = False):
    """Parse a Quantum ESPRESSO input files, '.in', '.pwi'.
    """
    atoms = read_espresso_in(fileobj)
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'rU')
    # parse namelist section and extract remaining lines
    data, card_lines = read_fortran_namelist(fileobj)
    input_data = {}
    for sec, paras in data.items():
        for key, value in paras.items():
            input_data[key] = value
    # get number of type
    ntyp = data['system']['ntyp']
    pseudopotentials = {}
    trimmed_lines = (line for line in card_lines
                     if line.strip() and not line[0] == '#')
    for line in trimmed_lines:
        if line.strip().startswith('ATOMIC_SPECIES'):
            for i in range(ntyp):
                species = next(trimmed_lines).split()
                pseudopotentials[species[0]] = species[2]
        if line.strip().startswith('K_POINTS'):
            kpts = next(trimmed_lines)
    # calc = Espresso(pseudopotentials=pseudopotentials, 
    #                 input_data = input_data,
    #                 kpts=kpts)
    return atoms, input_data, pseudopotentials, kpts, constraints

def sort_qe_input(parameters, package = 'PW'):
    """
    """
    from xespresso.input_parameters import qe_namespace, default_parameters
    import copy
    pw_parameters = qe_namespace[package]
    if 'input_data' not in parameters:
        parameters['input_data'] = {}
    sorted_parameters = copy.deepcopy(parameters)
    unuse_parameters = {}
    # section_names = ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'ATOMIC_SPECIES', 'K_POINTS', 'CELL_PARAMETERS', 'CONSTRAINTS', 'OCCUPATIONS', 'ATOMIC_VELECITIES', 'ATOMIC_FORCES']
    section_names = ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'INPUT_NTYP']
    for section in section_names:
        if section not in sorted_parameters['input_data']: sorted_parameters['input_data'][section] = {}
    for key, value in parameters.items():
        if key in ['pseudopotentials', 'kpts', 'kspacing', 'koffset', 'input_data', 'climbing_images', 'path_data', 'crystal_coordinates']:
            continue
        flag = False
        for section in section_names:
            if key.upper() == section and isinstance(value, dict):
                sorted_parameters['input_data'][section].update(value)
                flag = True
                del sorted_parameters[key]
                break
            if section.upper() == 'INPUT_NTYP': continue
            if key in pw_parameters[section]:
                sorted_parameters['input_data'][section][key] = value
                flag = True
                del sorted_parameters[key]
                break
        if not flag:
            unuse_parameters[key] = value
            del sorted_parameters[key]
    input_data = copy.deepcopy(sorted_parameters['input_data'])
    for key, value in input_data.items():
        flag = False
        if key.upper() in section_names:
            sorted_parameters['input_data'][key.upper()] = sorted_parameters['input_data'].pop(key)
            continue
        for section in section_names:
            if section.upper() == 'INPUT_NTYP':
                continue
            if key in pw_parameters[section]:
                sorted_parameters['input_data'][section][key] = value
                flag = True
                del sorted_parameters['input_data'][key]
                break
        if not flag:
            unuse_parameters[key] = value
            del sorted_parameters['input_data'][key]
    return sorted_parameters, unuse_parameters
    #
def check_qe_input(input_parameters, package = 'PW'):
    from xespresso.input_parameters import qe_namespace, default_parameters
    from xespresso.utils import check_type
    pw_parameters = qe_namespace[package]
    for section, parameters in input_parameters.items():
        if section == 'INPUT_NTYP':
            for key, subparas in parameters.items():
                for value in subparas.values():
                    check_type(key, value, pw_parameters)
        else:
            for key, value in parameters.items():
                check_type(key, value, pw_parameters)

def read_espresso_asei(fileobj, package = 'PW'):
    """Parse Quantum ESPRESSO input parameters
    """
    with open(fileobj, 'rb') as fp:
        atoms, parameters = pickle.load(fp)
    if package == 'PW':
        sorted_parameters, unuse_parameters = sort_qe_input(parameters)
        check_qe_input(sorted_parameters['input_data'], package)
    else:
        sorted_parameters = parameters
    return atoms, sorted_parameters
def write_espresso_asei(fileobj, atoms, parameters):
    """save Quantum ESPRESSO input parameters
    """
    with open(fileobj, 'wb') as fp:
        pickle.dump([atoms, parameters], fp)

def get_atomic_constraints(pwo, n_atoms):
    """
    Read constraints from QE output
    """
    constraints = []
    with open(pwo, 'r') as f:
        lines = f.readlines()
        indexes = []
        for i, line in enumerate(lines):
            if 'ATOMIC_POSITIONS' in line:
                indexes.append(i)
        if len(indexes) > 0:
            image_index = indexes[0] + 1
            for line in lines[image_index:image_index + n_atoms + 1]:
                split_line = line.split()
                if len(split_line) > 4:
                    constraint = (float(split_line[4]),
                                  float(split_line[5]),
                                  float(split_line[6]))
                else:
                    constraint = (None, None, None)
                constraints.append(constraint)
    constraints = get_constraint(constraints)
    return constraints

def get_constraint(constraint_idx):
    """
    Map constraints from QE output to FixAtoms or FixCartesian constraint
    """
    if not np.any(constraint_idx):
        return None

    a = [a for a, c in enumerate(constraint_idx) if np.all(c is not None)]
    mask = [[(ic + 1) % 2 for ic in c] for c in constraint_idx
            if np.all(c is not None)]

    if np.all(np.array(mask)) == 1:
        constraint = FixAtoms(a)
    else:
        constraint = FixCartesian(a, mask)
    return constraint
