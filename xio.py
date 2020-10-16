import os
import operator as op
import warnings
from collections import OrderedDict
from os import path

import numpy as np

from ase.atoms import Atoms
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)
from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
from ase.dft.kpoints import kpoint_convert
from ase.constraints import FixAtoms, FixCartesian
from ase.data import chemical_symbols, atomic_numbers
from ase.units import create_units

from ase.io.espresso import construct_namelist, grep_valence, SSSP_VALENCE


def write_espresso_in(filename, atoms, input_data=None, pseudopotentials=None,
                      kspacing=None, kpts=None, koffset=(0, 0, 0),
                      crystal_coordinates=False, **kwargs):
    """
    Create an input file for pw.x.

    Use set_initial_magnetic_moments to turn on spin, if ispin is set to 2
    with no magnetic moments, they will all be set to 0.0. Magnetic moments
    will be converted to the QE units (fraction of valence electrons) using
    any pseudopotential files found, or a best guess for the number of
    valence electrons.

    Units are not converted for any other input data, so use Quantum ESPRESSO
    units (Usually Ry or atomic units).

    Keys with a dimension (e.g. Hubbard_U(1)) will be incorporated as-is
    so the `i` should be made to match the output.

    Implemented features:

    - Conversion of :class:`ase.constraints.FixAtoms` and
                    :class:`ase.constraints.FixCartesian`.
    - `starting_magnetization` derived from the `mgmoms` and pseudopotentials
      (searches default paths for pseudo files.)
    - Automatic assignment of options to their correct sections.
    - Interpretation of ibrav (cell must exactly match the vectors defined
      in the QE docs).

    Not implemented:

    - Lists of k-points
    - Other constraints
    - Hubbard parameters
    - Validation of the argument types for input
    - Validation of required options
    - Reorientation for ibrav settings
    - Noncollinear magnetism

    Parameters
    ----------
    fd: file
        A file like object to write the input file to.
    atoms: Atoms
        A single atomistic configuration to write to `fd`.
    input_data: dict
        A flat or nested dictionary with input parameters for pw.x
    pseudopotentials: dict
        A filename for each atomic species, e.g.
        {'O': 'O.pbe-rrkjus.UPF', 'H': 'H.pbe-rrkjus.UPF'}.
        A dummy name will be used if none are given.
    kspacing: float
        Generate a grid of k-points with this as the minimum distance,
        in A^-1 between them in reciprocal space. If set to None, kpts
        will be used instead.
    kpts: (int, int, int) or dict
        If kpts is a tuple (or list) of 3 integers, it is interpreted
        as the dimensions of a Monkhorst-Pack grid.
        If ``kpts`` is set to ``None``, only the Γ-point will be included
        and QE will use routines optimized for Γ-point-only calculations.
        Compared to Γ-point-only calculations without this optimization
        (i.e. with ``kpts=(1, 1, 1)``), the memory and CPU requirements
        are typically reduced by half.
        If kpts is a dict, it will either be interpreted as a path
        in the Brillouin zone (*) if it contains the 'path' keyword,
        otherwise it is converted to a Monkhorst-Pack grid (**).
        (*) see ase.dft.kpoints.bandpath
        (**) see ase.calculators.calculator.kpts2sizeandoffsets
    koffset: (int, int, int)
        Offset of kpoints in each direction. Must be 0 (no offset) or
        1 (half grid offset). Setting to True is equivalent to (1, 1, 1).
    crystal_coordinates: bool
        Whether the atomic positions should be written to the QE input file in
        absolute (False, default) or relative (crystal) coordinates (True).

    """

    # Convert to a namelist to make working with parameters much easier
    # Note that the name ``input_data`` is chosen to prevent clash with
    # ``parameters`` in Calculator objects
    input_parameters = construct_namelist(input_data, **kwargs)

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
    if 'species' not in atoms.arrays:
        print('no species')
        atoms.arrays['species'] = atoms.get_chemical_symbols()
    if pseudopotentials is None:
        pseudopotentials = {}
    species_info = {}
    atomic_species_str = []
    ntyp = 0
    for i in range(len(atoms)):
        species = atoms.arrays['species'][i]
        if species not in species_info:
            ntyp += 1
            species_info[species] = {}
            species_info[species]['index'] = ntyp
            species_info[species]['mass'] = atoms[i].mass
    if 'input_ntyp' in input_data:
        for key, value in input_data['input_ntyp'].items():
            for species in value:
                mag_str = '{0}({1})'.format(key, species_info[species]['index'])
                input_parameters['system'][mag_str] = value[species]
    for species in species_info:
        pseudo = pseudopotentials.get(species, '{}_dummy.UPF'.format(species))
        for pseudo_dir in pseudo_dirs:
            if path.exists(path.join(pseudo_dir, pseudo)):
                valence = grep_valence(path.join(pseudo_dir, pseudo))
                break
        else:  # not found in a file
            valence = SSSP_VALENCE[atomic_numbers[species]] 
        species_info[species]['pseudo'] = pseudo,
        species_info[species]['valence'] = valence
        atomic_species_str.append(
            '{species} {mass} {pseudo}\n'.format(
                species=species, mass=species_info[species]['mass'],
                pseudo=pseudo))
    # Convert atoms into species.
    # Each different magnetic moment needs to be a separate type even with
    # the same pseudopotential (e.g. an up and a down for AFM).
    # if any magmom are > 0 or nspin == 2 then use species labels.
    # Rememeber: magnetisation uses 1 based indexes
    atomic_positions_str = []


    for i in range(len(atoms)):
        atom = atoms[i]
        species = atoms.arrays['species'][i]
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

    # Add computed parameters
    # different magnetisms means different types
    input_parameters['system']['ntyp'] = len(species_info)
    input_parameters['system']['nat'] = len(atoms)

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

    # Construct input file into this
    pwi = []

    # Assume sections are ordered (taken care of in namelist construction)
    # and that repr converts to a QE readable representation (except bools)
    for section in input_parameters:
        pwi.append('&{0}\n'.format(section.upper()))
        for key, value in input_parameters[section].items():
            if value is True:
                pwi.append('   {0:16} = .true.\n'.format(key))
            elif value is False:
                pwi.append('   {0:16} = .false.\n'.format(key))
            else:
                # repr format to get quotes around strings
                pwi.append('   {0:16} = {1!r:}\n'.format(key, value))
        pwi.append('/\n')  # terminate section
    pwi.append('\n')

    # Pseudopotentials
    pwi.append('ATOMIC_SPECIES\n')
    pwi.extend(atomic_species_str)
    pwi.append('\n')

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
        pwi.append('K_POINTS crystal_b\n')
        assert hasattr(kgrid, 'path') or 'path' in kgrid
        kgrid = kpts2ndarray(kgrid, atoms=atoms)
        pwi.append('%s\n' % len(kgrid))
        for k in kgrid:
            pwi.append('{k[0]:.14f} {k[1]:.14f} {k[2]:.14f} 0\n'.format(k=k))
        pwi.append('\n')
    elif isinstance(kgrid, str) and (kgrid == "gamma"):
        pwi.append('K_POINTS gamma\n')
        pwi.append('\n')
    else:
        pwi.append('K_POINTS automatic\n')
        pwi.append('{0[0]} {0[1]} {0[2]}  {1[0]:d} {1[1]:d} {1[2]:d}\n'
                   ''.format(kgrid, koffset))
        pwi.append('\n')

    # CELL block, if required
    if input_parameters['SYSTEM']['ibrav'] == 0:
        pwi.append('CELL_PARAMETERS angstrom\n')
        pwi.append('{cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
                   '{cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
                   '{cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}\n'
                   ''.format(cell=atoms.cell))
        pwi.append('\n')

    # Positions - already constructed, but must appear after namelist
    if crystal_coordinates:
        pwi.append('ATOMIC_POSITIONS crystal\n')
    else:
        pwi.append('ATOMIC_POSITIONS angstrom\n')
    pwi.extend(atomic_positions_str)
    pwi.append('\n')
    # DONE!
    with open(filename, 'w') as fd:
        fd.write(''.join(pwi))