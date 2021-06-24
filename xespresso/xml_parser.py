from operator import pos
import xml.etree.ElementTree as ET
from pprint import pprint
from ase import Atoms, symbols
from ase.atom import Atom
from ase.atoms import default
from ase.units import create_units
from ase.utils import iofunction
import numpy as np
from xespresso.utils import modify_text

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

different_names = {
'CONTROL': {'forces': 'tprnfor', 'stress': 'tstress'}
}
def xml_parser(xmlfile):
    """
    """
    pwi_parameters = {}
    # section_names = ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', 'ATOMIC_SPECIES', 'K_POINTS', 'CELL_PARAMETERS', 'CONSTRAINTS', 'OCCUPATIONS', 'ATOMIC_VELECITIES', 'ATOMIC_FORCES']
    section_names = ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL', ]
    for section_name in section_names:
        pwi_parameters[section_name] = {}
    #
    tree = ET.parse(xmlfile)  
    root = tree.getroot()
    input = root.find('./input')
    parameters = {}
    subs = ['control_variables', 'spin', 'bands', 'basis', 'electron_control', 'ion_control', 'cell_control']
    for sub in subs:
        parameters[sub] = {}
        for key in input.find('./%s'%sub):  
            parameters[sub][key.tag] = key.text
            parameters[sub].update(key.attrib)
    #
    pwi_parameters['CONTROL'].update(parameters['control_variables'])
    pwi_parameters['SYSTEM'].update(parameters['bands'])
    pwi_parameters['SYSTEM'].update(parameters['basis'])
    pwi_parameters['ELECTRONS'].update(parameters['electron_control'])
    #
    for section_name, paras in different_names.items():
        for key, value in paras.items():
            pwi_parameters[section_name][value] = pwi_parameters[section_name].pop(key)
    
    # atomic_species
    atomic_species, pseudopotentials, input_ntyp = get_atomic_species(input.find('./atomic_species'))
    atoms = get_atomic_structure(input.find('./atomic_structure'))
    kpts, koffset = get_kpoint_grid(input.find('./k_points_IBZ'))
    #
    # dftU
    dftU_parameters = get_dftU(input.find('./dft'), input_ntyp)
    pwi_parameters['SYSTEM'].update(dftU_parameters)
    #
    if parameters['spin']['lsda'].upper() == 'TRUE':
        pwi_parameters['SYSTEM'].update({'nspin': 2})
    #
    pwi_parameters = xml2pw(pwi_parameters)
    # Ha to Ry
    for key in ['forc_conv_thr', 'etot_conv_thr']:
        if key in pwi_parameters['CONTROL']:
            pwi_parameters['CONTROL'][key] = pwi_parameters['CONTROL'][key]*2
    for key in ['ecutrho', 'ecutwfc', 'ecutfock']:
        if key in pwi_parameters['SYSTEM']:
            pwi_parameters['SYSTEM'][key] = pwi_parameters['SYSTEM'][key]*2
    #   
    pwi_parameters['INPUT_NTYP']  = input_ntyp
    parameters = {
    'input_data': pwi_parameters,
    'pseudopotentials': pseudopotentials,
    'kpts': kpts,
    'koffset': koffset,
    # 'atomic_species': atomic_species,
    }
    return atoms, parameters
def get_atomic_structure(xml):
    """
    Parse atom.

    Parameters
    """
    nat = xml.attrib['nat']
    alat = xml.attrib['alat']
    atoms = Atoms()
    positions = []
    symbols = []
    all_species = []
    for item in xml.find('./atomic_positions'):
        species = item.attrib['name']
        all_species.append(species)
        symbol = species.split('_')[0].split('-')[0]
        if len(symbol) == 3:
            symbol = symbol[0:2]
        symbols.append(symbol)
        index = item.attrib['index']
        data = item.text.split()
        # These can be fractions and other expressions
        position = [float(data[0])*units['Bohr'], float(data[1])*units['Bohr'], float(data[2])*units['Bohr']]
        positions.append(position)
    positions = np.array(positions)
    # positions = positions**units['Bohr']
    atoms = Atoms(symbols = symbols, positions = positions)
    atoms.info['species'] = all_species
    # cell
    cell = []
    for item in xml.find('./cell'):
        data = item.text.split()
        cell.append([float(data[0]), float(data[1]), float(data[2])])
    cell = np.array(cell)
    cell = cell*units['Bohr']
    atoms.cell = cell
    atoms.pbc = True
    return atoms
def get_atomic_species(xml):
    """
    Parameters
    """
    default = ['mass', 'pseudo_file']
    atomic_species = []
    pseudopotentials = {}
    input_ntyp = {}
    for sub in xml.findall('./species'):
        species = {}
        name = sub.attrib['name']
        species['name'] = name
        for item in sub:
            species[item.tag] = item.text
            if item.tag not in default:
                if item.tag not in input_ntyp:
                    input_ntyp[item.tag] = {}
                input_ntyp[item.tag][name] = float(item.text)
        pseudopotentials[name] = species['pseudo_file']
        atomic_species.append(species)
    return atomic_species, pseudopotentials, input_ntyp
def get_dftU(xml, input_ntyp):
    """
    Parameters
    """
    dftU_parameters = {}
    if not xml.find('./dftU'):
     return dftU_parameters
    
    for item in xml.find('./dftU'):
        if item.tag == 'Hubbard_U':
            dftU_parameters['lda_plus_u'] = True
            if 'Hubbard_U' not in input_ntyp:
                input_ntyp['Hubbard_U'] = {}
            input_ntyp['Hubbard_U'][item.attrib['specie']] = float(item.text)/0.073498810939358 #Ry to eV
        else:
            dftU_parameters[item.tag] = item.text
    # print(input_ntyp)
    return dftU_parameters
def get_kpoint_grid(xml):
    """
    Parameters
    """
    data = xml[0].attrib
    kpts = (int(data['nk1']), int(data['nk2']), int(data['nk3']))
    koffset = (int(data['k1']), int(data['k2']), int(data['k3']))
    return kpts, koffset
def xml2pw(parameters, package = 'PW'):
    import copy
    from xespresso.input_parameters import qe_namespace
    new_parameters = copy.deepcopy(parameters)
    for section, paras in parameters.items():
        for key, value in paras.items():
            if key not in qe_namespace[package][section]:
                del new_parameters[section][key]
            else:
                # print(key, parameters[section][key], qe_namespace[package][section][key][0])
                # print(key, parameters[section][key], qe_namespace[package][section][key][1])
                new_parameters[section][key] = modify_text(parameters[section][key], type=qe_namespace[package][section][key][0])
                # if not compare_value(new_parameters[section][key], qe_namespace[package][section][key][1], qe_namespace[package][section][key][0]):
                    # del new_parameters[section][key]
    return new_parameters


if __name__ == "__main__":
    xmlfile = '/home/xing/ase/xespresso/restart/scf/fe+u/fe+u.xml'
    xmlfile = '/home/xing/ase/xespresso/restart/relax/h2o/h2o.xml'
    xml_data = xml_parser(xmlfile)
    print(xml_data['atoms'])
    print('pseudopotentials: ', xml_data['pseudopotentials'])
    pprint(xml_data['parameters'])
    