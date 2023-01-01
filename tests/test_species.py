from _common_helpers import set_envs


def test_spin(bulk_h):
    from xespresso import Espresso
    import numpy as np

    set_envs()
    atoms = bulk_h
    atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
    atoms.arrays["species"][0] = "H"
    atoms.arrays["species"][1] = "H1"
    print(atoms.arrays["species"])
    input_ntyp = {
        "starting_magnetization": {
            "H": 1.0,
            "H1": -1.0,
        }
    }
    pseudopotentials = {
        "H": "H.pbe-rrkjus_psl.1.0.0.UPF",
        "H1": "H.pbe-rrkjus_psl.1.0.0.UPF",
    }
    calc = Espresso(
        pseudopotentials=pseudopotentials,
        label="calculations/scf/h-afm",
        ecutwfc=30,
        occupations="smearing",
        degauss=0.02,
        nspin=2,
        input_data={"input_ntyp": input_ntyp},
        kpts=(4, 4, 4),
        debug=True,
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print("Energy  {0:1.5f}".format(e))
    assert np.isclose(e, -27.57523)


def test_dft_u():
    """_Example using 1) AFM 2) DFT+U"""
    from ase.build import bulk
    from xespresso import Espresso
    import numpy as np

    set_envs()
    atoms = bulk("Fe", cubic=True)
    atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))
    atoms.arrays["species"][1] = "Fe1"
    input_ntyp = {
        "starting_magnetization": {
            "Fe": 0.5,
            "Fe1": -0.5,
        },
        "Hubbard_U": {"Fe": 4.3, "Fe1": 4.3},
    }
    input_data = {
        "ecutwfc": 30.0,
        "occupations": "smearing",
        "degauss": 0.03,
        "nspin": 2,
        "lda_plus_u": True,
        "input_ntyp": input_ntyp,
        #
        "mixing_beta": 0.3,
        "conv_thr": 1.0e-8,
    }
    pseudopotentials = {
        "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
        "Fe1": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
    }
    calc = Espresso(
        pseudopotentials=pseudopotentials,
        label="calculations/scf/fe-afm",
        input_data=input_data,
        kpts=(4, 4, 4),
        debug=True,
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print("Energy {0:1.5f}".format(e))
    assert np.isclose(e, -6675.33857)
