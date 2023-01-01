from _common_helpers import set_envs


def test_hp():
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
    #
    from xespresso.post.hp import EspressoHp

    parameters_hp = {
        "nq1": 1,
        "nq2": 1,
        "nq3": 1,
        "equiv_type": 0,
        "conv_thr_chi": 1e-5,
        "perturb_only_atom": {"1": True},  # key(i)
    }
    hp = EspressoHp(
        parent_directory="calculations/scf/fe-afm", prefix="fe-afm", **parameters_hp
    )
    hp.run()
