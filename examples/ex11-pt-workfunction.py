from ase.build import fcc111
from xespresso import Espresso
from xespresso.post.pp import EspressoPp

atoms = fcc111("Pt", (1, 1, 4), vacuum=5)
atoms.translate([0, 0, 0.1 - min(atoms.positions[:, 2])])
pseudopotentials = {
    "Pt": "Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF",
}
calc = Espresso(
    label="calculations/scf/pt-fcc111",
    pseudopotentials=pseudopotentials,
    ecutwfc=50,
    occupations="smearing",
    degauss=0.03,
    kpts=(8, 8, 1),
    atoms=atoms,
    debug=True,
)
e = atoms.get_potential_energy()
print("Energy = {0:1.3f} eV".format(e))
# assert np.isclose(e, -606.94121029)
# ===============================================================
# start nscf calculation
fe = calc.get_fermi_level()
print("fermi level: ", fe)
# start nscf calculation
from xespresso.post.nscf import EspressoNscf

nscf = EspressoNscf(
    calc.directory,
    prefix=calc.prefix,
    occupations="tetrahedra",
    kpts=(8, 8, 1),
    debug=True,
)
nscf.run()
# ===============================================================
pp = EspressoPp(
    parent_directory=calc.directory,
    prefix=calc.prefix,
    plot_num=11,
    fileout="potential.cube",
    iflag=3,
    output_format=6,
    debug=True,
)
pp.run()
wf = calc.get_work_function(output="Pt-fcc111-work-function.png")
print(f"work function: {wf} eV")
