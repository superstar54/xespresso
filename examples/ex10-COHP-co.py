from ase.build import molecule
from xespresso import Espresso
from xespresso.cohp import COHP

# build CO molecule
atoms = molecule("CO")
atoms.center(5)
atoms.pbc = [True, True, True]
pseudopotentials = {
    "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
}
queue = None

# start scf calculation
calc = Espresso(
    label="scf/co",
    pseudopotentials=pseudopotentials,
    wf_collect=True,  # for lobster run, make sure adding this
    outdir=".",  # for lobster run, make sure adding this
    ecutwfc=30,
    occupations="smearing",
    degauss=0.03,
    kpts=(1, 1, 1),
    debug=True,
)
atoms.calc = calc
e = atoms.get_potential_energy()
print("Energy = {0:1.3f} eV".format(e))

# start COHP analysis
cohp = COHP(
    directory="scf/co",  # same as label in scf calculation
    prefix="co",  # same as scf calculation
    # the bond of interest, multiple bond can be set as [[1,2],[3,5]]
    index=[[1, 2]],
    command="your lobster path",  # path of lobster program
    # other arguments can also be parsed, such as follows
    # enter the energetic window in eV (relative to the Fermi level)
    COHPStartEnergy=-10,
    COHPEndEnergy=5
    # specify the basis functions per element manually
    # basisFunctions = ['Fe 3s 3p 3d 4s','S 3s 3p']
)
cohp.write_input()
cohp.run()

# read lobster results and plot
cohp.read_cohp()
cohp.plot_coop()
cohp.plot_cohp()

# read integrated lobster results
cohp.read_icohp()
cohp.icobi
cohp.icohp
cohp.icoop
