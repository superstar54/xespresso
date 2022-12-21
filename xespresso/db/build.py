import sys
import os
from ase.visualize import view
from xespresso import Espresso
from ase.db import connect
from ase.io.espresso import read_espresso_in
from xespresso.xio import read_espresso_input
import subprocess
import os.path
import time


def build_db(update=".", prefix="datas"):
    user = os.getenv("USER")
    print("Reading.....")
    dbfile = "%s.db" % prefix
    with connect(dbfile) as db:
        calc = Espresso()
        cwd = os.getcwd()
        for dire, j, y in os.walk(update):
            prefix = is_espresso(dire)
            if prefix:
                calc.prefix = prefix
            else:
                continue
            print(dire)
            label = dire
            inputfile = dire + "/%s.pwi" % calc.prefix
            calc.directory = os.path.join(cwd, dire)
            calc.prefix = prefix
            print(inputfile)
            atoms, input_data, pseudopotentials, kpts = read_espresso_input(inputfile)
            newcalc = Espresso(
                label=dire,
                prefix=prefix,
                pseudopotentials=pseudopotentials,
                input_data=input_data,
                kpts=kpts,
            )
            converged, fromfile, meg0 = calc.read_convergence()
            calc.results = {}
            try:
                calc.read_results()
                atoms = calc.results["atoms"]
            except Exception as e:
                print("=" * 30, "\n", dire, e)
            # atoms.calc = newcalc
            print("dire: ", dire, converged)
            nid = 0
            for row in db.select(label=label):
                dbid = row.id
                nid += 1
            if nid > 0:
                db.write(atoms, id=dbid, label=label, converged=meg0)
            else:
                db.write(atoms, label=label)
    print("Finished")


def is_espresso(path):
    """
    check espresso
    """
    dirs = os.listdir(path)
    for file in dirs:
        if file[-4:] == ".pwi":
            return file[0:-4]
    return False


# tools
def summary(updates=[], prefix="datas"):
    import pickle
    import pandas as pd

    columns = ["label", "atoms", "cell", "positions", "energy", "forces", "stress"]
    file = "%s.pickle" % prefix
    db = "%s.db" % prefix
    if os.path.exists(file):
        with open(file, "rb") as f:
            datas, df = pickle.load(f)
    else:
        datas = {}
        df = pd.DataFrame(columns=columns)
    calc = Espresso()
    print("Reading.....")
    for update in updates:
        cwd = os.getcwd()
        for i, j, y in os.walk(update):
            output = is_espresso(i)
            if output:
                os.chdir(i)
                print("dire:", i)
                calc.directory = cwd + "/" + i
                calc.prefix = output[0:-4]
                try:
                    calc.results = {}
                    calc.read_results()
                    # gap, p1, p2 = bandgap(calc)
                    # calc.results['gap'] = gap
                    # t = calc.get_time()
                    # calc.results['time'] = t
                    datas[i] = calc.results
                    atoms = calc.results["atoms"]
                    atoms.write(os.path.join(calc.directory, "%s.cif" % calc.prefix))
                    # results = ana(i, calc)
                    # df.loc[len(df)] = results
                except Exception as e:
                    print("=" * 30, "\n", i, e)
            os.chdir(cwd)
    with open(file, "wb") as f:
        pickle.dump([datas, df], f)
    print("Finished")


if __name__ == "__main__":
    summary(".")
