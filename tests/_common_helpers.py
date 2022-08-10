"""This module is only to be used together with tests
"""
import os


def set_envs():
    """Detect pw"""
    import os
    cwd = os.getcwd()
    os.environ['ESPRESSO_PSEUDO'] = os.path.join(cwd, "datas/pseudo")

def bulk_h():
    from ase.build import bulk
    h = bulk('Fe', cubic=True)
    h.set_chemical_symbols(['H', 'H'])
    return h

