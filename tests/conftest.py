import pytest


@pytest.fixture
def bulk_h():
    from ase.build import bulk

    h = bulk("Fe", cubic=True)
    h.set_chemical_symbols(["H", "H"])
    return h
