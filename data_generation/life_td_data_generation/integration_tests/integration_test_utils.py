from astropy.table import Table
from utils.io import Path, load


def test_load():
    [table] = load(["mes_mass_st"], location=Path().data)

    assert type(table) == type(Table())
