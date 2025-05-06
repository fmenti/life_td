from utils.io import load, Path
from astropy.table import Table
import numpy as np

def test_load():
    [table]=load(['mes_mass_st'],location=Path().data)
    
    assert type(table)==type(Table())
    






