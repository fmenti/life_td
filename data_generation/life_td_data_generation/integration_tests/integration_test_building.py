from astropy.table import Table
from building import create_multi_measurement_star_basic_table, create_best_para_star_basic
from utils.io import load, Path
import numpy as np
from sdata import empty_dict

def test_create_best_para_star_basic():
    #load data
    [multimes] = load(["temp_star_basic"],location=Path().additional_data)

    new_star_basic = create_best_para_star_basic(multimes)
    print(new_star_basic)

    # hm apparently it run through everything but didn't do output
    assert len(new_star_basic) < len(multimes)
    assert False



def test_create_multi_measurement_star_basic_table():
    # load data
    [star_basic] = load(["temp_star_basic"], location=Path().additional_data)
    [objects] = load(["objects"], location=Path().data)

    cat = empty_dict
    cat["star_basic"] = star_basic
    cat["objects"] = objects

    # not quite correct, should merge provider star basics instead of loading this

    multimes = create_multi_measurement_star_basic_table(cat)
    print(multimes)



    assert False







