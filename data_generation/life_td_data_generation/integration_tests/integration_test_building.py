from astropy.table import Table, unique
from building import create_multi_measurement_star_basic_table, create_best_para_star_basic, best_parameters_ingestion
from utils.io import load, Path
import numpy as np
from sdata import empty_dict, para_columns_dict

def test_best_parameters_ingestion():
    [multimes] = load(["temp_star_basic"], location=Path().additional_data)
    main_id = "2MASS J00282091+2249050"

    print(multimes["plx_value", "plx_qual", "sptype_string"][
              multimes["main_id"] == main_id])
    result = best_parameters_ingestion(
        multimes["main_id","plx_value","plx_qual",
                "plx_err",
                "plx_source_idref",
                "plx_ref"],
        unique(multimes,keys = "main_id"),
        "plx",
        ["plx_value","plx_qual",
                "plx_err",
                "plx_source_idref",
                "plx_ref"],
    )

    print(result["plx_value","plx_qual"][result["main_id"]==main_id])
    # so for a single parameter it works. issue is either in the combination of them or in how I call the ingestion function
    assert False

def test_create_best_para_star_basic():
    #load data
    [multimes] = load(["temp_star_basic"],location=Path().additional_data)
    main_id = "2MASS J00282091+2249050"

    print(multimes["plx_value","plx_qual","sptype_string"][multimes["main_id"]==main_id])

    new_star_basic = create_best_para_star_basic(multimes)
    print(new_star_basic)
    print(new_star_basic["plx_value","plx_qual","sptype_string"][new_star_basic["main_id"] == main_id])

    # hm apparently it run through everything but didn't do output
    assert len(new_star_basic) < len(multimes)
    assert new_star_basic["plx_value"][
        new_star_basic["main_id"] == main_id] == 42.8879
    assert new_star_basic["plx_qual"][
        new_star_basic["main_id"] == main_id] == "A"
    assert new_star_basic["sptype_string"][
        new_star_basic["main_id"] == main_id] == "L6.0"
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







