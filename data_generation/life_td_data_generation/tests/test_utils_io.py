import numpy as np
from astropy.table import Table
from sdata import empty_dict_wit_columns
from utils.io import (
    objecttostring,
    save,
    stringtoobject,
    string_to_object_whole_dict
)

def test_objecttostring():
    # data
    cat = Table(
        data=[
            ["dM3.51", "dM3:", "dM5.0"],
            ["G 227-48B", "* mu. Dra C", "FBS 1415+456"],
        ],
        names=["main_id", "sptype_string"],
        dtype=[object, object],
    )

    # function
    cat = objecttostring(cat)

    assert cat["main_id"].dtype == "<U6"
    assert cat["sptype_string"].dtype == "<U12"


def test_save():
    # data
    cat = Table(
        data=[
            ["dM3.51", "dM3:", "dM5.0"],
            ["G 227-48B", "* mu. Dra C", "FBS 1415+456"],
        ],
        names=["main_id", "sptype_string"],
        dtype=[object, object],
    )
    # function
    save([cat], "test")

    # assert
    assert cat["main_id"].dtype == object
    # original catalog does get changed


def test_stringtoobject():
    # data
    sptype = np.array(["dM3.51", "dM3:", "dM5.0"])
    main_id = np.array(["G 227-48B", "* mu. Dra C", "FBS 1415+456"])
    example_table = Table((main_id, sptype), names=("main_id", "sptype_string"))

    # function
    example_table = stringtoobject(example_table, number=100)

    # assert
    assert example_table["main_id"].dtype == object


def test_string_to_object_whole_dict():
    # data
    sptype = np.array(["dM3.51", "dM3:", "dM5.0"])
    main_id = np.array(["G 227-48B", "* mu. Dra C", "FBS 1415+456"])
    example_table = Table((main_id, sptype), names=("main_id", "sptype_string"))
    dictionary = empty_dict_wit_columns.copy()

    dictionary["star_basic"] = example_table
    # function
    dictionary = string_to_object_whole_dict(dictionary)

    # assert
    assert dictionary["star_basic"]["main_id"].dtype == object
    for tablename in dictionary.keys():
        if tablename == "star_basic":  # change to more once this works
            for columnname in dictionary[tablename].colnames:
                # no idea why there is a type error here, when there is none above
                assert dictionary[tablename][columnname].dtype == object
