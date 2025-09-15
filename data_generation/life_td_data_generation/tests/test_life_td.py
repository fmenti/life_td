import numpy as np
from astropy.table import Table, setdiff

# self created modules
from sdata import empty_dict_wit_columns
from utils.io import string_to_object_whole_dict, save
from life_td import load_cat

def test_load_cat():
    # data
    sptype = np.array(["dM3.51", "dM3:", "dM5.0"])
    main_id = np.array(["G 227-48B", "* mu. Dra C", "FBS 1415+456"])
    example_table = Table((main_id, sptype), names=("main_id", "sptype_string"))
    dictionary = empty_dict_wit_columns.copy()
    dictionary["star_basic"] = example_table
    dictionary = string_to_object_whole_dict(dictionary, 400)

    provider_name = "test_cat"
    save(
        list(dictionary.values()),
        [provider_name + "_" + element for element in list(dictionary.keys())],
    )

    # function
    loaded_dict = load_cat(provider_name)

    # assert
    for tablename in dictionary.keys():
        if len(dictionary[tablename]) > 0:  # change to more once this works
            for columnname in dictionary[tablename].colnames:
                assert (
                    len(setdiff(dictionary[tablename], loaded_dict[tablename]))
                    == 0
                )


# def test_partial_create():
#     partial_create(distance_cut_in_pc,create=[])
#     assert
