"""
This file takes an object from the database and (re-)moves it. E.g. changing a planet into a star object.
"""

import numpy as np
from astropy.table import Table, vstack


def remove_system(main_id, dictionary):
    # doing via main_id instead of object_id because I want to be able to do it before whole db is built
    # remove measurements from tables...

    tables = ["ident", "star_basic", "mes_binary", "mes_sep_ang"]
    for table in tables:
        dictionary[table] = dictionary[table][
            np.where(dictionary[table]["main_id"] != main_id)
        ]

    # look in detail at h_link, (best_h_link) stuff
    # -> adapt parameters of other objects where necessary. might need function below for finding companions

    # sources (via ref) -> can't just remove them, might be used in other parameter.
    # instead need to find the ones that are never used

    # lastly remove object in objects table
    dictionary["objects"] = dictionary["objects"][
        np.where(dictionary["objects"]["main_id"] != main_id)
    ]
    return dictionary


def remove_object(main_id, dictionary):
    # get object type from main id
    object_type = dictionary["objects"]["type"][
        np.where(dictionary["objects"]["main_id"] == main_id)
    ]

    # calls different functions for different object types:
    if object_type == "st":
        dictionary = remove_star(main_id, dictionary)
    elif object_type == "pl":
        dictionary = remove_planet(main_id, dictionary)
    elif object_type == "sy":
        dictionary = remove_system(main_id, dictionary)
    elif object_type == "di":
        dictionary = remove_disk(main_id, dictionary)

    return dictionary


def find_companions(boundary_objects, database_tables):
    for b_object_id in boundary_objects:
        parents = database_tables["best_h_link"]["parent_object_idref"][
            np.where(
                database_tables["best_h_link"]["child_object_idref"]
                == b_object_id
            )
        ]
        children = database_tables["best_h_link"]["child_object_idref"][
            np.where(
                database_tables["best_h_link"]["parent_object_idref"]
                == b_object_id
            )
        ]
        if len(parents) > 1:
            siblings = Table()
            for par in parents:
                par_children = database_tables["best_h_link"][
                    "child_object_idref"
                ][
                    np.where(
                        database_tables["best_h_link"]["parent_object_idref"]
                        == par
                    )
                ]
                # remove b_object_id
                b_object_siblings = par_children[
                    np.where(par_children["child_object_idref"] != b_object_id)
                ]
                b_object_siblings.rename_column(
                    "child_object_idref", "sibling_object_idref"
                )
                siblings = vstack(siblings, b_object_siblings)

        companions = list(parents["parent_object_idref"])
        companions = companions + list(children["child_object_idref"])
        companions = companions + list(siblings["sibling_object_idref"])
    return companions


def manage_boundary_objects(database_tables):
    # find boundary objects
    boundary_objects = database_tables["star_basic"]["object_idref"][
        np.where(database_tables["star_basic"]["dist_st_value"] > 30)
    ]

    # find their companions
    companions = find_companions(boundary_objects, database_tables)

    # assign boundary_flag to companions objects table
    database_tables["objects"]["boundary_flag"] = [
        False for j in range(len(boundary_objects))
    ]
    for companion in companions:
        database_tables["objects"]["boundary_flag"][
            np.where(database_tables["objects"]["object_id"] == companion)
        ] = True

    # remove all occurences of boundary objects
    # this needs to be made nicer with AI
    for table_name in table_names:
        if table_name == "objects":
            colname = "object_id"
        elif table_name == "best_h_link" or "h_link":
            # issue, need to go over two columns !!!
            # colname = 'child_object_idref'
            colname = "parent_object_idref"
        else:
            colanme = "object_idref"
        for b_object_id in boundary_objects:
            database_tables[table_name] = database_tables[table_name][
                np.where(database_tables[table_name][colname] != b_object_id)
            ]

    return database_tables


def test_manage_boundary_objects():
    # data
    database_tables = empty_dict_wit_columns
    # insert three objects of which one is outside of distance cut
    database_tables["sources"]["ref", "provider_name", "source_id"].add_row(
        ["test", "LIFE", 1]
    )

    database_tables["objects"]["type", "ids", "main_id", "object_id"].add_row(
        ["sy", "Boundary AB|Boundary", "Boundary AB", 1]
    )
    database_tables["objects"]["type", "ids", "main_id", "object_id"].add_row(
        ["st", "Boundary A|A", "Boundary A", 2]
    )
    database_tables["objects"]["type", "ids", "main_id", "object_id"].add_row(
        ["st", "Boundary B|B", "Boundary B", 3]
    )

    database_tables["provider"][
        "provider_name", "provider_url", "provider_bibcode", "provider_access"
    ].add_row(["LIFE", "", "test", "2024-12-13"])

    database_tables["ident"]["object_idref", "id", "id_source_idref"].add_row(
        [1, "Boundary AB", 1]
    )
    database_tables["ident"]["object_idref", "id", "id_source_idref"].add_row(
        [1, "Boundary", 1]
    )
    database_tables["ident"]["object_idref", "id", "id_source_idref"].add_row(
        [2, "Boundary A", 1]
    )
    database_tables["ident"]["object_idref", "id", "id_source_idref"].add_row(
        [2, "A", 1]
    )
    database_tables["ident"]["object_idref", "id", "id_source_idref"].add_row(
        [3, "Boundary B", 1]
    )
    database_tables["ident"]["object_idref", "id", "id_source_idref"].add_row(
        [3, "B", 1]
    )

    database_tables["best_h_link"][
        "child_object_idref",
        "parent_object_idref",
        "h_link_source_idref",
        "h_link_ref",
        "membership",
    ].add_row([2, 1, "test", 100])
    database_tables["best_h_link"][
        "child_object_idref",
        "parent_object_idref",
        "h_link_source_idref",
        "h_link_ref",
        "membership",
    ].add_row([3, 1, "test", 100])
    # most likely need to assignd differently as I don't have all the rows listed
    database_tables["star_basic"][
        "object_idref", "main_id", "dist_st_value"
    ].add_row([2, "Boundary A", 29.9])
    database_tables["star_basic"][
        "object_idref", "main_id", "dist_st_value"
    ].add_row([3, "Boundary B", 30.1])

    # add some more data to make sure it gets removed smoothly from the other tables too

    # function
    # database_tables = manage_boundary_objects(database_tables)

    # assert
    # make sure there is a flag companion_out_of_cut
    # make sure the objects and all its occurences gets removed from the db
