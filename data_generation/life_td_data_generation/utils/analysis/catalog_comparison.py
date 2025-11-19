"""
Helper functions for the creation and analysis of the LIFE Target
    Database.
"""

import astropy as ap  # votables
import matplotlib.pyplot as plt
import numpy as np  # arrays

from utils.io import stringtoobject


def testobject_dropout(test_objects, parent_sample, silent=False):
    """
    Find objects in test_objects that are not in parent_sample.

    :param test_objects: Array/column of object identifiers to test
    :param parent_sample: Array/column of parent sample identifiers
    :param silent: If True, suppress print output
    :return: Tuple of (dropout objects, test objects without dropout)
    """
    if len(test_objects) > 0:
        # Find objects in test_objects that are NOT in parent_sample
        drop_out = test_objects[
            np.where(
                np.isin(test_objects, parent_sample, invert=True)
            )
        ]

        # Find objects in test_objects that ARE in parent_sample
        # (i.e., objects that did NOT drop out)
        test_objects_without_dropout = test_objects[
            np.where(
                np.isin(test_objects, parent_sample, invert=False)
            )
        ]

        if not silent:
            print("The following objects are not part of the parent sample: \n")
            print(drop_out)
    else:
        print("test objects sample is empty")
        drop_out = np.array([])  # Return empty array instead of empty list
        test_objects_without_dropout = np.array([])

    return drop_out, test_objects_without_dropout

def type_system(cat_h,cat_o, lists_dict, main_id, name, verbose):
    """
    Process a system-type object and categorize it based on its children.

    This function checks if a system object (type='sy') has non planet or disk
    child objects in the hierarchical link table. If children are found, they
    are added to the 'children' list. If no children are found, the system is added
    to the 'system_without_child' list for tracking.

    This is part of the catalog comparison workflow that analyzes why certain
    objects may not be included in StarCat4.

    :param cat_h: Hierarchical link table containing parent-child relationships
            with columns 'parent_main_id' and 'child_main_id'.
    :type cat_h: astropy.table.table.Table
    :param lists_dict: Mutable dictionary containing categorized object lists.
            Expected keys include 'children' and 'system_without_child'.
    :type lists_dict: dict
    :param main_id: Main identifier of the system object being processed.
    :type main_id: str
    :param name: Name/identifier of the object for tracking purposes.
    :type name: str
    :param verbose: If True, prints diagnostic information about found children.
    :type verbose: bool
    :returns: None. Modifies lists_dict in place.
    :rtype: None
    """
    parent_clause = np.where(cat_h["parent_main_id"] == main_id)

    if len(cat_h[parent_clause]) > 0:
        if verbose:
            reason=f"system object but with found child object, {main_id}"
            #TBD have reason as something to be given back by the function
            # and do the same for the other functions
            print(reason)
            print(
                cat_h["child_main_id"][parent_clause]
            )
        n = 0
        for child in cat_h["child_main_id"][parent_clause]:
            if cat_o["type"][np.where(child == cat_o["main_id"])][0] != "pl":
                lists_dict["children"].append(child)
            else:
                n+=1
        if n==len(cat_h[parent_clause]):
            # only planet type children which we are not interested in here
            lists_dict["system_without_child"].append(name)
    else:
        lists_dict["system_without_child"].append(name)


def type_star(lists_dict, cat_h, cat_o, main_id, name, verbose):
    """
    Process a star-type object and categorize it based on its hierarchical relationships.

    This function analyzes a star object (type='st') to determine its multiplicity status
    by examining its parent-child relationships in the hierarchical link table. Stars are
    categorized as:
    - star_without_parent: Single stars with no parent system
    - single_child: Stars with exactly one stellar sibling
    - binary: Stars in a binary system (2 stellar siblings, no nested systems)
    - higher_order_multiple: Stars in systems with >2 stellar components or nested systems
    - multiple_parents: Stars that appear to have multiple parent systems (unusual case)

    This categorization helps explain why certain objects may not be included in StarCat4,
    particularly those in complex multiple systems that don't meet binary criteria.

    :param lists_dict: Mutable dictionary containing categorized object lists.
            Expected keys include 'star_without_parent', 'multiple_parents', 'single_child',
            'binary', 'higher_order_multiple', and 'siblings'.
    :type lists_dict: dict
    :param cat_h: Hierarchical link table containing parent-child relationships
            with columns 'parent_main_id' and 'child_main_id'.
    :type cat_h: astropy.table.table.Table
    :param cat_o: Objects table containing object types with columns 'main_id' and 'type'
            (where type can be 'st' for star or 'sy' for system).
    :type cat_o: astropy.table.table.Table
    :param main_id: Main identifier of the star object being processed.
    :type main_id: str
    :param name: Name/identifier of the object for tracking purposes.
    :type name: str
    :param verbose: If True, prints diagnostic information about parent relationships,
            siblings, and complex multiplicity cases.
    :type verbose: bool
    :returns: None. Modifies lists_dict in place by appending the object name to
            appropriate categorization lists.
    :rtype: None

    Note:
        The function counts stellar siblings (st_sib) and detects nested hierarchies
        (nestled=True when a sibling is a system). A star with exactly 2 stellar
        siblings and no nesting is classified as a binary; more complex configurations
        are classified as higher_order_multiple.
    """
    # if it has parents:
    if len(cat_h[np.where(cat_h["child_main_id"] == main_id)]) > 0:
        # print('system object but with found child object',main_id)
        # print(cat_h[np.where(cat_h['parent_main_id']==main_id)])
        if len(cat_h[np.where(cat_h["child_main_id"] == main_id)]) > 1:
            if verbose:
                print("star's parents")
                print(cat_h[np.where(cat_h["child_main_id"] == main_id)])
            lists_dict["multiple_parents"].append(name)
        elif len(cat_h[np.where(cat_h["child_main_id"] == main_id)]) == 1:
            st_sib = 0
            sib_name_list = []
            nestled = False
            parent = cat_h["parent_main_id"][
                np.where(cat_h["child_main_id"] == main_id)
            ]
            # wrong code here, parent can't be type st
            for sibling in cat_h["child_main_id"][
                np.where(cat_h["parent_main_id"] == parent)
            ]:
                # print(cat_h[np.where(cat_h['child_main_id'] == sibling)])
                if (
                    cat_o["type"][np.where(cat_o["main_id"] == sibling)][0]
                    == "sy"
                ):
                    nestled = True
                    if verbose:
                        print("sibling is system", sibling)
                elif (
                    cat_o["type"][np.where(cat_o["main_id"] == sibling)][0]
                    == "st"
                ):
                    st_sib += 1
                    sib_name_list.append(sibling)
            lists_dict["siblings"].append(sib_name_list)
            if st_sib == 1 and nestled == False:
                lists_dict["single_child"].append(name)
            if st_sib > 2 or nestled == True:
                lists_dict["higher_order_multiple"].append(name)
            if st_sib == 2 and nestled == False:
                lists_dict["binary"].append(name)
            if verbose and st_sib > 2:
                print(name, "has more than one sibling")
                print(
                    cat_h["child_main_id"][
                        np.where(cat_h["parent_main_id"] == parent)
                    ]
                )
    else:
        lists_dict["star_without_parent"].append(name)


def object_in_db(lists_dict, cat_h, cat_i, cat_o, name, verbose):
    """
    Route an object to appropriate type-specific analysis based on its database type.

    This function serves as a dispatcher that looks up an object by its identifier name,
    determines its type from the objects table, and delegates to the appropriate handler:
    - Systems (type='sy') are processed by type_system()
    - Stars (type='st') are processed by type_star()
    - Other types trigger a warning message

    This is used in the catalog comparison workflow to analyze why objects may not be
    included in StarCat4 by categorizing them based on their hierarchical structure
    and multiplicity.

    :param lists_dict: Mutable dictionary containing categorized object lists that will
        be populated by the type-specific handlers. Expected keys depend on object type
        (e.g., 'children', 'system_without_child', 'star_without_parent', 'binary', etc.).
    :type lists_dict: dict
    :param cat_h: Hierarchical link table containing parent-child relationships
        with columns 'parent_main_id' and 'child_main_id'.
    :type cat_h: astropy.table.table.Table
    :param cat_i: Identifier table used to resolve the given name to a main_id,
        with columns 'id' and 'main_id'.
    :type cat_i: astropy.table.table.Table
    :param cat_o: Objects table containing object types with columns 'main_id' and 'type'
        (where type can be 'st' for star, 'sy' for system, or other types).
    :type cat_o: astropy.table.table.Table
    :param name: Identifier/alias of the object to be analyzed (will be resolved to main_id).
    :type name: str
    :param verbose: If True, prints diagnostic information via the type-specific handlers.
    :type verbose: bool
    :returns: None. Delegates to type_system() or type_star() which modify lists_dict in place.
        Prints a warning for unrecognized object types.
    :rtype: None

    Note:
        This function assumes the identifier exists in cat_i. If the name is not found,
        an IndexError will be raised. Use this function as part of the detail_criteria()
        workflow which filters objects first.
    """
    main_id = cat_i["main_id"][np.where(cat_i["id"] == name)][0]
    if cat_o["type"][np.where(cat_o["main_id"] == main_id)][0] == "sy":
        type_system(cat_h, cat_o, lists_dict, main_id, name, verbose)
    elif cat_o["type"][np.where(cat_o["main_id"] == main_id)][0] == "st":
        type_star(lists_dict, cat_h, cat_o, main_id, name, verbose)
    else:
        print(
            "no sys or st \n",
            name,
            cat_o["main_id", "type"][np.where(cat_o["main_id"] == main_id)][0],
        )


def result(lists_dict, l):
    """
    Print a formatted summary of object categorization results.

    This function displays the categorization results from detail_criteria analysis,
    showing why objects may not be included in StarCat4. Results are organized by
    category with counts and object lists.

    :param lists_dict: Dictionary containing categorized object lists with keys like
            'system_without_child', 'star_without_parent', 'not_found', 'children',
            'multiple_parents', 'higher_order_multiple', 'siblings', 'single_child', 'binary'.
    :type lists_dict: dict
    :param l: Original list of input objects for context in the summary.
    :type l: list
    :returns: None. Prints formatted results to stdout.
    :rtype: None
    """
    # Convert all lists to regular Python strings for clean output
    for key in lists_dict:
        if key == "siblings":
            # Special handling for nested lists
            lists_dict[key] = [[str(item) for item in sublist] for sublist in
                               lists_dict[key]]
        else:
            lists_dict[key] = [str(item) for item in lists_dict[key]]
    print("\n \n Of the", len(l), "objects given:")
    print("Some are not in cat 4 because they are either:")
    print(
        "system_without_child",
        lists_dict["system_without_child"],
        len(lists_dict["system_without_child"]),
    )
    print(
        "star_without_parent",
        lists_dict["star_without_parent"],
        len(lists_dict["star_without_parent"]),
    )
    print("not found", lists_dict["not_found"], "\n")
    print(
        "Some were not included because they are systems and instead one should look at their child objects:"
    )
    print("children", lists_dict["children"], "\n")
    print("some are non trivial binaries:")
    print(
        "multiple parents",
        lists_dict["multiple_parents"],
        len(lists_dict["multiple_parents"]),
    )
    print(
        "higher_order_multiple",
        lists_dict["higher_order_multiple"],
        len(lists_dict["higher_order_multiple"]),
        "\n",
    )
    print("siblings", lists_dict["siblings"])
    if len(lists_dict["single_child"]) > 0:
        print(
            "single_child",
            lists_dict["single_child"],
            len(lists_dict["single_child"]),
        )
    print(
        "And the reminder have companions that don t fit the spectral type requirements"
    )
    print("trivial binary", lists_dict["binary"], len(lists_dict["binary"]))


def detail_criteria(database_tables, l, verbose=True):
    """
    This code scans for reasons, why the given objects would not be included in StarCat4

    """
    print(
        "tbd: file starcat4 analysis where this code is used to give specific reason why a given object was not included into the cat4"
    )
    # to do:
    # include get main id search
    # this code could use printing the corresponding objects why it was not included.
    # write tests
    cat_i = database_tables["ident"]
    cat_o = database_tables["objects"]
    cat_h = database_tables["best_h_link"]  # -> use best_h_link
    lists_dict = {
        "children": [],
        "not_found": [],
        "siblings": [],
        "multiple_parents": [],  # multiple components in wds
        "single_child": [],
        "binary": [],
        "higher_order_multiple": [],
        "star_without_parent": [],
        "system_without_child": [],
    }
    for name in l:
        if len(cat_i["main_id"][np.where(cat_i["id"] == name)]) > 0:
            object_in_db(lists_dict, cat_h, cat_i, cat_o, name, verbose)
        else:
            # print('object not found',name)
            lists_dict["not_found"].append(name)
    result(lists_dict, l)
    return


def object_contained(stars, cat, details=False):
    """
    Checks which of the entries in stars is also in cat.

    This function checks which of the star identifiers in stars are not
    present in cat. If details is True then the amount and individual objects
    which are in stars but not cat are printed. Returns only star identifiers
    in stars that are also present in cat.

    :param stars: numpy string array containing star identifiers
    :type stars:
    :param cat: numpy string array containing star identifiers
    :type cat:
    :param bool details: Defaults to False.
    :returns: numpy string array containing star identifiers
        which are also present in cat
    :rtype:
    """
    lost = stars[np.where(np.invert(np.in1d(stars, cat)))]
    if len(lost) > 0:
        print("This criterion was not met by:", len(lost))
        if details:
            print(lost)
    stars = stars[np.where(np.invert(np.in1d(stars, lost)))]
    return stars


def compare_catalogs(
    cat1, cat2, cat1_idname, cat2_idname, cat1_paranames, cat2_paranames
):
    """
    This function compares two catalogs.

    :param cat1:  to be compared to cat2
    :type cat1: astropy.table.table.Table
    :param cat2:  to be compared to cat1
    :type cat2: astropy.table.table.Table
    :param str cat1_idname: column name of cat 1 column containing object
        identifiers
    :param str cat2_idname: column name of cat 2 column containing object
        identifiers
    :param cat1_paranames: column names of cat 1 columns containing
        measurements to be compared to cat 2
    :type cat1_paranames: list(str)
    :param cat2_paranames: column names of cat 2 columns containing measurements
        to be compared to cat 1. Order of column names needs to be same as in cat1_parameters
    :type cat2_paranames: list(str)
    """
    print("lenght cat 1:", len(cat1))
    print("lenght cat 2:", len(cat2))
    common = cat1[np.where(np.in1d(cat1[cat1_idname], cat2[cat2_idname]))]
    print("number of common objects:", len(common))
    print("number of objects in cat 1 but not cat 2:", len(cat1) - len(common))
    print("number of objects in cat 2 but not cat 1:", len(cat2) - len(common))
    common2 = cat2[np.where(np.in1d(cat2[cat2_idname], cat1[cat1_idname]))]
    difference = ap.table.join(
        common, common2, keys_left=cat1_idname, keys_right=cat2_idname
    )
    # still need to do some detail work on colnames. I can call all paranames of 1 with 1 in end
    for i in range(len(cat1_paranames)):
        print(f"Distribution of parameter {cat1_paranames[i]}")
        if cat1_paranames[i] != cat2_paranames[i]:
            print(f" respectively {cat2_paranames[i]}")

        # removing objects with masked or nan values:
        print("tbd remove objects with masked or nan values")
        # histogram with three subplots (cat1,cat2,common)
        fig, ((ax1, ax2, ax3)) = plt.subplots(
            1,
            3,
            figsize=[16, 5],
            gridspec_kw={"hspace": 0, "wspace": 0},
            sharey=True,
        )

        ax1.set_xlabel(cat1_paranames[i])
        ax1.set_ylabel("Number of objects")
        ax1.hist(cat1[cat1_paranames[i]], histtype="bar", log=True)
        # ax1.set_xscale('log')

        ax2.set_xlabel(cat2_paranames[i])
        ax2.hist(cat2[cat2_paranames[i]], histtype="bar", log=True)

        ax3.set_xlabel(cat1_paranames[i])
        ax3.hist(common[cat1_paranames[i]], histtype="bar", log=True)
        plt.show()

        plt.figure
        plt.xlabel(cat1_paranames[i])
        plt.ylabel("Number of objects")
        plt.hist(
            cat1[cat1_paranames[i]],
            histtype="bar",
            log=True,
            color="C0",
            alpha=0.8,
            label="cat1",
        )
        plt.hist(
            cat2[cat2_paranames[i]],
            histtype="bar",
            log=True,
            color="C1",
            alpha=0.8,
            label="cat2",
        )
        plt.hist(
            common[cat1_paranames[i]],
            histtype="bar",
            log=True,
            color="C2",
            alpha=0.8,
            label="common",
        )
        plt.legend()
        plt.show()
        print("cat1, cat2 and common objects with values of cat1")

        # mean
        print(
            "cat 1 mean and std:",
            np.mean(cat1[cat1_paranames[i]]),
            np.std(cat1[cat1_paranames[i]]),
        )
        # std (The standard deviation is the square root of the average of
        #   the squared deviations from the mean, i.e., std = sqrt(mean(x)), where x = abs(a - a.mean())**2)
        print(
            "cat 2 mean and std:",
            np.mean(cat2[cat2_paranames[i]]),
            np.std(cat2[cat2_paranames[i]]),
        )
        print(
            "common objects cat 1 values mean and std:",
            np.mean(common[cat1_paranames[i]]),
            np.std(common[cat1_paranames[i]]),
        )
        print(
            "common objects cat 2 values mean and std:",
            np.mean(common2[cat2_paranames[i]]),
            np.std(common2[cat2_paranames[i]]),
        )
        print(f"Distribution of difference of parameter {cat1_paranames[i]}")
        if cat1_paranames[i] != cat2_paranames[i]:
            print(f" respectively {cat2_paranames[i]}")
            difference[f"diff_{cat1_paranames[i]}"] = (
                difference[cat1_paranames[i]] - difference[cat2_paranames[i]]
            )
        else:
            difference[f"diff_{cat1_paranames[i]}"] = (
                difference[cat1_paranames[i] + "_1"]
                - difference[cat2_paranames[i] + "_2"]
            )
        print(
            "differences of common objects values mean and std:",
            np.mean(difference[f"diff_{cat1_paranames[i]}"]),
            np.std(difference[f"diff_{cat1_paranames[i]}"]),
        )
        plt.figure
        plt.xlabel(f"diff_{cat1_paranames[i]}")
        plt.ylabel("Number of objects")
        plt.hist(
            difference[f"diff_{cat1_paranames[i]}"], histtype="bar", log=True
        )
        plt.show()
    return


def create_common(cat1, cat2):
    """
    Takes two catalogs, removes objects that are not in the other one.

    :param cat1: First catalog.
    :type cat1: astropy.table.table.Table
    :param cat2: Second catalog.
    :type cat2: astropy.table.table.Table
    :returns: Touple of cat1 containing only objects that are also in
        cat2 and cat2 containing only objects that are also in cat1.
    :rtype: list(astropy.table.table.Table)
    """
    common_stars = np.intersect1d(list(cat1["main_id"]), list(cat2["main_id"]))
    common_cat1 = cat1[np.where(np.in1d(cat1["main_id"], common_stars))]
    common_cat2 = cat2[np.where(np.in1d(cat2["main_id"], common_stars))]
    common_cat2 = stringtoobject(common_cat2, number=1000)
    common_cat1 = stringtoobject(common_cat1, number=1000)
    return common_cat1, common_cat2
