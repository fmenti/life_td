"""Generate WDS provider data for the database."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np  # arrays
from astropy.table import Table, join, unique, vstack
from provider.assign_quality_funcs import assign_quality
from provider.utils import (
    create_provider_table,
    create_sources_table,
    distance_cut,
    ids_from_ident,
    query,
)
from sdata import empty_dict
from utils.io import load, save


def _print_test_matches(
    label: str, values: Sequence[Any], candidates: Table, column: str
) -> None:
    """Print matching test objects for a given candidate column."""
    if len(values) == 0:
        return
    print(
        label,
        values[np.where(np.isin(values, candidates[column]))],
    )


def load_wds_helptab() -> Table:
    """Load the cached WDS helper table.

    :returns: Cached WDS helper table.
    :rtype: astropy.table.Table
    """
    print(" loading...")
    [wds_helptab] = load(["wds_helptab"])
    # currently temp=True not giving same result because
    # wds['system_main_id'][j] are '' and not masked
    for col in ["system_main_id", "primary_main_id", "secondary_main_id"]:
        wds_helptab[col][np.where(wds_helptab[col] == "")] = np.ma.masked

    # tbd: add provider_access of last query
    return wds_helptab


def assign_names(wds_helptab: Table) -> Table:
    """Assign WDS system/component names.

    :param wds_helptab: Raw WDS helper table.
    :type wds_helptab: astropy.table.Table
    :returns: Helper table with name columns populated.
    :rtype: astropy.table.Table
    """
    for j in range(len(wds_helptab)):
        wds_name = wds_helptab["wds_name"][j]
        wds_comp = wds_helptab["wds_comp"][j]

        if wds_comp == "":  # trivial binaries
            wds_helptab["system_name"][j] = f"WDS J{wds_name}AB"
            # AB added since SIMBAD calls trivial binary system AB too.
            wds_helptab["primary"][j] = f"WDS J{wds_name}A"
            wds_helptab["secondary"][j] = f"WDS J{wds_name}B"
            continue

        # Higher-order multiples.
        wds_helptab["system_name"][j] = f"WDS J{wds_name}{wds_comp}"
        if len(wds_comp) == 2:
            wds_helptab["primary"][j] = f"WDS J{wds_name}{wds_comp[0]}"
            wds_helptab["secondary"][j] = f"WDS J{wds_name}{wds_comp[1]}"
        elif "," in wds_comp:
            comp_a, comp_b = wds_comp.split(",", maxsplit=1)
            wds_helptab["primary"][j] = f"WDS J{wds_name}{comp_a}"
            wds_helptab["secondary"][j] = f"WDS J{wds_name}{comp_b}"
        else:
            print("not sure how to handle:", wds_comp)

    return wds_helptab


def look_at_test_objects_after_name_assignment(
    test_objects: np.ndarray,
    wds_helptab: Table,
) -> None:
    """Print which test objects match the derived WDS name columns."""
    _print_test_matches(
        "in wds as system_name", test_objects, wds_helptab, "system_name"
    )
    _print_test_matches(
        "in wds as primary", test_objects, wds_helptab, "primary"
    )
    _print_test_matches(
        "in wds as secondary", test_objects, wds_helptab, "secondary"
    )


def look_at_test_objects_after_wds_creation(
    test_objects: np.ndarray,
    wds_helptab: Table,
) -> None:
    """Print which test objects match the main-id columns after creation."""
    if len(test_objects) == 0:
        return
    print("test objects matching main-id columns after creation:")
    print(wds_helptab["system_main_id", "primary_main_id", "secondary_main_id"])
    _print_test_matches(
        "in wds as system_main_id", test_objects, wds_helptab, "system_main_id"
    )
    _print_test_matches(
        "in wds as primary_main_id",
        test_objects,
        wds_helptab,
        "primary_main_id",
    )
    _print_test_matches(
        "in wds as secondary_main_id",
        test_objects,
        wds_helptab,
        "secondary_main_id",
    )


def create_wds_helptab(
    adql_query: list[str],
    test_objects: np.ndarray,
    wds: dict[str, Table],
) -> Table:
    """Query WDS and build the helper table.

    :param adql_query: Query list; the first query is used here.
    :type adql_query: list[str]
    :param test_objects: Objects used for debug prints.
    :type test_objects: numpy.ndarray
    :param wds: Provider table dictionary.
    :type wds: dict[str, astropy.table.Table]
    :returns: WDS helper table.
    :rtype: astropy.table.Table
    """
    print(" querying VizieR for WDS...")
    wds_helptab = query(wds["provider"]["provider_url"][0], adql_query[0])

    # Match WDS objects with SIMBAD-derived objects to enforce the distance cut.
    for col in ["sim_wds_id", "system_name", "primary", "secondary"]:
        wds_helptab[col] = wds_helptab["wds_name"].astype(object)

    wds_helptab = assign_names(wds_helptab)
    look_at_test_objects_after_name_assignment(test_objects, wds_helptab)
    # an alternative would be to query simbad for the main id and then cut
    # by distance
    # this however takes way longer as it joins 150'000 elements
    #    wds=fetch_main_id(wds,colname='wds_full_name',name='main_id',oid=False)
    #    wds=distance_cut(wds,colname='wds_full_name',main_id=True)
    print(" performing distance cut...")
    wds_system_cut = distance_cut(
        wds_helptab, colname="system_name", main_id=False
    )
    wds_system_cut.rename_column("main_id", "system_main_id")

    wds_primary_cut = distance_cut(
        wds_helptab, colname="primary", main_id=False
    )
    wds_secondary_cut = distance_cut(
        wds_helptab, colname="secondary", main_id=False
    )

    [sim_h_link] = load(["sim_h_link"])
    # joining parent object
    wds_primary_cut = join(
        wds_primary_cut,
        sim_h_link["main_id", "parent_main_id"],
        keys="main_id",
        join_type="left",
    )
    wds_primary_cut.rename_columns(
        ["main_id", "parent_main_id"],
        ["primary_main_id", "system_main_id"],
    )

    wds_secondary_cut = join(
        wds_secondary_cut,
        sim_h_link["main_id", "parent_main_id"],
        keys="main_id",
        join_type="left",
    )
    wds_secondary_cut.rename_columns(
        ["main_id", "parent_main_id"],
        ["secondary_main_id", "system_main_id"],
    )
    # here some empty ones when child is known in simbad but no parent. in
    # this case would I want to assign system_name in system main_id? do it
    # later
    wds_helptab = vstack([wds_system_cut, wds_primary_cut])
    wds_helptab = vstack([wds_helptab, wds_secondary_cut])
    look_at_test_objects_after_wds_creation(test_objects, wds_helptab)
    save([wds_helptab], ["wds_helptab"])
    return wds_helptab


def create_wds_helpertable(
    temp: bool,
    test_objects: Sequence[str],
) -> tuple[Table, dict[str, Table]]:
    """Create the WDS helper table and provider metadata.

    :param temp: If ``True``, load the cached helper table.
    :type temp: bool
    :param test_objects: Debug object identifiers.
    :type test_objects: collections.abc.Sequence[str]
    :returns: Helper table and provider dictionary.
    :rtype: tuple[astropy.table.Table, dict[str, astropy.table.Table]]
    """
    wds: dict[str, Table] = empty_dict.copy()
    wds["provider"] = create_provider_table(
        "WDS",
        "http://tapvizier.u-strasbg.fr/TAPVizieR/tap",
        "2001AJ....122.3466M",
    )

    adql_query = [
        """SELECT WDS  as wds_name,
                  Comp as wds_comp,
                  sep1 as wds_sep1,
                  sep2 as wds_sep2,
                  Obs1 as wds_obs1,
                  Obs2 as wds_obs2
           FROM "B/wds/wds" """
    ]

    print("Creating ", wds["provider"]["provider_name"][0], " tables ...")
    test_objects_arr = np.array(test_objects)

    if temp:
        wds_helptab = load_wds_helptab()
    else:
        wds_helptab = create_wds_helptab(adql_query, test_objects_arr, wds)

    wds_helptab["system_main_id"] = wds_helptab["system_main_id"].astype(object)
    wds_helptab["system_name"] = wds_helptab["system_name"].astype(object)
    return wds_helptab, wds


def create_ident_and_h_link_table(
    wds_helptab: Table,
    wds: dict[str, Table],
    test_objects: Sequence[str],
) -> tuple[Table, Table]:
    """Create identifier and hierarchical-link tables.

    :param wds_helptab: WDS helper table.
    :type wds_helptab: astropy.table.Table
    :param wds: Provider dictionary.
    :type wds: dict[str, astropy.table.Table]
    :param test_objects: Debug object identifiers.
    :type test_objects: collections.abc.Sequence[str]
    :returns: Identifier and hierarchical-link tables.
    :rtype: tuple[astropy.table.Table, astropy.table.Table]
    """
    wds_ident = Table(
        names=["main_id", "id"], dtype=[object, object], masked=True
    )
    wds_h_link = Table(
        names=["main_id", "parent_main_id"], dtype=[object, object]
    )

    id_cols_1 = [
        "system_name",
        "system_main_id",
        "system_main_id",
        "primary",
        "primary_main_id",
        "primary_main_id",
        "secondary",
        "secondary_main_id",
        "secondary_main_id",
    ]
    id_cols_2 = [
        "system_name",
        "system_main_id",
        "system_name",
        "primary",
        "primary_main_id",
        "primary",
        "secondary",
        "secondary_main_id",
        "secondary",
    ]

    empty = Table(names=["main_id"], dtype=[object], masked=True)
    for id1, id2 in zip(id_cols_1, id_cols_2):
        temp = empty.copy()
        temp["main_id"] = wds_helptab[id1].astype(object)
        temp["id"] = wds_helptab[id2].astype(object)
        wds_ident = vstack([wds_ident, temp])

    relation_cols_1 = [
        "primary",
        "primary",
        "primary_main_id",
        "primary_main_id",
        "secondary",
        "secondary",
        "secondary_main_id",
        "secondary_main_id",
    ]
    relation_cols_2 = [
        "system_name",
        "system_main_id",
        "system_name",
        "system_main_id",
        "system_name",
        "system_main_id",
        "system_name",
        "system_main_id",
    ]
    for id1, id2 in zip(relation_cols_1, relation_cols_2):
        temp = empty.copy()
        temp["main_id"] = wds_helptab[id1].astype(object)
        temp["parent_main_id"] = wds_helptab[id2].astype(object)
        wds_h_link = vstack([wds_h_link, temp])
    # delete all rows containing masked entries
    wds_ident.remove_rows(wds_ident["main_id"].mask.nonzero()[0])
    wds_ident.remove_rows(wds_ident["id"].mask.nonzero()[0])
    wds_h_link.remove_rows(wds_h_link["main_id"].mask.nonzero()[0])
    wds_h_link.remove_rows(wds_h_link["parent_main_id"].mask.nonzero()[0])

    wds_ident = unique(wds_ident)
    wds_h_link = unique(wds_h_link)

    # delete entries where id instead of main_id used
    not_identical_rows_id = wds_ident["id"][
        np.where(wds_ident["main_id"] != wds_ident["id"])
    ]
    remove = np.isin(wds_ident["main_id"], not_identical_rows_id)
    wds_ident.remove_rows(remove)

    # for h_link replacing instead of deleting because there can be cases where
    # the information I need is only available this way
    # e.g. simbad query on system name results in system main_id. simbad query
    # on primary gives primary_main_id. but
    # simbad does not have those two objects connected through h_link. therefore
    # I only have [primary,system_main_id] but
    # would want to have [primary_main_id,system_main_id]

    # where h_link main_id not in ident_main_id
    not_main_id = np.invert(
        np.isin(wds_h_link["main_id"], wds_ident["main_id"])
    )

    if len(test_objects) > 0:
        test_objects = np.array(test_objects)
        print(
            "number of test objects that are in h_link main_id \n",
            test_objects[
                np.where(np.isin(test_objects, wds_h_link["main_id"]))
            ],
        )
        print(
            "number of test objects that are in h_link parent_main_id \n",
            test_objects[
                np.where(np.isin(test_objects, wds_h_link["parent_main_id"]))
            ],
        )
        print(
            "number of test objects that are in main_id of ident table \n",
            test_objects[np.where(np.isin(test_objects, wds_ident["main_id"]))],
        )
        print(
            "number of test objects that are in main_id of h_link but not "
            "ident table \n",
            test_objects[
                np.where(
                    np.isin(test_objects, wds_h_link["main_id"][not_main_id])
                )
            ],
        )

    for temp in wds_h_link["main_id"][not_main_id]:
        wds_h_link["main_id"][np.where(wds_h_link["main_id"] == temp)] = (
            wds_ident["main_id"][np.where(wds_ident["id"] == temp)]
        )

    not_parent_main_id = np.invert(
        np.isin(wds_h_link["parent_main_id"], wds_ident["main_id"])
    )

    for temp in wds_h_link["parent_main_id"][not_parent_main_id]:
        if len(wds_ident["main_id"][np.where(wds_ident["id"] == temp)]) == 1:
            wds_h_link["parent_main_id"][
                np.where(wds_h_link["parent_main_id"] == temp)
            ] = wds_ident["main_id"][np.where(wds_ident["id"] == temp)]

    wds_h_link = unique(wds_h_link)

    wds_ident["id_ref"] = [
        wds["provider"]["provider_bibcode"][0] for _ in range(len(wds_ident))
    ]
    wds_h_link["h_link_ref"] = [
        wds["provider"]["provider_bibcode"][0] for _ in range(len(wds_h_link))
    ]
    return wds_ident, wds_h_link


def create_objects_table(
    wds_helptab: Table,
    wds: dict[str, Table],
    test_objects: Sequence[str],
) -> Table:
    """Create the object table.

    Objects start as systems, then objects without children are reclassified
    as stars.

    :param wds_helptab: WDS helper table.
    :type wds_helptab: astropy.table.Table
    :param wds: Provider dictionary.
    :type wds: dict[str, astropy.table.Table]
    :param test_objects: Debug object identifiers.
    :type test_objects: collections.abc.Sequence[str]
    :returns: Object table.
    :rtype: astropy.table.Table
    """
    wds_objects = Table(names=["main_id", "ids"], dtype=[object, object])
    wds_objects = ids_from_ident(wds["ident"]["main_id", "id"], wds_objects)
    wds_objects["type"] = ["sy" for _ in range(len(wds_objects))]

    has_no_children = np.invert(
        np.isin(wds_objects["main_id"], wds["h_link"]["parent_main_id"])
    )
    wds_objects["type"][has_no_children] = "st"

    if len(test_objects) > 0:
        test_objects = np.array(test_objects)
        print(
            "number of test objects that are in objects main_id \n",
            test_objects[
                np.where(np.isin(test_objects, wds_objects["main_id"]))
            ],
        )
    return wds_objects


def create_mes_binary_table(
    wds_helptab: Table,
    wds: dict[str, Table],
    test_objects: Sequence[str],
) -> Table:
    """Create the binary measurement table.

    :param wds_helptab: WDS helper table.
    :type wds_helptab: astropy.table.Table
    :param wds: Provider dictionary.
    :type wds: dict[str, astropy.table.Table]
    :param test_objects: Debug object identifiers.
    :type test_objects: collections.abc.Sequence[str]
    :returns: Binary measurement table.
    :rtype: astropy.table.Table
    """
    wds_mes_binary = wds["objects"]["main_id", "type"]
    wds_mes_binary.rename_column("type", "binary_flag")
    wds_mes_binary["binary_flag"] = wds_mes_binary["binary_flag"].astype(object)
    wds_mes_binary["binary_flag"] = ["True" for _ in range(len(wds_mes_binary))]
    wds_mes_binary["binary_ref"] = [
        wds["provider"]["provider_bibcode"][0]
        for _ in range(len(wds_mes_binary))
    ]
    wds_mes_binary = assign_quality(
        wds_mes_binary,
        "binary_qual",
        special_mode="wds_binary",
    )

    if len(test_objects) > 0:
        test_objects = np.array(test_objects)
        print(
            "number of test objects that are in mes_binary main_id \n",
            test_objects[
                np.where(np.isin(test_objects, wds_mes_binary["main_id"]))
            ],
        )
    return wds_mes_binary


def create_mes_sep_ang_table(
    wds_helptab: Table,
    wds: dict[str, Table],
    test_objects: Sequence[str],
) -> Table:
    """Create the angular-separation measurement table.

    :param wds_helptab: WDS helper table.
    :type wds_helptab: astropy.table.Table
    :param wds: Provider dictionary.
    :type wds: dict[str, astropy.table.Table]
    :param test_objects: Debug object identifiers.
    :type test_objects: collections.abc.Sequence[str]
    :returns: Angular-separation measurement table.
    :rtype: astropy.table.Table
    """
    wds_mes_sep_ang0 = join(
        wds_helptab[
            "system_name", "wds_sep1", "wds_obs1", "wds_sep2", "wds_obs2"
        ],
        wds["ident"]["main_id", "id"],
        keys_left="system_name",
        keys_right="id",
    )

    wds_mes_sep_ang1 = wds_mes_sep_ang0["main_id", "wds_sep1", "wds_obs1"]
    wds_mes_sep_ang1.rename_columns(
        ["wds_sep1", "wds_obs1"],
        ["sep_ang_value", "sep_ang_obs_date"],
    )
    wds_mes_sep_ang1 = assign_quality(
        wds_mes_sep_ang1,
        "sep_ang_qual",
        special_mode="wds_sep1",
    )

    wds_mes_sep_ang2 = wds_mes_sep_ang0["main_id", "wds_sep2", "wds_obs2"]
    wds_mes_sep_ang2.rename_columns(
        ["wds_sep2", "wds_obs2"],
        ["sep_ang_value", "sep_ang_obs_date"],
    )
    wds_mes_sep_ang2 = assign_quality(
        wds_mes_sep_ang2,
        "sep_ang_qual",
        special_mode="wds_sep2",
    )

    wds_mes_sep_ang = vstack([wds_mes_sep_ang1, wds_mes_sep_ang2])
    wds_mes_sep_ang["sep_ang_ref"] = [
        wds["provider"]["provider_bibcode"][0]
        for _ in range(len(wds_mes_sep_ang))
    ]

    wds_mes_sep_ang.remove_columns(
        wds_mes_sep_ang["sep_ang_value"].mask.nonzero()[0]
    )
    # uniqueness where obs date not known
    if len(wds_mes_sep_ang["sep_ang_obs_date"].mask.nonzero()[0]) > 0:
        unique_unknown_obs_date = unique(
            wds_mes_sep_ang[
                np.where(wds_mes_sep_ang["sep_ang_obs_date"].mask.nonzero()[0])
            ],
            keys=["main_id", "sep_ang_value"],
        )
        unique_known_obs_date = unique(
            wds_mes_sep_ang[
                np.where(
                    np.invert(
                        wds_mes_sep_ang["sep_ang_obs_date"].mask.nonzero()[0]
                    )
                )
            ],
            keys=["main_id", "sep_ang_value", "sep_ang_obs_date"],
        )
        wds_mes_sep_ang = vstack(
            [unique_unknown_obs_date, unique_known_obs_date]
        )
    else:
        wds_mes_sep_ang = unique(wds_mes_sep_ang)

    if len(test_objects) > 0:
        test_objects = np.array(test_objects)
        print(
            "number of test objects that are in mes_sep_ang main_id \n",
            test_objects[
                np.where(np.isin(test_objects, wds_mes_sep_ang["main_id"]))
            ],
        )

    return wds_mes_sep_ang


def create_wds_sources_table(wds: dict[str, Table]) -> Table:
    """Create the sources table.

    :param wds: Provider dictionary.
    :type wds: dict[str, astropy.table.Table]
    :returns: Sources table.
    :rtype: astropy.table.Table
    """
    tables = [wds["provider"], wds["ident"]]
    ref_columns = [["provider_bibcode"], ["id_ref"]]
    return create_sources_table(
        tables,
        ref_columns,
        wds["provider"]["provider_name"][0],
    )


def provider_wds(
    temp: bool = False, test_objects: Sequence[str] | None = None
) -> dict[str, Table]:
    """Build all WDS provider tables.

    :param temp: If ``True``, load cached helper data instead of querying.
    :type temp: bool
    :param test_objects: Debug object identifiers. Defaults to an empty list.
    :type test_objects: collections.abc.Sequence[str] | None
    :returns: Provider table dictionary.
    :rtype: dict[str, astropy.table.Table]
    """
    if test_objects is None:
        test_objects = []

    wds_helptab, wds = create_wds_helpertable(temp, test_objects)
    wds["ident"], wds["h_link"] = create_ident_and_h_link_table(
        wds_helptab, wds, test_objects
    )
    wds["objects"] = create_objects_table(wds_helptab, wds, test_objects)
    wds["mes_binary"] = create_mes_binary_table(wds_helptab, wds, test_objects)
    wds["mes_sep_ang"] = create_mes_sep_ang_table(
        wds_helptab, wds, test_objects
    )
    wds["sources"] = create_wds_sources_table(wds)

    save(list(wds.values()), ["wds_" + element for element in list(wds.keys())])
    return wds
