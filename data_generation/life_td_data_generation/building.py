"""Combines the data from the individual data providers."""

from itertools import islice

import numpy as np
from astropy.table import (
    Column,
    MaskedColumn,
    Row,
    Table,
    column,
    join,
    unique,
    vstack,
)
from provider.utils import nullvalues, replace_value
from sdata import empty_dict, empty_dict_wit_columns, paras_dict
from utils.io import Path, save


def idsjoin(cat: Table, column_ids1: str, column_ids2: str) -> Table:
    """
    Merges two identifier columns into one, removing duplicates.

    - Merge identifiers across both columns per row.
    - Remove duplicate identifiers within the same row.
    - Return the result in a row-wise manner, where all identifiers for that
      row are sorted.

    :param cat: Astropy Table containing two identifier columns.
    :type cat: Table
    :param column_ids1: Name of the first identifier column.
    :type column_ids1: str
    :param column_ids2: Name of the second identifier column.
    :type column_ids2: str
    :return: Table with a unified 'ids' column containing merged, unique
        identifiers.
    :rtype: Table
    """
    # Step 1: Replace masked/empty values in both columns with empty strings
    ids1 = (
        cat[column_ids1].filled("")
        if isinstance(cat[column_ids1], MaskedColumn)
        else cat[column_ids1]
    )
    ids2 = (
        cat[column_ids2].filled("")
        if isinstance(cat[column_ids2], MaskedColumn)
        else cat[column_ids2]
    )

    # Step 2: Vectorized merging of identifiers using set operations
    merged_ids: list[str] = []
    for val1, val2 in zip(ids1, ids2):
        # Split strings on '|'
        ids1_list = val1.split("|") if val1 not in (None, "") else []
        ids2_list = val2.split("|") if val2 not in (None, "") else []

        # Merge and remove duplicates/empty strings
        unique_ids = set(ids1_list + ids2_list) - {""}
        merged_ids.append("|".join(sorted(unique_ids)))

    # Step 3: Add the merged identifiers column to the table
    cat["ids"] = Column(data=merged_ids, dtype=object)
    return cat


def assign_type(cat: Table, i: int) -> str:
    """
    Assigns an object type based on two potential type columns.

    Priority is given to 'type_2' unless it is masked or 'None'.

    :param cat: The table containing type columns.
    :type cat: Table
    :param i: Index of the row to process.
    :type i: int
    :return: The assigned type string.
    :rtype: str
    """
    # Check if 'type_2' is masked or equal to 'None', then fall back to
    # 'type_1', otherwise use 'type_2'
    val_2 = cat["type_2"][i]
    if isinstance(val_2, np.ma.core.MaskedConstant) or val_2 == "None":
        cat["type"][i] = cat["type_1"][i]
    else:
        cat["type"][i] = val_2
    return cat["type"][i]


def objectmerging(cat: Table) -> Table:
    """
    Merges the data of each object given in the different providers.

    The object is the same physical one but the data is provided by
    different providers and merged into one entry.

    :param cat: Table containing multiple entries for the same objects.
    :type cat: Table
    :returns: Table with unique object entries.
    :rtype: Table
    """
    cat = idsjoin(cat, "ids_1", "ids_2")
    cat.remove_columns(["ids_1", "ids_2"])

    # Initializing merged type column if it doesn't exist
    if "type" not in cat.colnames:
        cat["type"] = Column(dtype=object, length=len(cat))
        cat["type_1"] = cat["type_1"].astype(object)
        cat["type_2"] = cat["type_2"].astype(object)
        for i in range(len(cat)):
            cat["type"][i] = assign_type(cat, i)
        cat.remove_columns(["type_1", "type_2"])
    return cat


def assign_source_idref(
    cat: Table, sources: Table, paras: list[str], provider: str
) -> Table:
    """
    Joins source identifiers to parameters in a catalog table.

    For each parameter that has a reference column, this function:
    1. Validates and handles existing source ID columns.
    2. Processes null values in reference columns.
    3. Links reference data with source identifiers.
    4. Manages masked parameter values.

    :param cat: Table with empty para_source_id columns.
    :type cat: Table
    :param sources: Table containing reference data.
    :type sources: Table
    :param paras: List of parameters to process.
    :type paras: list[str]
    :param provider: Name of the data provider.
    :type provider: str
    :returns: Catalog table with parameter source IDs added.
    :rtype: Table
    """
    for para in paras:
        ref_column = para + "_ref"
        source_id_col = f"{para}_source_idref"
        value_column = para + "_value"

        # Skip if no reference column exists
        if ref_column not in cat.colnames:
            continue

        # Check for existing source ID column
        if source_id_col in cat.colnames:
            print(
                f"warning, {source_id_col} already in table. "
                "something went wrong with loading"
            )
            cat.remove_column(source_id_col)

        # Replace null values in reference column
        cat = nullvalues(cat, ref_column, "?")

        # Join with sources table to get source IDs for this provider
        source_subset = sources["ref", "source_id"][
            np.where(sources["provider_name"] == provider)
        ]
        cat = join(
            cat,
            source_subset,
            keys_left=ref_column,
            keys_right="ref",
            join_type="left",
        )

        # Rename generic source_id column to parameter-specific name
        cat.rename_column("source_id", source_id_col)
        cat.remove_columns("ref")

        # Handle masked values in parameter column: assign 999999 to ref ID
        if value_column in cat.colnames:
            if isinstance(cat[value_column], column.MaskedColumn):
                for i in cat[value_column].mask.nonzero()[0]:
                    cat[source_id_col][i] = 999999

    return cat


def merge_table(cat1: Table, cat2: Table) -> Table:
    """
    Merges two tables.

    If one table is empty, vstack is used as join would fail.

    :param cat1: First table to merge.
    :type cat1: Table
    :param cat2: Second table to merge.
    :type cat2: Table
    :returns: Merged table.
    :rtype: Table
    """
    if len(cat1) == 0 or len(cat2) == 0:
        return vstack([cat1, cat2])
    return join(cat1, cat2)


def best_para_id(mes_table: Table) -> Table:
    """
    Selects the best identifier for each object based on reference priority.

    :param mes_table: Measurement table for identifiers.
    :type mes_table: Table
    :returns: Table with prioritized identifier rows.
    :rtype: Table
    """
    grouped_mes_table = mes_table.group_by("id_ref")

    # 1. Making simbad identifiers the default best parameters
    simbad_ref = "2000A&AS..143....9W"
    mask = grouped_mes_table.groups.keys["id_ref"] == simbad_ref
    best_para_table = grouped_mes_table.groups[mask]

    # 2. Adding identifiers that are not in the best_para_table yet.
    # Higher quality priority order for provider identifier references.
    # TBD: use id_ref as variable from provider_bibcode
    #        instad of constant""")
    priority_refs = [
        "2022A&A...664A..21Q",
        "2016A&A...595A...1G",
        "priv. comm.",
        "2020A&C....3100370A",
        "2001AJ....122.3466M",
    ]

    for ref in priority_refs:
        mask = grouped_mes_table.groups.keys["id_ref"] == ref
        all_ref_ids = grouped_mes_table.groups[mask]

        # removing those already in best_para_table
        new_ids = all_ref_ids[
            np.where(
                np.invert(np.isin(all_ref_ids["id"], best_para_table["id"]))
            )
        ]
        best_para_table = vstack([best_para_table, new_ids])

    best_para_table.remove_column("id_ref")
    return best_para_table


def best_para_membership(mes_table: Table) -> Table:
    """
    Selects the best membership measurement for parent-child pairs.

    Chooses the row with the maximum membership value.

    :param mes_table: Measurement table for membership.
    :type mes_table: Table
    :returns: Table with best membership entries.
    :rtype: Table
    """
    para = "membership"
    best_para_table = mes_table[:0].copy()
    grouped_mes_table = mes_table.group_by(
        ["child_object_idref", "parent_object_idref"]
    )
    indices = grouped_mes_table.groups.indices

    for i in range(len(indices) - 1):
        idx_start, idx_end = indices[i], indices[i + 1]
        group = grouped_mes_table[idx_start:idx_end]

        if len(group) == 1:
            best_para_table.add_row(group[0])
        else:
            with_value = group[np.where(group[para] != 999999)]
            if len(with_value) > 0:
                max_val = max(with_value[para])
                # Add first row matching the max value
                for row in group:
                    if row[para] == max_val:
                        best_para_table.add_row(row)
                        break  # make sure not multiple of same max
                        # value are added
            else:
                # if none of the objects has a membership entry
                # then pick just first one
                best_para_table.add_row(group[0])

    return best_para_table


def _get_parameter_columns(para: str) -> list[str]:
    """
    Helper function to determine which columns to include based on parameter.

    :param para: Parameter name.
    :type para: str
    :returns: List of column names to include.
    :rtype: list[str]
    """
    cols = ["main_id"]
    if para == "binary":
        cols.extend([f"{para}_flag", f"{para}_qual", f"{para}_source_idref"])
    elif para == "mass_pl":
        cols.extend(
            [
                f"{para}_value",
                f"{para}_rel",
                f"{para}_err_max",
                f"{para}_err_min",
                f"{para}_qual",
                f"{para}_sini_flag",
                f"{para}_source_idref",
            ]
        )
    elif para == "sep_ang":
        cols.extend(
            [
                f"{para}_value",
                f"{para}_err",
                f"{para}_obs_date",
                f"{para}_qual",
                f"{para}_source_idref",
            ]
        )
    else:
        cols.extend(
            [
                f"{para}_value",
                f"{para}_err",
                f"{para}_qual",
                f"{para}_source_idref",
            ]
        )
    return cols


def _find_best_quality_measurement(group: Table, para: str) -> Row | None:
    """
    Helper function to find the highest quality measurement in a group.

    Quality levels: A > B > C > D > E > ?.

    :param group: Group of measurements for a single object.
    :type group: Table
    :param para: Parameter name.
    :type para: str
    :returns: Row with highest quality measurement or None.
    :rtype: Row or None
    """
    quality_levels = ["A", "B", "C", "D", "E", "?"]
    qual_column = f"{para}_qual"

    for quality in quality_levels:
        for row in group:
            if row[qual_column] == quality:
                return row
    return None


def best_para(para: str, mes_table: Table) -> Table:
    """
    Selects the highest quality measurement for each object in the table.

    :param para: Parameter name (e.g., 'mass', 'id').
    :type para: str
    :param mes_table: Table containing measurements.
    :type mes_table: Table
    :returns: Table with highest quality rows for each unique object.
    :rtype: Table
    """
    # Special case handlers
    if para == "id":
        return best_para_id(mes_table)
    if para == "membership":
        return best_para_membership(mes_table)

    # Define columns based on parameter type
    columns = _get_parameter_columns(para)

    # Select only needed columns and create empty result table
    mes_table = mes_table[columns]
    best_para_table = mes_table[:0].copy()

    # Group by main_id and process each group
    grouped_mes_table = mes_table.group_by("main_id")

    for group in grouped_mes_table.groups:
        best_measurement = _find_best_quality_measurement(group, para)
        if best_measurement is not None:
            best_para_table.add_row(best_measurement)

    return best_para_table


def best_parameters_ingestion(
    cat_mes: Table,
    cat_basic: Table,
    para: str,
    columns: list[str] | None = None,
) -> Table:
    """
    Updates a basic table with the best measurements from a measurement table.

    :param cat_mes: Table containing multiple measurements.
    :type cat_mes: Table
    :param cat_basic: Table to be updated.
    :type cat_basic: Table
    :param para: Parameter name.
    :type para: str
    :param columns: Columns to remove from cat_basic before joining.
    :type columns: list[str] or None
    :returns: Updated basic table.
    :rtype: Table
    """
    best_para_cat_mes = best_para(para, cat_mes)
    if columns:
        cat_basic.remove_columns(columns)
    return join(cat_basic, best_para_cat_mes, join_type="left")


def provider_data_merging(
    cat: dict[str, Table],
    table_name: str,
    prov_tables_dict: dict[str, dict[str, Table]],
    o_merging: bool = False,
    para_match: bool = False,
) -> dict[str, Table]:
    """
    Merges the data from the different providers for a specific table.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :param table_name: Name of the table to build/merge.
    :type table_name: str
    :param prov_tables_dict: Dictionary mapping providers to their tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :param o_merging: Whether to perform object-level merging (ids and types).
    :type o_merging: bool
    :param para_match: Whether to perform source ID matching.
    :type para_match: bool
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    print(f"Building {table_name} table ...")
    for prov_name, prov_data in prov_tables_dict.items():
        if para_match:
            cat = matching_parameters(
                cat, prov_name, prov_tables_dict, table_name
            )

        if table_name not in cat or len(cat[table_name]) == 0:
            cat[table_name] = prov_data[table_name]
        else:
            cat = join_different_provider_data(
                cat, o_merging, prov_name, prov_tables_dict, table_name
            )
    return cat


def join_different_provider_data(
    cat: dict[str, Table],
    o_merging: bool,
    prov_name: str,
    prov_tables_dict: dict[str, dict[str, Table]],
    table_name: str,
) -> dict[str, Table]:
    """
    Joins data from a specific provider into the main table.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :param o_merging: Whether to perform object-level merging.
    :type o_merging: bool
    :param prov_name: Name of the provider.
    :type prov_name: str
    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :param table_name: Name of the table to merge.
    :type table_name: str
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    prov_table = prov_tables_dict[prov_name][table_name]
    if len(prov_table) > 0:
        if o_merging:
            cat[table_name] = join(
                cat[table_name],
                prov_table,
                keys="main_id",
                join_type="outer",
            )
            cat[table_name] = objectmerging(cat[table_name])
        else:
            cat[table_name] = join(
                cat[table_name],
                prov_table,
                join_type="outer",
            )
    return cat


def matching_parameters(
    cat: dict[str, Table],
    prov_name: str,
    prov_tables_dict: dict[str, dict[str, Table]],
    table_name: str,
) -> dict[str, Table]:
    """
    Redefines source reference columns with their corresponding IDs.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :param prov_name: Name of the provider.
    :type prov_name: str
    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :param table_name: Name of the table to process.
    :type table_name: str
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    prov_table = prov_tables_dict[prov_name][table_name]
    if len(prov_table) > 0:
        paras = paras_dict.copy()
        provider_name = prov_tables_dict[prov_name]["provider"][
            "provider_name"
        ][0]
        prov_tables_dict[prov_name][table_name] = assign_source_idref(
            prov_table, cat["sources"], paras[table_name], provider_name
        )
    return cat


def unify_null_values(cat: dict[str, Table]) -> dict[str, Table]:
    """
    Unifies null values ('N', 'N/A') to '?' across specific tables/columns.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    print("Unifying null values...")
    table_keys = [
        "star_basic",
        "planet_basic",
        "disk_basic",
        "mes_mass_pl",
        "mes_teff_st",
        "mes_radius_st",
        "mes_mass_st",
        "mes_binary",
    ]
    columns_map = [
        [
            "coo_qual",
            "coo_gal_qual",
            "plx_qual",
            "dist_st_qual",
            "sep_ang_qual",
            "teff_st_qual",
            "radius_st_qual",
            "binary_flag",
            "binary_qual",
            "mass_st_qual",
            "sptype_qual",
            "class_temp",
            "class_temp_nr",
        ],
        ["mass_pl_qual"],
        ["rad_qual", "rad_rel"],
        ["mass_pl_qual"],
        ["teff_st_qual"],
        ["radius_st_qual"],
        ["mass_st_qual"],
        ["binary_qual"],
    ]
    for key, cols in zip(table_keys, columns_map):
        if key in cat:
            for col in cols:
                cat[key] = replace_value(cat[key], col, "N", "?")
                cat[key] = replace_value(cat[key], col, "N/A", "?")
    return cat


def build_sources_table(
    prov_tables_dict: dict[str, dict[str, Table]],
) -> dict[str, Table]:
    """
    Initializes the catalog and builds the unique sources table.

    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :returns: Dictionary containing the initialized sources table.
    :rtype: dict[str, Table]
    """
    cat = empty_dict.copy()
    empty = empty_dict_wit_columns.copy()

    # for the sources and objects joins tables from different prov_tables_dict
    cat = provider_data_merging(cat, "sources", prov_tables_dict)

    # Adding empty template columns and keeping unique entries
    cat["sources"] = vstack([cat["sources"], empty["sources"]])
    cat["sources"] = unique(cat["sources"], silent=True)
    cat["sources"]["source_id"] = [j + 1 for j in range(len(cat["sources"]))]

    return cat


def build_objects_table(
    cat: dict[str, Table], prov_tables_dict: dict[str, dict[str, Table]]
) -> dict[str, Table]:
    """
    Builds the objects table and assigns unique object IDs.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    cat = provider_data_merging(
        cat, "objects", prov_tables_dict, o_merging=True
    )
    cat["objects"]["object_id"] = [j + 1 for j in range(len(cat["objects"]))]

    # At one point I would like to be able to merge objects with main_id
    # NAME Proxima Centauri b and Proxima Centauri b
    return cat


def build_provider_table(
    cat: dict[str, Table], prov_tables_dict: dict[str, dict[str, Table]]
) -> dict[str, Table]:
    """
    Builds the provider table.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    empty = empty_dict_wit_columns.copy()
    cat = provider_data_merging(cat, "provider", prov_tables_dict)

    # I do this to get those columns that are empty in the data
    cat["provider"] = vstack([cat["provider"], empty["provider"]])
    return cat


def _handle_object_id_linking(table: Table, objects: Table) -> Table:
    """
    Helper to join object_id from objects table into another table.

    :param table: Table needing object_idref.
    :type table: Table
    :param objects: Main objects table.
    :type objects: Table
    :returns: Table with object_idref added.
    :rtype: Table
    """
    if "object_idref" in table.colnames:
        table.remove_column("object_idref")
    table = join(table, objects["object_id", "main_id"], join_type="left")
    table.rename_column("object_id", "object_idref")
    return table


def _process_h_link(cat: dict[str, Table]) -> dict[str, Table]:
    """
    Processes h_link table to link parent/child IDs and find best memberships.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    h_link = cat["h_link"]
    objects = cat["objects"]

    # Link child_object_idref
    if "child_object_idref" in h_link.colnames:
        h_link.remove_column("child_object_idref")
    h_link = join(
        h_link,
        objects["object_id", "main_id"],
        keys="main_id",
        join_type="left",
    )
    h_link.rename_columns(
        ["object_id", "main_id"], ["child_object_idref", "child_main_id"]
    )

    # Link parent_object_idref
    if "parent_object_idref" in h_link.colnames:
        h_link.remove_column("parent_object_idref")
    # Only keep links where parent is also in our objects table
    h_link = join(
        h_link,
        objects["object_id", "main_id"],
        keys_left="parent_main_id",
        keys_right="main_id",
    )
    h_link.remove_column("main_id")
    h_link.rename_column("object_id", "parent_object_idref")

    cat["h_link"] = h_link
    cat["best_h_link"] = best_para("membership", h_link)
    return cat


def _process_basic_tables(cat: dict[str, Table]) -> dict[str, Table]:
    """
    Specialized processing for star_basic and planet_basic.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    objects = cat["objects"]

    # Star basic: include all objects typed as stars or systems
    stars = objects["object_id", "main_id"][np.where(objects["type"] == "st")]
    systems = objects["object_id", "main_id"][np.where(objects["type"] == "sy")]
    temp = vstack([stars, systems])
    temp.rename_column("object_id", "object_idref")

    cat["star_basic"] = join(
        cat["star_basic"],
        temp,
        join_type="outer",
        keys=["object_idref", "main_id"],
    )

    # Planet basic: include all objects typed as planets
    planets = objects["object_id", "main_id"][np.where(objects["type"] == "pl")]
    planets.rename_column("object_id", "object_idref")
    cat["planet_basic"] = planets
    return cat


def build_rest_of_tables(
    cat: dict[str, Table], prov_tables_dict: dict[str, dict[str, Table]]
) -> dict[str, Table]:
    """
    Builds all remaining tables (basic, measurements, etc.) and performs links.

    :param cat: Dictionary of cumulative tables.
    :type cat: dict[str, Table]
    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :returns: Updated dictionary of cumulative tables.
    :rtype: dict[str, Table]
    """
    empty_dict_cols = empty_dict_wit_columns.copy()

    # Skip first 3 tables (sources, objects, provider)
    for table_name in islice(cat, 3, None):
        cat = provider_data_merging(
            cat, table_name, prov_tables_dict, para_match=True
        )

        cat[table_name] = vstack([cat[table_name], empty_dict_cols[table_name]])
        cat[table_name] = cat[table_name].filled()

        # Link object_idrefs if applicable
        if (
            "object_idref" in cat[table_name].colnames
            and len(cat[table_name]) > 0
        ):
            cat[table_name] = _handle_object_id_linking(
                cat[table_name], cat["objects"]
            )

        # Specialized handling by table name
        if table_name == "ident":
            cat[table_name] = best_para("id", cat[table_name])
        elif table_name == "h_link":
            cat = _process_h_link(cat)
        elif table_name == "planet_basic":
            cat = _process_basic_tables(cat)
        elif table_name == "mes_teff_st":
            cat["star_basic"] = best_parameters_ingestion(
                cat[table_name],
                cat["star_basic"],
                "teff_st",
                [
                    "teff_st_value",
                    "teff_st_err",
                    "teff_st_qual",
                    "teff_st_source_idref",
                    "teff_st_ref",
                ],
            )
        elif table_name == "mes_radius_st":
            cat["star_basic"] = best_parameters_ingestion(
                cat[table_name],
                cat["star_basic"],
                "radius_st",
                [
                    "radius_st_value",
                    "radius_st_err",
                    "radius_st_qual",
                    "radius_st_source_idref",
                    "radius_st_ref",
                ],
            )
        elif table_name == "mes_mass_st":
            cat["star_basic"] = best_parameters_ingestion(
                cat[table_name],
                cat["star_basic"],
                "mass_st",
                [
                    "mass_st_value",
                    "mass_st_err",
                    "mass_st_qual",
                    "mass_st_source_idref",
                    "mass_st_ref",
                ],
            )
        elif table_name == "mes_mass_pl":
            cat["planet_basic"] = best_parameters_ingestion(
                cat[table_name], cat["planet_basic"], "mass_pl"
            )
        elif table_name == "mes_binary":
            cat["star_basic"] = best_parameters_ingestion(
                cat[table_name],
                cat["star_basic"],
                "binary",
                [
                    "binary_flag",
                    "binary_qual",
                    "binary_source_idref",
                    "binary_ref",
                ],
            )
        elif table_name == "mes_sep_ang":
            cat["star_basic"] = best_parameters_ingestion(
                cat[table_name],
                cat["star_basic"],
                "sep_ang",
                [
                    "sep_ang_value",
                    "sep_ang_err",
                    "sep_ang_obs_date",
                    "sep_ang_qual",
                    "sep_ang_source_idref",
                    "sep_ang_ref",
                ],
            )

        cat[table_name] = cat[table_name].filled()
        if len(cat[table_name]) == 0:
            print(f"warning: empty table {table_name}")
        else:
            cat[table_name] = unique(cat[table_name], silent=True)

    return cat


def build_tables(
    prov_tables_dict: dict[str, dict[str, Table]],
) -> dict[str, Table]:
    """
    Orchestrates the building of all database tables.

    :param prov_tables_dict: Dictionary of provider tables.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :returns: Dictionary of built tables.
    :rtype: dict[str, Table]
    """
    cat = build_sources_table(prov_tables_dict)
    cat = build_objects_table(cat, prov_tables_dict)
    cat = build_provider_table(cat, prov_tables_dict)
    cat = build_rest_of_tables(cat, prov_tables_dict)
    return cat


def building(prov_tables_dict: dict[str, dict[str, Table]]) -> dict[str, Table]:
    """
    Builds the complete LIFE database from provider tables and saves it.

    :param prov_tables_dict: Dictionary containing data from providers
        Simbad, Grant Kennedy, Exo-MerCat, Gaia and WDS.
    :type prov_tables_dict: dict[str, dict[str, Table]]
    :returns: Dictionary of processed tables.
    :rtype: dict[str, Table]
    """
    cat = build_tables(prov_tables_dict)

    # Ensure star_basic has no masked entries after multi-measurement ingestions
    cat["star_basic"] = cat["star_basic"].filled()

    cat = unify_null_values(cat)

    # TBD: Add exact object distance cut. So far for correct treatment
    #       of boundary objects 10% additional distance cut used""")

    print("Saving data...")
    save(
        list(cat.values()),
        list(cat.keys()),
        location=Path().data,
    )
    return cat
