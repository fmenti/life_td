"""
Generates the data for the database for each of the data providers separately.
"""

from collections.abc import Sequence

import numpy as np  # arrays
from astropy.table import MaskedColumn, Table, join, setdiff, unique, vstack
from provider.assign_quality_funcs import assign_quality
from provider.utils import (
    OidCreator,
    create_provider_table,
    create_sources_table,
    fetch_main_id,
    nullvalues,
    query,
    replace_value,
)
from sdata import empty_dict
from utils.io import save

# ---------------define queries--------------------------------------
adql_queries_moduls = {
    "select_statement": """SELECT b.main_id,b.ra AS coo_ra,b.dec AS coo_dec,
        b.coo_err_angle, b.coo_err_maj, b.coo_err_min,b.oid,
        b.coo_bibcode AS coo_ref, b.coo_qual,b.sp_type AS sptype_string,
        b.sp_qual AS sptype_qual, b.sp_bibcode AS sptype_ref,
        b.plx_err, b.plx_value, b.plx_bibcode AS plx_ref,b.plx_qual,
        h_link.membership, h_link.parent AS parent_oid,
        h_link.link_bibcode AS h_link_ref, a.otypes,ids.ids,
        f.I as mag_i_value, f.J as mag_j_value, f.K as mag_k_value
        """,
    "tables_statement": """
    FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    LEFT JOIN allfluxes AS f ON b.oid=f.oidref

    """,
}


def main_adql_queries(plx_cut: float) -> str:
    """
    Build the main ADQL query for SIMBAD with a parallax cut.

    :param plx_cut: Minimum parallax in milliarcseconds to include.
    :type plx_cut: float
    :returns: Complete ADQL string for the TAP service.
    :rtype: str
    """
    return (
        adql_queries_moduls["select_statement"]
        + adql_queries_moduls["tables_statement"]
        + "WHERE b.plx_value >="
        + str(plx_cut)
    )


adql_upload_queries = {
    "sy_without_plx_but_child_with_upload": adql_queries_moduls[
        "select_statement"
    ]
    + adql_queries_moduls["tables_statement"]
    + """JOIN TAP_UPLOAD.t1 ON b.oid=t1.parent_oid
            WHERE (b.plx_value IS NULL) AND (otype='**..')""",
    "pl_without_plx_but_host_with_upload": adql_queries_moduls[
        "select_statement"
    ]
    + adql_queries_moduls["tables_statement"]
    + """JOIN TAP_UPLOAD.t1 ON b.oid=t1.oid
            WHERE (b.plx_value IS NULL) AND (otype='Pl..')""",
    "ids_from_upload": """SELECT id, t1.*
                          FROM ident
                          JOIN TAP_UPLOAD.t1 ON oidref = t1.oid""",
}


def create_simbad_helpertable(
    distance_cut_in_pc: float,
    test_objects: Sequence[str] | None,
) -> tuple[Table, dict[str, Table]]:
    """
    Create the main SIMBAD helper table.

    :param distance_cut_in_pc: Distance up to which stars are included (pc).
    :type distance_cut_in_pc: float
    :param test_objects: Optional list of object names to trace through the
        filtering pipeline for debugging.
    :type test_objects: list[str] or None
    :returns: Helper table and a dictionary of database tables.
    :rtype: (astropy.table.Table, dict[str, astropy.table.Table])
    """
    plx_in_mas_cut = 1000.0 / distance_cut_in_pc
    # making cut a bit bigger for correct treatment of objects on boundary
    plx_cut = plx_in_mas_cut - plx_in_mas_cut / 10.0

    sim = empty_dict.copy()
    sim["provider"] = create_provider_table(
        "SIMBAD",
        "http://simbad.u-strasbg.fr:80/simbad/sim-tap",
        "2000A&AS..143....9W",
    )

    # ------------------querrying----------------------------------------
    # perform query for objects with in distance given
    sim_helptab = query(
        sim["provider"]["provider_url"][0], main_adql_queries(plx_cut)
    )
    save([sim_helptab], ["sim_helptab_query"])
    # querries parent and children objects with no parallax value
    parents_without_plx = query(
        sim["provider"]["provider_url"][0],
        adql_upload_queries["sy_without_plx_but_child_with_upload"],
        [sim_helptab],
    )
    save([parents_without_plx], ["parents_without_plx_query"])
    children_without_plx = query(
        sim["provider"]["provider_url"][0],
        adql_upload_queries["pl_without_plx_but_host_with_upload"],
        [sim_helptab],
    )
    save([children_without_plx], ["children_without_plx_query"])

    test_objects = np.array(test_objects)
    if len(test_objects) > 0:
        print(
            "in sim through plx query",
            test_objects[
                np.where(np.isin(test_objects, sim_helptab["main_id"]))
            ],
        )
        print(
            "in sim through child plx query",
            test_objects[
                np.where(np.isin(test_objects, parents_without_plx["main_id"]))
            ],
        )
        print(
            "in sim through parent plx query",
            test_objects[
                np.where(np.isin(test_objects, children_without_plx["main_id"]))
            ],
        )

    # adding of no_parallax objects to rest of simbad query objects

    sim_helptab = vstack([sim_helptab, parents_without_plx])
    sim_helptab = vstack([sim_helptab, children_without_plx])

    print(" sorting object types...")

    # sorting from object type into star, system and planet type
    sim_helptab["type"] = ["None" for i in range(len(sim_helptab))]
    sim_helptab["type"] = sim_helptab["type"].astype(object)
    sim_helptab["binary_flag"] = np.array(
        ["False" for i in range(len(sim_helptab))], dtype=object
    )
    to_remove_list = []
    removed_otypes = []
    for i in range(len(sim_helptab)):
        # planets
        if "Pl" in sim_helptab["otypes"][i]:
            sim_helptab["type"][i] = "pl"
        # stars
        elif "*" in sim_helptab["otypes"][i]:
            # system containing multiple stars
            if "**" in sim_helptab["otypes"][i]:
                sim_helptab["type"][i] = "sy"
                sim_helptab["binary_flag"][i] = "True"
            # individual stars
            else:
                sim_helptab["type"][i] = "st"
        else:
            removed_otypes.append(sim_helptab["otypes"][i])
            # most likely single brown dwarfs
            # -> wait, they have BD* and get caught above...
            # so I am not as restrictive as I thought I am here.
            # shouldn't matter as they shouldn't make it into the catalog
            # because of the spectral type criteria
            # storing information for later removal from table called simbad
            to_remove_list.append(i)
    # removing any objects that are neither planet, star nor system in type
    if to_remove_list != []:
        print(
            "removing",
            len(removed_otypes),
            " objects that had object types:",
            list(set(removed_otypes)),
        )
        print(
            "example object of them:", sim_helptab["main_id"][to_remove_list[0]]
        )
        sim_helptab.remove_rows(to_remove_list)

    if len(test_objects) > 0:
        print(
            "in sim through otype criteria",
            test_objects[
                np.where(np.isin(test_objects, sim_helptab["main_id"]))
            ],
        )

    return sim_helptab, sim


def stars_in_multiple_system(
    cat: Table,
    sim_h_link: Table,
    all_objects: Table,
) -> Table:
    """
    Assign object type "st" to a subset of systems.

    Assign the object type "st" to those entries in a system table that are
    members of multiple systems but do not have stellar children.

    :param cat: Table alias containing objects of type "sy".
    :type cat: astropy.table.Table
    :param sim_h_link: Table with columns parent_main_id, main_id and
        h_link_ref (parent-child pairs).
    :type sim_h_link: astropy.table.Table
    :param all_objects: Table copy containing columns main_id and type.
        Rows are all objects with child objects and their children.
    :type all_objects: astropy.table.Table
    :returns: The input table with selected rows having their type set to "st".
    :rtype: astropy.table.Table
    """

    # Initiate parent table with consistent names.
    def split_into_sy_w_and_wo_child(
        sim_h_link: Table,
        cat: Table,
        all_objects: Table,
    ) -> tuple[Table, Table]:
        parents = sim_h_link["parent_main_id", "main_id", "h_link_ref"][:]
        parents.rename_column("main_id", "child_main_id")
        parents.rename_column("parent_main_id", "main_id")

        # Objects of type system without children.
        sy_wo_child = setdiff(
            cat["main_id", "type", "sptype_string"][:],
            parents[:],
            keys=["main_id"],
        )

        # Objects of type system that have children.
        sy_w_child = join(
            parents[:],
            cat["main_id", "type", "sptype_string"][:],
            keys=["main_id"],
        )

        # Add information about those children.
        all_objects.rename_columns(
            ["type", "main_id"], ["child_type", "child_main_id"]
        )
        sy_w_child = join(
            sy_w_child[:],
            all_objects["child_type", "child_main_id"][:],
            keys=["child_main_id"],
            join_type="left",
        )
        return sy_wo_child, sy_w_child

    sy_wo_child, sy_w_child = split_into_sy_w_and_wo_child(
        sim_h_link, cat, all_objects
    )

    # Get those that have planetary children or no children at all.
    def get_sy_wo_child_st(
        sy_wo_child: Table,
        sy_w_child: Table,
    ) -> Table:
        sy_w_child_pl = sy_w_child[np.where(sy_w_child["child_type"] == "pl")]
        if len(sy_w_child_pl) == 0:
            return sy_wo_child

        # Join with list of systems without children.
        sy_wo_child_st = vstack([sy_wo_child[:], sy_w_child_pl[:]])
        sy_wo_child_st.remove_column("child_type")
        return sy_wo_child_st

    sy_wo_child_st = get_sy_wo_child_st(sy_wo_child, sy_w_child)

    def special_treatment_spectral_type(sy_wo_child_st: Table) -> Table:
        # No "+" in sptype_string because that indicates binarity.
        temp = [len(i.split("+")) == 1 for i in sy_wo_child_st["sptype_string"]]
        mask = np.array(temp)
        indices = list(np.where(mask)[0])
        return sy_wo_child_st[:][indices]

    single_sptype = special_treatment_spectral_type(sy_wo_child_st)

    # Reassign type to "st".
    cat["type"][np.where(np.isin(cat["main_id"], single_sptype["main_id"]))] = [
        "st"
        for _ in range(
            len(
                cat[np.where(np.isin(cat["main_id"], single_sptype["main_id"]))]
            )
        )
    ]
    return cat


def creating_helpertable_stars(
    sim_helptab: Table,
) -> Table:
    """
    Create a helper table for stars.

    Removes planets, de-duplicates by main_id, and returns the table that
    contains only systems and stars.

    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict[str, astropy.table.Table]
    :returns: Helper table for star/system entries.
    :rtype: astropy.table.Table
    """
    temp_stars = sim_helptab[np.where(sim_helptab["type"] != "pl")]
    # Removing double objects (in there due to multiple parents).
    stars = Table(unique(temp_stars, keys="main_id"), copy=True)
    return stars


def expanding_helpertable_stars(
    sim_helptab: Table,
    sim: dict[str, Table],
    stars: Table,
) -> Table:
    """
    Expand the star helper table with additional fields and quality flags.

    Updates the multiplicity type for certain systems, sets binary flags,
    normalizes reference fields, and assigns quality flags.

    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict[str, astropy.table.Table]
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.Table
    :returns: Expanded star helper table.
    :rtype: astropy.table.Table
    """
    # Update multiplicity object type: no children and sptype without "+"
    # implies the type needs to be "st".
    stars[np.where(stars["type"] == "sy")] = stars_in_multiple_system(
        stars[np.where(stars["type"] == "sy")],
        sim["h_link"][:],
        sim_helptab["main_id", "type"][:],
    )

    # Binary_flag "True" for all stars with parents. That is, stars[main_id]
    # in sim_h_link[child_main_id].
    stars["binary_flag"][
        np.where(np.isin(stars["main_id"], sim["h_link"]["main_id"]))
    ] = [
        "True"
        for _ in range(
            len(
                stars[
                    np.where(
                        np.isin(stars["main_id"], sim["h_link"]["main_id"])
                    )
                ]
            )
        )
    ]

    # Change null value of plx_qual.
    stars["plx_qual"] = stars["plx_qual"].astype(object)
    stars = replace_value(stars, "plx_qual", "", stars["plx_qual"].fill_value)

    # Initialize and fill photometry references if values exist.
    for band in ["i", "j", "k"]:
        ref_col = f"mag_{band}_ref"
        val_col = f"mag_{band}_value"

        stars[ref_col] = MaskedColumn(
            dtype=object, length=len(stars), mask=[True] * len(stars)
        )
        mask_has_val = np.where(stars[val_col].mask == False)

        stars[ref_col][mask_has_val] = [
            sim["provider"]["provider_bibcode"][0]
            for _ in range(len(stars[ref_col][mask_has_val]))
        ]

    # Spectral type reference: fill empty with provider bibcode where value set.
    non_empty_sptype = np.where(stars["sptype_string"] != "")
    stars[non_empty_sptype] = replace_value(
        stars[non_empty_sptype],
        "sptype_ref",
        "",
        sim["provider"]["provider_bibcode"][0],
    )

    # Make the remaining ones into masked entries.
    stars["sptype_ref"] = MaskedColumn(stars["sptype_ref"])
    stars["sptype_ref"].mask[np.where(stars["sptype_ref"] == "")] = [
        True for _ in range(len(stars[np.where(stars["sptype_ref"] == "")]))
    ]

    # Handle parallax and coordinate references similarly.
    for colname, colval in zip(
        ["plx_ref", "coo_ref"],
        ["plx_value", "coo_ra"],
    ):
        stars[colname] = MaskedColumn(stars[colname])
        mask_has_val = np.where(stars[colval].mask == False)
        stars[mask_has_val] = replace_value(
            stars[mask_has_val],
            colname,
            "",
            sim["provider"]["provider_bibcode"][0],
        )
        # Make the remaining ones into masked entries.
        stars[colname].mask[np.where(stars[colname] == "")] = [
            True for _ in range(len(stars[np.where(stars[colname] == "")]))
        ]

    # Binarity columns.
    stars["binary_ref"] = [
        sim["provider"]["provider_bibcode"][0] for _ in range(len(stars))
    ]
    stars = assign_quality(stars, "binary_qual", special_mode="sim_binary")
    return stars


def create_ident_table(sim_helptab: Table, sim: dict[str, Table]) -> Table:
    """
    Create the identifier table.

    Upload the object IDs to SIMBAD to fetch the corresponding identifiers
    (id) and add a reference column.

    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict[str, astropy.table.Table]
    :returns: Identifier table with id and id_ref columns.
    :rtype: astropy.table.Table
    """
    sim_ident = query(
        sim["provider"]["provider_url"][0],
        adql_upload_queries["ids_from_upload"],
        [sim_helptab["oid", "main_id"][:].copy()],
    )
    sim_ident["id_ref"] = [
        sim["provider"]["provider_bibcode"][0] for _ in range(len(sim_ident))
    ]
    sim_ident["id_ref"] = sim_ident["id_ref"].astype(object)
    sim_ident.remove_column("oid")
    return sim_ident


def create_h_link_table(
    sim_helptab: Table,
    sim: dict[str, Table],
    stars: Table,
) -> Table:
    """
    Create the hierarchical link (parent-child) table.

    Filters to hierarchical multiples (i.e., removes cluster/association
    parents), resolves parent oids to main_ids, normalizes null values,
    and adds reference information.

    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict[str, astropy.table.Table]
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.Table
    :returns: Hierarchical link table with parent_main_id, main_id,
        membership, and h_link_ref.
    :rtype: astropy.table.Table
    """
    sim_h_link = sim_helptab[
        "main_id", "parent_oid", "h_link_ref", "membership"
    ]

    # Keep only links whose parents are stellar/system objects (not clusters).
    sim_h_link = sim_h_link[
        np.where(np.isin(sim_h_link["parent_oid"], stars["oid"]))
    ]

    # Resolve parent oid to main_id.
    sim_h_link = fetch_main_id(
        sim_h_link, OidCreator(name="parent_main_id", colname="parent_oid")
    )
    sim_h_link.remove_column("parent_oid")

    # Convert to int and normalize nulls to a sentinel value for joins.
    sim_h_link["membership"] = sim_h_link["membership"].astype(int)
    sim_h_link = nullvalues(sim_h_link, "membership", 999999)

    # Ensure reference is set.
    sim_h_link = replace_value(
        sim_h_link, "h_link_ref", "", sim["provider"]["provider_bibcode"][0]
    )
    sim_h_link = unique(sim_h_link)
    return sim_h_link


def create_objects_table(sim_helptab: Table, stars: Table) -> Table:
    """
    Create the objects table.

    Combines unique planet rows with star/system rows, keeping main_id,
    ids, and type.

    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.Table
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.Table
    :returns: Objects table with main_id, ids, and type.
    :rtype: astropy.table.Table
    """
    # Planets subset, unique by main_id.
    temp_sim_planets = sim_helptab["main_id", "ids", "type"][
        np.where(sim_helptab["type"] == "pl")
    ]
    sim_planets = Table(unique(temp_sim_planets, keys="main_id"), copy=True)

    # Stack planets and stars/systems.
    sim_objects = vstack(
        [sim_planets["main_id", "ids", "type"], stars["main_id", "ids", "type"]]
    )
    sim_objects["ids"] = sim_objects["ids"].astype(object)
    # tbd: add identifier simbad main_id without leading * and whitespaces
    return sim_objects


def create_sim_sources_table(stars: Table, sim: dict[str, Table]) -> Table:
    """
    Create the sources table.

    Collects and normalizes reference columns from provider, stars,
    h_link, and ident tables into a unified sources table.

    :param sim: Dictionary of database table names and tables.
    :type sim: dict[str, astropy.table.Table]
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.Table
    :returns: Sources table with deduplicated references.
    :rtype: astropy.table.Table
    """
    tables = [sim["provider"], stars, sim["h_link"], sim["ident"]]
    # Columns that contain reference information per input table.
    ref_columns = [
        ["provider_bibcode"],
        [
            "coo_ref",
            "plx_ref",
            "mag_i_ref",
            "mag_j_ref",
            "mag_k_ref",
            "binary_ref",
            "sptype_ref",
        ],
        ["h_link_ref"],
        ["id_ref"],
    ]
    sim_sources = create_sources_table(
        tables, ref_columns, sim["provider"]["provider_name"][0]
    )
    return sim_sources


def create_star_basic_table(stars: Table) -> Table:
    """
    Create the table with basic stellar data.

    Converts certain columns to strings to support downstream joins that
    rely on fixed-width string dtype compatibility.

    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.Table
    :returns: Basic stellar data table.
    :rtype: astropy.table.Table
    """
    sim_star_basic = stars[
        "main_id",
        "coo_ra",
        "coo_dec",
        "coo_err_angle",
        "coo_err_maj",
        "coo_err_min",
        "coo_qual",
        "coo_ref",
        "mag_i_value",
        "mag_i_ref",
        "mag_j_value",
        "mag_j_ref",
        "mag_k_value",
        "mag_k_ref",
        "sptype_string",
        "sptype_qual",
        "sptype_ref",
        "plx_value",
        "plx_err",
        "plx_qual",
        "plx_ref",
    ]
    # Change type from object to string for later join functions.
    sim_star_basic["sptype_string"] = sim_star_basic["sptype_string"].astype(
        str
    )
    sim_star_basic["sptype_qual"] = sim_star_basic["sptype_qual"].astype(str)
    sim_star_basic["sptype_ref"] = sim_star_basic["sptype_ref"].astype(str)
    return sim_star_basic


def provider_simbad(
    distance_cut_in_pc: float,
    test_objects: Sequence[str] | None = None,
) -> dict[str, Table]:
    """
    Obtain and arrange SIMBAD data.

    Runs the SIMBAD queries, builds helper tables, expands star data with
    quality flags and references, and assembles the final set of SIMBAD
    provider tables.

    :param distance_cut_in_pc: Distance up to which stars are included (pc).
    :type distance_cut_in_pc: float
    :param test_objects: Optional list of object names to trace through the
        filtering pipeline for debugging.
    :type test_objects: list[str] or None
    :returns: Dictionary of SIMBAD tables (sources, provider, objects, ident,
        h_link, star_basic, mes_binary).
    :rtype: dict[str, astropy.table.Table]
    """
    # Normalize optional test_objects to a list for downstream functions.
    _test_objects = list(test_objects) if test_objects is not None else []

    sim_helptab, sim = create_simbad_helpertable(
        distance_cut_in_pc, _test_objects
    )
    stars = creating_helpertable_stars(sim_helptab)
    sim["ident"] = create_ident_table(sim_helptab, sim)
    sim["h_link"] = create_h_link_table(sim_helptab, sim, stars)
    stars = expanding_helpertable_stars(sim_helptab, sim, stars)
    sim["objects"] = create_objects_table(sim_helptab, stars)
    sim["sources"] = create_sim_sources_table(stars, sim)
    sim["star_basic"] = create_star_basic_table(stars)
    sim["mes_binary"] = stars[
        "main_id", "binary_flag", "binary_qual", "binary_ref"
    ]

    save(list(sim.values()), ["sim_" + element for element in list(sim.keys())])
    return sim
