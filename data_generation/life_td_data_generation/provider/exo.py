"""
Generates the data for the database for the provider Exo-Mercat.
Distance cut is performed by joining on life_db simbad objects.
"""

import numpy as np  # arrays
from astropy.io import ascii
from astropy.table import Column, Table, join, setdiff, unique, vstack
from provider.assign_quality_funcs import assign_quality
from provider.utils import (
    IdentifierCreator,
    create_provider_table,
    create_sources_table,
    distance_cut,
    fetch_main_id,
    ids_from_ident,
    lower_quality,
    nullvalues,
    query,
    replace_value,
)
from sdata import empty_dict
from utils.io import Path, load, save, stringtoobject


def query_exomercat(adql_query: str, exo: dict) -> Table:
    """
    Query the Exo-MerCat TAP service for the helper table.

    :param str adql_query: ADQL query string executed on the Exo-MerCat
        TAP service.
    :param exo: Mutable dictionary that will receive the provider table.
        Expected to be the shared provider dict used through this module.
    :type exo: dict[str, Table]
    :returns: The queried Exo-MerCat helper table.
    :rtype: astropy.table.table.Table
    """
    exo["provider"] = create_provider_table(
        "Exo-MerCat",
        "http://archives.ia2.inaf.it/vo/tap/projects",
        "2020A&C....3100370A",
    )  # change bibcode once paper on ads
    exo_helptab = query(exo["provider"]["provider_url"][0], adql_query)
    return exo_helptab


def load_exomercat(exo: dict) -> Table:
    """
    Load a previously exported Exo-MerCat helper table from local storage.

    This is used as a fallback when the remote TAP service is not reachable.

    :param exo: Mutable dictionary that will receive the provider table.
        Expected to be the shared provider dict used through this module.
    :type exo: dict[str, Table]
    :returns: The locally loaded Exo-MerCat helper table.
    :rtype: astropy.table.table.Table
    """
    exo["provider"] = create_provider_table(
        "Exo-MerCat",
        "http://archives.ia2.inaf.it/vo/tap/projects",
        "2020A&C....3100370A",
        "2024-12-13",
    )
    exo_helptab = ascii.read(
        Path().additional_data + "exo-mercat13-12-2024_v2.0.csv"
    )
    exo_helptab = stringtoobject(exo_helptab, 3000)
    exo["provider"]["provider_access"] = ["2024-12-13"]
    print("loading exomercat version from", exo["provider"]["provider_access"])
    return exo_helptab


def query_or_load_exomercat() -> tuple[dict, Table]:
    """
    Get the Exo-MerCat helper table either via TAP or from local backup.

    Tries the live TAP query first. If that fails (e.g. server down), a
    shipped CSV snapshot is loaded instead. The returned provider dict is
    the shared container used throughout the Exo provider pipeline.

    :returns: Tuple of (provider dict, helper table).
    :rtype: tuple[dict[str, Table], astropy.table.table.Table]
    """
    exo = empty_dict.copy()
    adql_query = """SELECT *
                  FROM exomercat.exomercat"""

    try:
        exo_helptab = query_exomercat(adql_query, exo)
    except Exception:
        exo_helptab = load_exomercat(exo)
    return exo, exo_helptab


def create_object_main_id(exo_helptab: Table) -> Table:
    """
    Build canonical main identifiers for host and planet objects.

    Constructs:
    - host_main_id: host name augmented with a binary component (if present).
    - planet_main_id: host_main_id + letter (planet label).

    The inputs can contain masked values (astropy masked columns). This
    function preserves those semantics and falls back to the 'host'
    column if 'main_id' is masked.

    :param exo_helptab: Raw Exo-MerCat helper table with columns 'main_id',
        'host', 'binary' and 'letter'.
    :type exo_helptab: astropy.table.table.Table
    :returns: The same table with 'host_main_id' and 'planet_main_id' added.
    :rtype: astropy.table.table.Table
    """
    # Initialize target columns with proper length and dtype.
    exo_helptab["planet_main_id"] = Column(  # object dtype for mixed strings
        dtype=object, length=len(exo_helptab)
    )
    # host_main_id starts with the original main_id (may include mask).
    exo_helptab["host_main_id"] = exo_helptab["main_id"]

    for i in range(len(exo_helptab)):
        # If 'main_id' is present, use it, else fallback to 'host'.
        # Note: astropy masked values use MaskedConstant.
        if type(exo_helptab["main_id"][i]) != np.ma.core.MaskedConstant:
            hostname = exo_helptab["main_id"][i]
        else:
            hostname = exo_helptab["host"][i]

        # If a binary code is present, append it (e.g., 'A', 'B', ...).
        if type(exo_helptab["binary"][i]) != np.ma.core.MaskedConstant:
            exo_helptab["host_main_id"][i] = (
                hostname + " " + exo_helptab["binary"][i]
            )
        else:
            exo_helptab["host_main_id"][i] = hostname

        # Planet main id is host_main_id + planet letter (e.g., 'b', 'c').
        exo_helptab["planet_main_id"][i] = (
            exo_helptab["host_main_id"][i] + " " + exo_helptab["letter"][i]
        )
    return exo_helptab


def _distance_cut_and_strip(exo_helptab: Table) -> Table:
    """
    Apply distance cut (via SIMBAD join) and strip whitespace from IDs.

    Whitespace stripping is performed after the distance_cut call to avoid
    issues with masked values during the join.

    :param exo_helptab: Helper table after planet/host ids were created.
    :type exo_helptab: astropy.table.table.Table
    :returns: Table filtered to the distance-limited sample with trimmed
        'planet_main_id', 'main_id', and 'exomercat_name'.
    :rtype: astropy.table.table.Table
    """
    exo_helptab = distance_cut(exo_helptab, "main_id")

    for i in range(len(exo_helptab)):
        exo_helptab["planet_main_id"][i] = exo_helptab["planet_main_id"][i].strip()
        exo_helptab["main_id"][i] = exo_helptab["main_id"][i].strip()
        exo_helptab["exomercat_name"][i] = exo_helptab["exomercat_name"][i].strip()

    return exo_helptab


def _fetch_sim_names(exo_helptab: Table) -> Table:
    """
    Attach SIMBAD planet names to the helper table using a left join.

    A separate SIMBAD lookup is performed on the planet main ids to
    resolve the SIMBAD main id that should be used. The join is left-typed
    to preserve rows not present in SIMBAD.
     TBD: check if I loose them already in the distance_cut?

    :param exo_helptab: Distance-cut and stripped helper table.
    :type exo_helptab: astropy.table.table.Table
    :returns: Table with an added 'sim_planet_main_id' column (may be empty
        for rows not found in SIMBAD).
    :rtype: astropy.table.table.Table
    """
    # Note: include 'host_main_id' to keep a Table (not Column) in fetch call.
    sim_matches = fetch_main_id(
        exo_helptab["planet_main_id", "host_main_id"],
        IdentifierCreator(  # id_creator config for SIMBAD-planet lookup
            name="sim_planet_main_id", colname="planet_main_id"
        ),
    )
    exo_helptab = join(
        exo_helptab,
        sim_matches["sim_planet_main_id", "planet_main_id"],
        keys="planet_main_id",
        join_type="left",
    )
    return exo_helptab


def _compute_removed_objects_and_save(before: Table, after: Table) -> None:
    """
    Persist the list of Exo-MerCat objects removed during distance cut.

    Compares two helper tables by 'exomercat_name' and writes differences
    to a sidecar file for later inspection.

    :param before: Helper table before the distance cut.
    :type before: astropy.table.table.Table
    :param after: Helper table after the distance cut and SIMBAD join.
    :type after: astropy.table.table.Table
    :returns: Nothing. Writes the 'exomercat_removed_objects' table via save().
    :rtype: None
    """
    before["exomercat_name"] = before["exomercat_name"].astype(object)
    removed_objects = setdiff(before, after, keys=["exomercat_name"])
    save([removed_objects], ["exomercat_removed_objects"])
    return



def create_exo_helpertable() -> tuple[dict, Table]:
    """
    Create the Exo-MerCat helper table used downstream.

    Steps:
    1. Query or load Exo-MerCat base data (provider dict + raw table).
    2. Build host/planet canonical ids (host_main_id, planet_main_id).
    3. Apply distance cut (via SIMBAD join) and strip whitespace.
    4. Fetch SIMBAD planet names and left-join them.
    5. Save a table of removed objects for auditing.
    6. Return provider dict and finalized helper table.

    :returns: Tuple of (provider dict, helper table).
    :rtype: tuple[dict[str, Table], astropy.table.table.Table]
    """
    exo, exo_helptab = query_or_load_exomercat()
    exo_helptab = create_object_main_id(exo_helptab)

    # Keep pre-cut copy to report removed objects later.
    exo_helptab_before_distance_cut = exo_helptab

    # Distance cut and id trimming.
    exo_helptab = _distance_cut_and_strip(exo_helptab)

    # Attach SIMBAD main ids for planets via left join.
    exo_helptab = _fetch_sim_names(exo_helptab)

    # Report dropped rows (by exomercat_name) for traceability.
    _compute_removed_objects_and_save(
        exo_helptab_before_distance_cut, exo_helptab
    )

    return exo, exo_helptab


def create_ident_table(exo_helptab: Table, exo: dict[str, Table]) -> tuple[Table, Table]:
    """
    Create the identifier table for Exo-MerCat objects.

    For each row in the helper table this function emits two identifier rows:
    - (main_id, main_id, id_ref) where id_ref comes from SIMBAD if
      sim_planet_main_id is available, otherwise from Exo-MerCat.
    - (main_id, exomercat_name, exomercat_bibcode) as the Exo-MerCat alias.

    If sim_planet_main_id is not empty, planet_main_id in the helper table is
    replaced with that SIMBAD value. This matches downstream expectations
    (e.g. uniqueness, join behavior).

    :param exo_helptab: Helper table with 'planet_main_id', 'sim_planet_main_id'
        and 'exomercat_name' columns.
    :type exo_helptab: astropy.table.table.Table
    :param exo: Provider dict containing the 'provider' table with
        'provider_bibcode'.
    :type exo: dict[str, astropy.table.table.Table]
    :returns: Tuple of (identifier table, possibly updated helper table).
    :rtype: tuple[astropy.table.table.Table, astropy.table.table.Table]
    """
    exo_ident = Table(
        names=["main_id", "id", "id_ref"], dtype=[object, object, object]
    )

    for i in range(len(exo_helptab)):
        # Use SIMBAD name if present, otherwise keep Exo-MerCat one.
        if exo_helptab["sim_planet_main_id"][i] != "":
            main_id = exo_helptab["sim_planet_main_id"][i]
            exo_helptab["planet_main_id"][i] = exo_helptab[
                "sim_planet_main_id"
            ][i]
            # SIMBAD bibcode from local snapshot
            [simbad_ref] = load(["sim_provider"])
            idref = simbad_ref["provider_bibcode"][0]
        else:
            main_id = exo_helptab["planet_main_id"][i]
            # in else because otherwise will result in double sim main id
            # from sim provider in case of sim being main id
            # Exo-MerCat bibcode
            idref = exo["provider"]["provider_bibcode"][0]

        # main_id as an identifier (self-id)
        exo_ident.add_row([main_id, main_id, idref])
        # Exo-MerCat name as alias
        exo_ident.add_row(
            [
                main_id,
                exo_helptab["exomercat_name"][i],
                exo["provider"]["provider_bibcode"][0],
            ]
        )

    exo_ident = unique(exo_ident)
    return exo_ident, exo_helptab


def create_objects_table(exo: dict[str, Table]) -> Table:
    """
    Create the objects table from identifiers.

    Builds a unique list of objects from the identifier table and assigns
    object type 'pl' (planet) for all entries produced by the Exo-MerCat
    provider.

    :param exo: Provider dict that contains 'ident' table with columns
        'main_id' and 'id'.
    :type exo: dict[str, astropy.table.table.Table]
    :returns: Objects table with columns 'main_id', 'ids', and 'type'.
    :rtype: astropy.table.table.Table
    """
    # Note: hosts not in SIMBAD are not added (tbd in future).
    exo_objects = Table(names=["main_id", "ids"], dtype=[object, object])
    exo_objects = ids_from_ident(exo["ident"]["main_id", "id"], exo_objects)
    exo_objects["type"] = ["pl" for _ in range(len(exo_objects))]
    return exo_objects


def deal_with_mass_nullvalues(exo_helptab, cols):
    for colname in cols:
        exo_helptab = nullvalues(exo_helptab, colname, 1e20, verbose=False)
        exo_helptab = replace_value(exo_helptab, colname, np.inf, 1e20)
    return exo_helptab


def create_para_exo_mes_mass_pl(exo_helptab, para, sini_flag):
    table = exo_helptab[
        "planet_main_id",
        para,
        para + "_max",
        para + "_min",
        para + "_url",
        para + "_pl_qual",
        "bestmass_provenance",
    ]
    table.rename_columns(
        ["planet_main_id", para, para + "_url", para + "_max", para + "_min"],
        [
            "main_id",
            "mass_pl_value",
            "mass_pl_ref",
            "mass_pl_err_max",
            "mass_pl_err_min",
        ],
    )
    if para == "msini":
        table.rename_column(para + "_pl_qual", "mass_pl_qual")
    table["mass_pl_sini_flag"] = [sini_flag for j in range(len(table))]
    # remove null values
    table = table[np.where(table["mass_pl_value"] != 1e20)]
    return table


def betterthan(qual1, qual2):
    # is it usefull to do all this work with the quality when I could just assign b for bestmass and c for otherone?
    # nearly done, finish this then do other one as well. call one exomercatqual other dbqual
    quals = ["A", "B", "C", "D", "E", "?"]
    for i in [0, 1, 2, 3, 4]:
        if qual1 == quals[i] and qual2 in quals[i + 1 :]:
            result = True
            return result
        else:
            result = False
    return result


def bestmass_better_qual(bestmass, qual_msini, qual_mass):
    if bestmass == "Mass":
        if betterthan(qual_msini, qual_mass) or qual_msini == qual_mass:
            qual_msini = lower_quality(qual_msini)
            qual_msini, qual_mass = bestmass_better_qual(
                bestmass, qual_msini, qual_mass
            )
    elif bestmass == "Msini":
        if betterthan(qual_mass, qual_msini) or qual_msini == qual_mass:
            qual_mass = lower_quality(qual_mass)
            qual_msini, qual_mass = bestmass_better_qual(
                bestmass, qual_msini, qual_mass
            )
    return qual_msini, qual_mass


def assign_new_qual(exo_mes_mass_pl, main_id, flag, new_qual):
    temp = [np.where(exo_mes_mass_pl["main_id"] == main_id)]
    for i in temp[0][0]:
        if exo_mes_mass_pl["mass_pl_sini_flag"][i] == flag:
            exo_mes_mass_pl["mass_pl_qual"][i] = new_qual
    return exo_mes_mass_pl


def align_quality_with_bestmass(exo_mes_mass_pl):
    grouped_mass_pl = exo_mes_mass_pl.group_by("main_id")
    ind = grouped_mass_pl.groups.indices
    for i in range(len(ind) - 1):  # range(len(ind)-1):
        if ind[i + 1] - ind[i] == 2:
            bestmass = grouped_mass_pl["bestmass_provenance"][ind[i]]
            if grouped_mass_pl["mass_pl_sini_flag"][ind[i]] == "True":
                qual_msini = grouped_mass_pl["mass_pl_qual"][ind[i]]
                qual_mass = grouped_mass_pl["mass_pl_qual"][ind[i] + 1]
            else:
                qual_msini = grouped_mass_pl["mass_pl_qual"][ind[i] + 1]
                qual_mass = grouped_mass_pl["mass_pl_qual"][ind[i]]
            new_qual_msini, new_qual_mass = bestmass_better_qual(
                bestmass, qual_msini, qual_mass
            )
            # assign them the right place
            # I think I need a where here
            if new_qual_msini != qual_msini:
                main_id = grouped_mass_pl["main_id"][ind[i]]
                exo_mes_mass_pl = assign_new_qual(
                    exo_mes_mass_pl, main_id, "True", new_qual_msini
                )
            if new_qual_mass != qual_mass:
                main_id = grouped_mass_pl["main_id"][ind[i]]
                exo_mes_mass_pl = assign_new_qual(
                    exo_mes_mass_pl, main_id, "False", new_qual_mass
                )

    return exo_mes_mass_pl


def create_mes_mass_pl_table(exo_helptab):
    """
    Creates planetary mass measurement table.

    :param exo_helptab: Exo-Mercat helper table.
    :type exo_helptab: astropy.table.table.Table
    :returns: Planetary mass measurement table.
    :rtype: astropy.table.table.Table
    """
    # include other parameters r, a, e, i, p, status
    cols = ["mass_max", "mass_min", "mass", "msini", "msini_max", "msini_min"]
    exo_helptab = deal_with_mass_nullvalues(exo_helptab, cols)

    exo_helptab = assign_quality(exo_helptab, special_mode="exo")

    exo_mes_mass_pl1 = create_para_exo_mes_mass_pl(exo_helptab, "mass", "False")
    exo_mes_mass_pl2 = create_para_exo_mes_mass_pl(exo_helptab, "msini", "True")

    exo_mes_mass_pl = vstack([exo_mes_mass_pl1, exo_mes_mass_pl2])
    exo_mes_mass_pl = unique(exo_mes_mass_pl)
    exo_mes_mass_pl = align_quality_with_bestmass(exo_mes_mass_pl)
    return exo_mes_mass_pl


def create_h_link_table(exo_helptab: Table, exo: dict[str, Table]) -> Table:
    """
    Create the hierarchical link table for planet-to-host relations.

    The link connects child 'main_id' (planet) to 'parent_main_id' (host).
    A uniform reference column is filled using the Exo-MerCat provider bibcode.

    :param exo_helptab: Helper table with 'planet_main_id' and 'host_main_id'.
    :type exo_helptab: astropy.table.table.Table
    :param exo: Provider dict containing 'provider' table with
        'provider_bibcode'.
    :type exo: dict[str, astropy.table.table.Table]
    :returns: h_link table with 'main_id', 'parent_main_id', 'h_link_ref'.
    :rtype: astropy.table.table.Table
    """
    exo_h_link = exo_helptab["planet_main_id", "host_main_id"]
    exo_h_link.rename_columns(
        ["planet_main_id", "host_main_id"], ["main_id", "parent_main_id"]
    )
    exo_h_link["h_link_ref"] = [
        exo["provider"]["provider_bibcode"][0] for _ in range(len(exo_h_link))
    ]
    return exo_h_link


def create_exo_sources_table(exo):
    """
    Creates sources table.

    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Sources table.
    :rtype: astropy.table.table.Table
    """
    tables = [exo["provider"], exo["h_link"], exo["ident"], exo["mes_mass_pl"]]
    # define header name of columns containing references data
    ref_columns = [
        ["provider_bibcode"],
        ["h_link_ref"],
        ["id_ref"],
        ["mass_pl_ref"],
    ]
    exo_sources = create_sources_table(
        tables, ref_columns, exo["provider"]["provider_name"][0]
    )
    return exo_sources


def provider_exo():
    """
    This function obtains and arranges exomercat data.

    :returns: Dictionary with names and astropy tables containing
        reference data, object data, identifier data, object to object
        relation data, basic planetary data and planetary mass measurement
        data.
    :rtype: dict(str,astropy.table.table.Table)
    """
    exo, exo_helptab = create_exo_helpertable()

    exo["ident"], exo_helptab = create_ident_table(exo_helptab, exo)
    exo["objects"] = create_objects_table(exo)
    exo["mes_mass_pl"] = create_mes_mass_pl_table(exo_helptab)
    exo["h_link"] = create_h_link_table(exo_helptab, exo)
    exo["sources"] = create_exo_sources_table(exo)

    save(list(exo.values()), ["exo_" + element for element in list(exo.keys())])
    return exo
