"""
Generates the data for the database for the provider LIFE.
"""

import numpy as np  # arrays
from astropy import coordinates, units
from astropy.io import ascii, votable
from astropy.table import Column, MaskedColumn, join, unique
from provider.assign_quality_funcs import assign_quality
from provider.utils import (
    create_provider_table,
    create_sources_table,
    initiate_columns,
    replace_value,
)
from sdata import empty_dict
from utils.io import Path, load, save


def extract_lum_class(nr, sptype):
    """
    Extracts luminocity class.

    :param int nr: index number
    :param str sptype: spectral type
    :returns: luminocity class
    :rtype: str
    """
    lum_class = sptype[nr]
    if len(sptype) > nr + 1 and sptype[nr + 1] in ["I", "V"]:
        lum_class = sptype[nr : nr + 2]
        if len(sptype) > nr + 2 and sptype[nr + 2] in ["I", "V"]:
            lum_class = sptype[nr : nr + 3]
    return lum_class


def assign_null_values(table, index):
    """
    Sets table entry to '?' for given index in specific columns.

    :param table: Table containing columns class_temp, class_temp_nr,
        class_lum and class_ref.
    :type table: astropy.table.table.Table
    :param int index: Index of entry.
    """
    table["class_temp"][index] = "?"
    table["class_temp_nr"][index] = "?"
    table["class_lum"][index] = "?"
    table["class_ref"][index] = "?"
    return table


def deal_with_leading_d_sptype(table, index):
    """
    Deals with old annotation of leading d in spectraltype representing dwarf star.

    :param table: Table containing column 'sptype_string'.
    :type table: astropy.table.table.Table
    :param int index: Index of entry.
    :returns: Entry with leading d removed.
    :rtype: str
    """
    sptype = table["sptype_string"][index]
    if len(sptype) > 0:
        if sptype[0] == "d":
            table["class_lum"][index] = "V"
            sptype = sptype.strip("d")
    return sptype


def deal_with_middle_minus(sptype_string):
    """
    Removes intermediate minus sign from sptype_string.

    :param str sptype_string:
    :returns: Entry with middle - removed.
    :rtype: str
    """
    if len(sptype_string) > 0:
        sp_type = sptype_string.split("-")
        if len(sp_type) > 1:
            sptype_string = "".join(sp_type)
    return sptype_string


def decimal_sptype(i, sptype, table):
    table["class_temp_nr"][i] = sptype[1:4]
    if len(sptype) > 4 and sptype[4] in ["I", "V"]:
        table["class_lum"][i] = extract_lum_class(4, sptype)
    else:
        table["class_lum"][i] = "V"  # assumption
    return table


def assign_diff_lum_classes(i, sptype, table):
    if len(sptype) > 1 and sptype[1] in [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
    ]:
        table["class_temp_nr"][i] = sptype[1]
        # distinguishing between objects like K5V and K5.5V
        if len(sptype) > 2 and sptype[2] == ".":
            table = decimal_sptype(i, sptype, table)
        elif len(sptype) > 2 and sptype[2] in ["I", "V"]:
            table["class_lum"][i] = extract_lum_class(2, sptype)
        else:
            table["class_lum"][i] = "V"  # assumption
    else:
        table["class_lum"][i] = "V"  # assumption
    return table


def sptype_string_to_class(table, ref):
    """
    Extracts stellar parameters from spectral type string one.

    This function extracts the temperature class, temperature class number
    and luminocity class information from the spectral type string (e.g.
    M5V to M, 5 and V). It stores that information in the for this purpose
    generated new columns. Only objects of temperature class O, B, A, F,
    G, K, and M are processed. Only objects of luminocity class IV, V and VI
    are processed.

    :param table: Table containing spectral type information in
        the column sptype_string.
    :type table: astropy.table.table.Table
    :param str ref: Designates origin of data.
    :returns: Table like temp with additional columns class_temp,
        class_temp_nr, class_lum and class_ref.
    :rtype: astropy.table.table.Table
    """
    table = initiate_columns(
        table,
        ["class_temp", "class_temp_nr", "class_lum", "class_ref"],
        [object, object, object, object],
        [True, True, True, True],
    )

    for i in range(len(table)):
        # sorting out objects like M5V+K7V
        sptype = deal_with_leading_d_sptype(table, i)
        sptype = deal_with_middle_minus(sptype)

        if (
            len(sptype.split("+")) == 1
            and
            # sorting out entries like ''
            len(sptype) > 0
            and
            # sorting out brown dwarfs i.e. T1V
            sptype[0] in ["O", "B", "A", "F", "G", "K", "M"]
        ):
            # assigning temperature class and reference
            table["class_temp"][i] = sptype[0]
            table["class_ref"][i] = ref
            table = assign_diff_lum_classes(i, sptype, table)
        else:
            table = assign_null_values(table, i)
        table["sptype_string"][i] = sptype
    return table


def realspectype(cat):
    """
    Removes rows not containing main sequence stars.

    Removes rows of cat where elements in column named 'sim_sptype' are
    either '', 'nan' or start with another letter than the main sequence
    spectral type classification.

    :param cat: Table containing class_temp and class_lum column
    :type cat: astropy.table.table.Table
    :returns: Table, param cat with undesired rows removed
    :rtype: astropy.table.table.Table
    """
    ms_tempclass = np.array(["O", "B", "A", "F", "G", "K", "M"])
    ms_temp = cat[np.where(np.isin(cat["class_temp"], ms_tempclass))]

    ms_lumclass = np.array(["V"])
    ms = ms_temp[np.where(np.isin(ms_temp["class_lum"], ms_lumclass))]

    return ms


def modeled_param():
    """
    Loads and cleans up model file.

    Loads the table of Eric E. Mamajek containing stellar parameters
    modeled from spectral types. Cleans up the columns for spectral
    type, effective temperature radius and mass.

    :returns: Table of the 4 parameters as columns
    :rtype: astropy.table.table.Table
    """
    eem_table_all_columns = ascii.read(
        Path().additional_data + "Mamajek2022-04-16.csv"
    )
    eem_table = eem_table_all_columns["SpT", "Teff", "R_Rsun", "Msun"]
    eem_table.rename_columns(["R_Rsun", "Msun"], ["Radius", "Mass"])
    eem_table = replace_value(eem_table, "Radius", " ...", "nan")
    eem_table = replace_value(eem_table, "Mass", " ...", "nan")
    eem_table = replace_value(eem_table, "Mass", " ....", "nan")
    eem_table["Teff"].unit = units.K
    eem_table["Radius"].unit = units.R_sun
    eem_table["Mass"].unit = units.M_sun
    votable.writeto(
        votable.from_table(eem_table),
        f"{Path().additional_data}model_param.xml",
    )  # saving votable
    return eem_table


def match_sptype(
    cat,
    sptypestring="mp_specmatch",
    teffstring="mod_Teff",
    rstring="mod_R",
    mstring="mod_M",
):
    """
    Assigns modeled parameter values.

    Matches the spectral types with the ones in Mamajek's table and
    includes the modeled effective Temperature,
    stellar radius and stellar mass into the catalog.

    :param cat: astropy table containing spectral type information
    :type cat: astropy.table.table.Table
    :param str sptypestring: Column name where the spectral
        type information is located
    :param str teffstring: Column name for the stellar effective
        temperature column
    :param str rstring: Column name for the stellar radius column
    :param str mstring: Column name for the stellar mass column
    :returns: Table cat with added new columns for
        effective temperature, radius and mass filled with model values
    :rtype: astropy.table.table.Table
    """
    model_param = modeled_param()  # Load Mamajek's model table

    # Initialize columns with proper units and masked arrays
    num_rows = len(cat)
    cat[teffstring] = MaskedColumn(
        mask=np.full(num_rows, True), length=num_rows, unit=units.K
    )
    cat[rstring] = MaskedColumn(
        mask=np.full(num_rows, True), length=num_rows, unit=units.R_sun
    )
    cat[mstring] = MaskedColumn(
        mask=np.full(num_rows, True), length=num_rows, unit=units.M_sun
    )

    # Process each spectral type in the catalog
    for j in range(num_rows):
        sptype = cat[sptypestring][j]

        # Handle empty spectral types
        if sptype == "":
            cat[sptypestring][j] = "None"
            continue

        # Remove old 'd' notation (dwarf = main sequence star)
        sptype = sptype.strip("d")
        cat[sptypestring][j] = sptype

        # Try to match spectral type with model
        matched = _match_spectral_type_to_model(sptype, model_param)

        if matched:
            cat[teffstring][j] = matched["Teff"]
            cat[rstring][j] = matched["Radius"]
            cat[mstring][j] = matched["Mass"]

    return cat


def _match_spectral_type_to_model(sptype, model_param):
    """
    Helper function to match a spectral type string to Mamajek's model.

    :param str sptype: Spectral type string (e.g., 'G2V', 'K5.5V')
    :param model_param: Model parameter table from Mamajek
    :returns: Dictionary with matched parameters or None if no match found
    :rtype: dict or None
    """
    # First try exact match on first 3 characters (e.g., 'G2V')
    for i in range(len(model_param["SpT"])):
        if model_param["SpT"][i][:3] == sptype[:3]:
            return {
                "Teff": model_param["Teff"][i],
                "Radius": model_param["Radius"][i],
                "Mass": model_param["Mass"][i]
            }

    # For half-subtypes (e.g., 'K5.5V'), try matching first 4 characters
    # The model doesn't cover all spectral types with .5 accuracy
    if len(sptype) >= 4 and sptype[2:4] == ".5":
        for i in range(len(model_param["SpT"])):
            if model_param["SpT"][i][:4] == sptype[:4]:
                return {
                    "Teff": model_param["Teff"][i],
                    "Radius": model_param["Radius"][i],
                    "Mass": model_param["Mass"][i]
                }

    # No match found
    return None


def spec(cat):
    """
    Runs the spectral type related functions realspectype and match_sptype.

    It also removes all empty columns of the effective temperature, removes
    rows that are not main sequence, removes rows with binary subtype and
    non unique simbad name.

    :param cat: astropy table containing columns named
        'main_id', class_temp_nr, class_temp and,class_lum
    :type cat: astropy.table.table.Table
    :returns: Catalog of mainsequence stars with unique
        simbad names, no binary subtypes and modeled parameters.
    :rtype: astropy.table.table.Table
    """
    # Do I even need realspectype function? I can just take cat where class_temp not empty
    cat = realspectype(cat)
    # model_param=votable.parse_single_table(\
    # f"catalogs/model_param.xml").to_table()
    cat = cat[np.where(cat["class_temp_nr"] != 0)]
    cat["specmatch_temp_nr"] = cat["class_temp_nr"]
    for i, temp_nr in enumerate(cat["specmatch_temp_nr"]):
        if temp_nr[1:3] == ".0":
            cat["specmatch_temp_nr"][i] = temp_nr[0]
    cat["mp_specmatch"] = (
        cat["class_temp"] + cat["specmatch_temp_nr"] + cat["class_lum"]
    )
    cat = match_sptype(cat)
    cat.remove_rows([np.where(cat["mod_Teff"].mask == True)])
    cat.remove_rows([np.where(np.isnan(cat["mod_Teff"]))])
    cat = unique(cat, keys="main_id")
    return cat


def create_star_basic_table():
    """
    Creates basic stellar data table.

    :returns: Dictionary of database table names and tables with
        filled basic stellar data table.
    :rtype: dict(str,astropy.table.table.Table)
    """
    life = empty_dict.copy()
    life["provider"] = create_provider_table(
        "LIFE", "www.life-space-mission.com", "2022A&A...664A..21Q"
    )

    # galactic coordinates:  transformed from simbad ircs coordinates using astropy
    [life_star_basic] = load(["sim_star_basic"])
    ircs_coord = coordinates.SkyCoord(
        ra=life_star_basic["coo_ra"],
        dec=life_star_basic["coo_dec"],
        frame="icrs",
    )
    gal_coord = ircs_coord.galactic
    life_star_basic["coo_gal_l"] = gal_coord.l.deg * units.degree
    life_star_basic["coo_gal_b"] = gal_coord.b.deg * units.degree
    life_star_basic["dist_st_value"] = 1000.0 / life_star_basic["plx_value"]
    life_star_basic["dist_st_value"] = np.round(
        life_star_basic["dist_st_value"], 2
    )
    # null value treatment: plx_value has masked entries therefore distance_values too
    # ref:
    life_star_basic["dist_st_ref"] = MaskedColumn(
        dtype=object,
        length=len(life_star_basic),
        mask=[True for j in range(len(life_star_basic))],
    )
    life_star_basic["dist_st_ref"][
        np.where(life_star_basic["dist_st_value"].mask == False)
    ] = [
        life["provider"]["provider_name"][0]
        for j in range(
            len(
                life_star_basic["dist_st_ref"][
                    np.where(life_star_basic["dist_st_value"].mask == False)
                ]
            )
        )
    ]
    # can I do the same transformation with the errors? -> try on some examples and compare to simbad ones
    life_star_basic["coo_gal_err_angle"] = [
        -1 for j in range(len(life_star_basic))
    ]
    life_star_basic["coo_gal_err_maj"] = [
        -1 for j in range(len(life_star_basic))
    ]
    life_star_basic["coo_gal_err_min"] = [
        -1 for j in range(len(life_star_basic))
    ]
    life_star_basic = assign_quality(life_star_basic, "coo_gal_qual")
    life_star_basic["main_id"] = life_star_basic["main_id"].astype(str)
    # source
    # transformed from simbad ircs coordinates using astropy
    life_star_basic["coo_gal_ref"] = Column(
        dtype=object, length=len(life_star_basic)
    )
    life_star_basic["coo_gal_ref"] = life["provider"]["provider_name"][0]
    # for all entries since coo_gal column not masked column

    life_star_basic["coo_gal_ref"] = life_star_basic["coo_gal_ref"].astype(str)
    life_star_basic = life_star_basic[
        "main_id",
        "coo_gal_l",
        "coo_gal_b",
        "coo_gal_err_angle",
        "coo_gal_err_maj",
        "coo_gal_err_min",
        "coo_gal_qual",
        "coo_gal_ref",
        "dist_st_value",
        "dist_st_ref",
        "sptype_string",
    ]

    life_star_basic = sptype_string_to_class(
        life_star_basic, life["provider"]["provider_name"][0]
    )
    life["star_basic"] = life_star_basic
    return life


def create_life_helpertable(life):
    """
    Creates helper table.

    :param life: Dictionary of database table names and tables.
    :type life: dict(str,astropy.table.table.Table)
    :returns: Helper table.
    :rtype: astropy.table.table.Table
    """
    # applying model from E. E. Mamajek on SIMBAD spectral type

    [sim_objects] = load(["sim_objects"], stringtoobjects=False)

    stars = sim_objects[np.where(sim_objects["type"] == "st")]
    life_helptab = join(stars, life["star_basic"])
    life_helptab = spec(
        life_helptab[
            "main_id",
            "sptype_string",
            "class_lum",
            "class_temp",
            "class_temp_nr",
        ]
    )
    # if I take only st objects from sim_star_basic I don't loose objects during realspectype
    return life_helptab


def create_mes_teff_st_table(life_helptab):
    """
    Creates stellar effective temperature table.

    :param life_helptab: Life helper table.
    :type life_helptab: astropy.table.table.Table
    :returns: Stellar effective temperature table.
    :rtype: astropy.table.table.Table
    """
    life_mes_teff_st = life_helptab["main_id", "mod_Teff"]
    life_mes_teff_st.rename_column("mod_Teff", "teff_st_value")
    life_mes_teff_st = assign_quality(life_mes_teff_st, "teff_st_qual", "model")
    life_mes_teff_st["teff_st_ref"] = [
        "2013ApJS..208....9P" for i in range(len(life_mes_teff_st))
    ]
    return life_mes_teff_st


def create_mes_radius_st_table(life_helptab):
    """
    Creates stellar radius table.

    :param life_helptab: Life helper table.
    :type life_helptab: astropy.table.table.Table
    :returns: Stellar radius table.
    :rtype: astropy.table.table.Table
    """
    life_mes_radius_st = life_helptab["main_id", "mod_R"]
    life_mes_radius_st.rename_column("mod_R", "radius_st_value")
    life_mes_radius_st = assign_quality(
        life_mes_radius_st, "radius_st_qual", "model"
    )
    life_mes_radius_st["radius_st_ref"] = [
        "2013ApJS..208....9P" for i in range(len(life_mes_radius_st))
    ]
    return life_mes_radius_st


def create_mes_mass_st_table(life_helptab):
    """
    Creates stellar mass table.

    :param life_helptab: Life helper table.
    :type life_helptab: astropy.table.table.Table
    :returns: Stellar mass table.
    :rtype: astropy.table.table.Table
    """
    life_mes_mass_st = life_helptab["main_id", "mod_M"]
    life_mes_mass_st.rename_column("mod_M", "mass_st_value")
    life_mes_mass_st = assign_quality(life_mes_mass_st, "mass_st_qual", "model")
    life_mes_mass_st["mass_st_ref"] = [
        "2013ApJS..208....9P" for i in range(len(life_mes_mass_st))
    ]

    # specifying stars cocerning multiplicity
    # main sequence simbad object type: MS*, MS? -> luminocity class
    # Interacting binaries and close CPM systems: **, **?
    return life_mes_mass_st


def create_life_sources_table(life):
    """
    Creates sources table.

    :param life: Dictionary of database table names and tables.
    :type life: dict(str,astropy.table.table.Table)
    :returns: Sources table.
    :rtype: astropy.table.table.Table
    """
    tables = [
        life["provider"],
        life["star_basic"],
        life["mes_teff_st"],
        life["mes_radius_st"],
        life["mes_mass_st"],
    ]
    # define header name of columns containing references data
    ref_columns = [
        ["provider_bibcode"],
        ["coo_gal_ref"],
        ["teff_st_ref"],
        ["radius_st_ref"],
        ["mass_st_ref"],
    ]
    life_sources = create_sources_table(
        tables, ref_columns, life["provider"]["provider_name"][0]
    )
    return life_sources


def provider_life():
    """
    Loads SIMBAD data and postprocesses it.

    Postprocessing enables to provide more useful information. It uses a model
    from Eric E. Mamajek to predict temperature, mass and radius from the simbad
    spectral type data.

    :returns: Dictionary with names and astropy tables containing
        reference data, provider data, basic stellar data, stellar effective
        temperature data, stellar radius data and stellar mass data.
    :rtype: dict(str,astropy.table.table.Table)
    """
    life = create_star_basic_table()
    life_helptab = create_life_helpertable(life)
    # removing this column because I had to adapt it where there was a leadin d entry but change not useful for db just for
    # life parameter creation
    life["star_basic"].remove_column("sptype_string")
    life["mes_teff_st"] = create_mes_teff_st_table(life_helptab)
    life["mes_radius_st"] = create_mes_radius_st_table(life_helptab)
    life["mes_mass_st"] = create_mes_mass_st_table(life_helptab)
    life["sources"] = create_life_sources_table(life)

    save(
        list(life.values()),
        ["life_" + element for element in list(life.keys())],
    )
    return life
