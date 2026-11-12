"""
Generates the data for the database for the provider gaia.
"""

import numpy as np  # arrays
from astropy.table import (
    Table,
    join,
    setdiff,
    vstack,
    unique,
    Column
)
from provider.assign_quality_funcs import assign_quality
from provider.utils import (
    create_provider_table,
    create_sources_table,
    ids_from_ident,
    sync_or_async_query,
    replace_value,
    distance_cut,
)
from pyvo.dal import TAPService
from sdata import empty_dict

# self created modules
from utils.io import save, load

gaia_ref = "2016A&A...595A...1G"


def create_gaia_helpertable(distance_cut_in_pc,temp):
    """
    Creates helper table.

    If the gaia server asynchronous query is not working a temporary method to
    use the synchronous query is used.

    :param float distance_cut_in_pc: Distance up to which stars are included.
    :param temp: If ``True``, load cached helper data instead of querying.
    :type temp: bool
    :returns: Helper table.
    :rtype: astropy.table.table.Table
    """
    plx_in_mas_cut = 1000.0 / distance_cut_in_pc
    # making cut a bit bigger for correct treatment of objects on boundary
    plx_cut = plx_in_mas_cut - plx_in_mas_cut / 10.0

    gaia = empty_dict.copy()
    gaia["provider"] = create_provider_table(
        "Gaia", "https://gea.esac.esa.int/tap-server/tap", gaia_ref
    )

    if temp:
        print(" loading...")
        [gaia_helptab] = load(["gaia_helptab_query"])
    else:
        # query
        adql_query = """
        SELECT s.source_id ,p.mass_flame, p.radius_flame,
            p.teff_gspphot, p.teff_gspspec, m.nss_solution_type, p.age_flame,
            p.teff_gspspec_lower, p.teff_gspspec_upper, p.flags_gspspec,
            s.ra AS coo_ra, s.dec AS coo_dec
        FROM gaiadr3.gaia_source as s
            JOIN gaiadr3.astrophysical_parameters as p ON s.source_id=p.source_id
                LEFT JOIN gaiadr3.nss_two_body_orbit as m ON s.source_id=m.source_id
        WHERE s.parallax >=""" + str(plx_in_mas_cut)


        gaia_helptab = sync_or_async_query(
                gaia["provider"]["provider_url"][0], adql_query,
                name = "gaia_helptab_query"
            )


    gaia_helptab.rename_columns(
        ["mass_flame", "radius_flame"], ["mass_st_value", "radius_st_value"]
    )
    gaia_helptab["gaia_id"] = [
        "Gaia DR3 " + str(gaia_helptab["source_id"][j])
        for j in range(len(gaia_helptab))
    ]
    gaia_helptab["ref"] = [
        gaia["provider"]["provider_bibcode"][0] for j in range(len(gaia_helptab))
    ]
    return gaia_helptab, gaia


def create_ident_table(gaia_helptab):
    """
    Creates identifier table.

    :param gaia_helptab: Gaia helper table.
    :type gaia_helptab: astropy.table.table.Table
    :returns: Identifier table and gaia helper table.
    :rtype: astropy.table.table.Table, astropy.table.table.Table
    """
    # ---------------gaia_ident-----------------------
    gaia_helptab["gaia_id"] = gaia_helptab["gaia_id"].astype(object)
    gaia_sim_idmatch = distance_cut(
        gaia_helptab["gaia_id", "ref"],
        colname="gaia_id",main_id = False,
    )
    # should be gaia_id, main_id, ref minus 40 objects that have only gaia_id
    gaia_ident = gaia_sim_idmatch.copy()
    gaia_ident.rename_columns(["gaia_id", "ref"], ["id", "id_ref"])
    gaia_ident["id_ref"] = gaia_ident["id_ref"].astype(object)
    # creating simbad main_id ident rows
    sim_main_id_ident = Table()
    sim_main_id_ident["main_id"] = gaia_ident["main_id"]
    sim_main_id_ident["id"] = gaia_ident["main_id"]
    sim_ref = "2000A&AS..143....9W"
    sim_main_id_ident["id_ref"] = Column(data = [
        sim_ref for _ in range(len(sim_main_id_ident))
        ], dtype=object)
    condition = np.isin(gaia_ident["id"],
                                 sim_main_id_ident["main_id"])
    gaia_ident["id_ref"][condition] = (
        [sim_ref for _ in range(len(gaia_ident[condition]))])
    gaia_ident = vstack([gaia_ident, sim_main_id_ident])

    gaia_ident = unique(gaia_ident)

    # now need to add the 40 objects that have only gaia_identifiers
    # for setdiff need both columns to be same type
    for col in ["gaia_id","ref"]:
        gaia_sim_idmatch[col] = gaia_sim_idmatch[col].astype(object)
        gaia_helptab[col] =gaia_helptab[col].astype(object)
    gaia_only_id = setdiff(
        gaia_helptab["gaia_id", "ref"], gaia_sim_idmatch["gaia_id", "ref"]
    )
    gaia_only_id["main_id"] = gaia_only_id["gaia_id"]
    gaia_only_id.rename_columns(["gaia_id", "ref"], ["id", "id_ref"])
    gaia_only_id = gaia_only_id[gaia_only_id["id_ref"] == gaia_ref]
    # for vstack need both columns to be same type
    for col in gaia_ident.colnames:
        gaia_ident[col] = gaia_ident[col].astype(object)
    gaia_ident = vstack([gaia_ident, gaia_only_id])
    # add main_id to gaia table
    gaia_helptab = join(
        gaia_ident["main_id", "id"],
        gaia_helptab,
        keys_left="id",
        keys_right="gaia_id",
    )
    gaia_helptab.remove_column("id")
    return gaia_ident, gaia_helptab


def create_objects_table(gaia_helptab, gaia):
    """
    Creates objects table.

    :param gaia_helptab: Gaia helper table.
    :type gaia_helptab: astropy.table.table.Table
    :param gaia: Dictionary of database table names and tables.
    :type gaia: dict(str,astropy.table.table.Table)
    :returns: Objects table.
    :rtype: astropy.table.table.Table
    """
    # -----------------gaia_objects------------------
    gaia_objects = Table(names=["main_id", "ids"], dtype=[object, object])
    gaia_objects = ids_from_ident(gaia["ident"]["main_id", "id"], gaia_objects)
    gaia_objects["type"] = ["st" for j in range(len(gaia_objects))]
    gaia_objects["main_id"] = gaia_objects["main_id"].astype(object)
    gaia_objects = join(
        gaia_objects,
        gaia_helptab["main_id", "nss_solution_type"],
        join_type="left",
    )
    gaia_objects["type"][np.where(gaia_objects["nss_solution_type"] != "")] = [
        "sy"
        for j in range(
            len(
                gaia_objects["type"][
                    np.where(gaia_objects["nss_solution_type"] != "")
                ]
            )
        )
    ]
    gaia_objects.remove_column("nss_solution_type")
    return gaia_objects


def create_mes_binary_table(gaia):
    """
    Creates binarity measurement table.

    :param gaia: Dictionary of database table names and tables.
    :type gaia: dict(str,astropy.table.table.Table)
    :returns: Binarity measurement table.
    :rtype: astropy.table.table.Table
    """
    gaia_mes_binary = gaia["objects"]["main_id", "type"]
    # tbd add binary flag True to children of system objects once I get h_link
    # info from gaia
    gaia_mes_binary.rename_column("type", "binary_flag")
    gaia_mes_binary["binary_flag"] = gaia_mes_binary["binary_flag"].astype(
        object
    )
    gaia_mes_binary = replace_value(
        gaia_mes_binary, "binary_flag", "sy", "True"
    )
    gaia_mes_binary = replace_value(
        gaia_mes_binary, "binary_flag", "st", "False"
    )
    gaia_mes_binary["binary_ref"] = [
        gaia["provider"]["provider_bibcode"][0] for j in range(len(gaia_mes_binary))
    ]
    gaia_mes_binary = assign_quality(
        gaia_mes_binary, "binary_qual", special_mode="gaia_binary"
    )
    # if necessary lower binary_qual for binary_flag = False to level of simbad.
    return gaia_mes_binary


def create_mes_teff_st_table(gaia_helptab):
    """
    Creates stellar effective temperature table.

    :param gaia_helptab: Gaia helper table.
    :type gaia_helptab: astropy.table.table.Table
    :returns: Stellar effective temperature table.
    :rtype: astropy.table.table.Table
    """
    gaia_mes_teff_st = gaia_helptab["main_id", "teff_gspphot"]
    gaia_mes_teff_st["ref"] = [
        gaia_helptab["ref"][j] + " GSP-Phot" for j in range(len(gaia_helptab))
    ]

    # remove masked rows
    # gsp_phot
    gaia_mes_teff_st.remove_rows(
        gaia_mes_teff_st["teff_gspphot"].mask.nonzero()[0]
    )
    gaia_mes_teff_st.rename_columns(
        ["teff_gspphot", "ref"], ["teff_st_value", "teff_st_ref"]
    )
    gaia_mes_teff_st = assign_quality(
        gaia_mes_teff_st, "teff_st_qual", special_mode="teff_st_phot"
    )

    # gsp_spec
    gaia_mes_teff_st_spec = gaia_helptab[
        "main_id", "teff_gspspec", "flags_gspspec"
    ]
    gaia_mes_teff_st_spec["ref"] = [
        gaia_helptab["ref"][j] + " GSP-Spec" for j in range(len(gaia_helptab))
    ]
    gaia_mes_teff_st_spec.remove_rows(
        gaia_mes_teff_st_spec["teff_gspspec"].mask.nonzero()[0]
    )
    gaia_mes_teff_st_spec = assign_quality(
        gaia_mes_teff_st_spec, special_mode="teff_st_spec"
    )
    gaia_mes_teff_st_spec.rename_columns(
        ["teff_gspspec", "ref"], ["teff_st_value", "teff_st_ref"]
    )

    # combine phot and spec
    gaia_mes_teff_st = vstack([gaia_mes_teff_st, gaia_mes_teff_st_spec])
    gaia_mes_teff_st = gaia_mes_teff_st[
        "main_id", "teff_st_value", "teff_st_qual", "teff_st_ref"
    ]
    return gaia_mes_teff_st


def create_mes_radius_st_table(gaia_helptab):
    """
    Creates stellar radius table.

    :param gaia_helptab: Gaia helper table.
    :type gaia_helptab: astropy.table.table.Table
    :returns: Stellar radius table.
    :rtype: astropy.table.table.Table
    """
    gaia_mes_radius_st = gaia_helptab["main_id", "radius_st_value", "ref"]
    gaia_mes_radius_st.remove_rows(
        gaia_mes_radius_st["radius_st_value"].mask.nonzero()[0]
    )
    gaia_mes_radius_st = assign_quality(
        gaia_mes_radius_st, "radius_st_qual", "radius_st_flame"
    )
    gaia_mes_radius_st["radius_st_ref"] = [
        gaia_mes_radius_st["ref"][j] + " FLAME"
        for j in range(len(gaia_mes_radius_st))
    ]
    gaia_mes_radius_st.remove_column("ref")
    return gaia_mes_radius_st


def create_mes_mass_st_table(gaia_helptab):
    """
    Creates stellar mass table.

    :param gaia_helptab: Gaia helper table.
    :type gaia_helptab: astropy.table.table.Table
    :returns: Stellar masse table.
    :rtype: astropy.table.table.Table
    """
    gaia_mes_mass_st = gaia_helptab["main_id", "mass_st_value", "ref"]
    gaia_mes_mass_st.remove_rows(
        gaia_mes_mass_st["mass_st_value"].mask.nonzero()[0]
    )
    gaia_mes_mass_st = assign_quality(
        gaia_mes_mass_st, "mass_st_qual", "mass_st_flame"
    )
    gaia_mes_mass_st["mass_st_ref"] = [
        gaia_mes_mass_st["ref"][j] + " FLAME"
        for j in range(len(gaia_mes_mass_st))
    ]
    gaia_mes_mass_st.remove_column("ref")
    return gaia_mes_mass_st

def create_gaia_star_basic_table(gaia_helptab,gaia):
    gaia_star_basic = gaia_helptab["main_id", "coo_ra", "coo_dec"]
    empty_column_object_type = Column(data=[
        "" for _ in range(len(gaia_star_basic))], dtype=object)

    gaia_star_basic["coo_ref"] = empty_column_object_type
    gaia_star_basic["coo_qual"] = empty_column_object_type
    # now wherever coo_ra is not masked, fill with provider bibcode
    mask_has_val = np.where(gaia_star_basic["coo_ra"].mask == False)
    gaia_star_basic[mask_has_val] = replace_value(
        gaia_star_basic[mask_has_val],
        "coo_ref",
        "",
        gaia["provider"]["provider_bibcode"][0],
    )
    gaia_star_basic[mask_has_val] = replace_value(
        gaia_star_basic[mask_has_val],
        "coo_qual",
        "",
        "B",
    )

    return gaia_star_basic


def create_gaia_sources_table(gaia):
    """
    Creates sources table.

    :param gaia: Dictionary of database table names and tables.
    :type gaia: dict(str,astropy.table.table.Table)
    :returns: Sources table.
    :rtype: astropy.table.table.Table
    """
    tables = [
        gaia["provider"],
        gaia["ident"],
        gaia["mes_teff_st"],
        gaia["mes_radius_st"],
        gaia["mes_mass_st"],
        gaia["mes_binary"],
    ]
    # define header name of columns containing references data
    ref_columns = [
        ["provider_bibcode"],
        ["id_ref"],
        ["teff_st_ref"],
        ["radius_st_ref"],
        ["mass_st_ref"],
        ["binary_ref"],
    ]
    gaia_sources = create_sources_table(
        tables, ref_columns, gaia["provider"]["provider_name"][0]
    )
    return gaia_sources


def provider_gaia(distance_cut_in_pc,temp: bool = False):
    """
    Obtains and arranges gaia data.

    :param float distance_cut_in_pc: Distance up to which stars are included.
    :param temp: If ``True``, load cached helper data instead of querying.
    :type temp: bool
    :returns: Dictionary with names and astropy tables containing
        reference data, provider data, object data, identifier data,
        stellar effective temperature, radius, mass and binarity data.
    :rtype: dict(str,astropy.table.table.Table)
    """
    gaia_helptab, gaia = create_gaia_helpertable(distance_cut_in_pc,temp)
    gaia["ident"], gaia_helptab = create_ident_table(gaia_helptab)
    gaia["objects"] = create_objects_table(gaia_helptab, gaia)
    gaia["mes_binary"] = create_mes_binary_table(gaia)
    gaia["mes_teff_st"] = create_mes_teff_st_table(gaia_helptab)
    gaia["mes_radius_st"] = create_mes_radius_st_table(gaia_helptab)
    gaia["mes_mass_st"] = create_mes_mass_st_table(gaia_helptab)
    gaia["star_basic"] = create_gaia_star_basic_table(gaia_helptab,gaia)
    gaia["sources"] = create_gaia_sources_table(gaia)

    save(
        list(gaia.values()),
        ["gaia_" + element for element in list(gaia.keys())],
    )
    return gaia
