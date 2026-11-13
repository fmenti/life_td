"""
Static table and parameter definitions.

This module keeps table schemas in one place and derives the public dictionaries
used elsewhere from those schemas.  The goal is to avoid defining the same table
names, column names, dtypes, and parameter groups multiple times.
"""

from __future__ import annotations

from astropy.table import Table

sim_backup_server = "https://simbad.idia.ac.za/simbad/sim-tap/"
sim_server = "http://simbad.u-strasbg.fr:80/simbad/sim-tap"
vizier_server = "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"
vizier_backup_server = "https://vizier.cfa.harvard.edu/viz-bin/VizieR"


FILTER_ABBREV = ["i", "j", "u", "g", "k", "u_sdss"]
PROVIDERS = ["sim", "sdb", "wds", "exo", "life", "gaia"]


def _mag_parameter_names() -> list[str]:
    return [f"mag_{filter_name}" for filter_name in FILTER_ABBREV]


def _mag_columns_for_parameter(parameter: str) -> list[str]:
    return [
        f"{parameter}_value",
        f"{parameter}_err",
        f"{parameter}_qual",
        f"{parameter}_source_idref",
        f"{parameter}_ref",
        f"{parameter}_sys",
    ]


def _parameter_spec(
    primary: str | None,
    quality: str | None,
    columns: list[str],
    dtype: list[type],
) -> dict[str, str | None | list[str] | list[type]]:
    if len(columns) != len(dtype):
        raise ValueError(
            f"Parameter specification for {primary!r} has "
            f"{len(columns)} columns but {len(dtype)} dtypes."
        )

    return {
        "primary": primary,
        "quality": quality,
        "columns": columns,
        "dtype": dtype,
    }


PARAMETER_SPECS = {
    "coo": _parameter_spec(
        primary="coo_ra",
        quality="coo_qual",
        columns=[
            "coo_ra",
            "coo_dec",
            "coo_err_angle",
            "coo_err_maj",
            "coo_err_min",
            "coo_qual",
            "coo_source_idref",
            "coo_ref",
        ],
        dtype=[float, float, float, float, float, object, int, object],
    ),
    "coo_gal": _parameter_spec(
        primary="coo_gal_l",
        quality="coo_gal_qual",
        columns=[
            "coo_gal_l",
            "coo_gal_b",
            "coo_gal_err_angle",
            "coo_gal_err_maj",
            "coo_gal_err_min",
            "coo_gal_qual",
            "coo_gal_source_idref",
            "coo_gal_ref",
        ],
        dtype=[float, float, float, float, float, object, int, object],
    ),
    "plx": _parameter_spec(
        primary="plx_value",
        quality="plx_qual",
        columns=[
            "plx_value",
            "plx_err",
            "plx_qual",
            "plx_source_idref",
            "plx_ref",
        ],
        dtype=[float, float, object, int, object],
    ),
    "dist_st": _parameter_spec(
        primary="dist_st_value",
        quality="dist_st_qual",
        columns=[
            "dist_st_value",
            "dist_st_err",
            "dist_st_qual",
            "dist_st_source_idref",
            "dist_st_ref",
        ],
        dtype=[float, float, object, int, object],
    ),
    "sptype": _parameter_spec(
        primary="sptype_string",
        quality="sptype_qual",
        columns=[
            "sptype_string",
            "sptype_err",
            "sptype_qual",
            "sptype_source_idref",
            "sptype_ref",
        ],
        dtype=[object, float, object, int, object],
    ),
    "class": _parameter_spec(
        primary="class_temp",
        quality="class_qual",
        columns=[
            "class_temp",
            "class_temp_nr",
            "class_lum",
            "class_qual",
            "class_source_idref",
            "class_ref",
        ],
        dtype=[object, object, object, object, int, object],
    ),
    "teff_st": _parameter_spec(
        primary="teff_st_value",
        quality="teff_st_qual",
        columns=[
            "teff_st_value",
            "teff_st_err",
            "teff_st_qual",
            "teff_st_source_idref",
            "teff_st_ref",
        ],
        dtype=[float, float, object, int, object],
    ),
    "radius_st": _parameter_spec(
        primary="radius_st_value",
        quality="radius_st_qual",
        columns=[
            "radius_st_value",
            "radius_st_err",
            "radius_st_qual",
            "radius_st_source_idref",
            "radius_st_ref",
        ],
        dtype=[float, float, object, int, object],
    ),
    "mass_st": _parameter_spec(
        primary="mass_st_value",
        quality="mass_st_qual",
        columns=[
            "mass_st_value",
            "mass_st_err",
            "mass_st_qual",
            "mass_st_source_idref",
            "mass_st_ref",
        ],
        dtype=[float, float, object, int, object],
    ),
    "binary": _parameter_spec(
        primary="binary_flag",
        quality="binary_qual",
        columns=[
            "binary_flag",
            "binary_qual",
            "binary_source_idref",
            "binary_ref",
        ],
        dtype=[object, object, int, object],
    ),
    "sep_ang": _parameter_spec(
        primary="sep_ang_value",
        quality="sep_ang_qual",
        columns=[
            "sep_ang_value",
            "sep_ang_err",
            "sep_ang_obs_date",
            "sep_ang_qual",
            "sep_ang_source_idref",
            "sep_ang_ref",
        ],
        dtype=[float, float, int, object, int, object],
    ),
    "mass_pl": _parameter_spec(
        primary="mass_pl_value",
        quality="mass_pl_qual",
        columns=[
            "mass_pl_value",
            "mass_pl_err_max",
            "mass_pl_err_min",
            "mass_pl_qual",
            "mass_pl_source_idref",
            "mass_pl_ref",
        ],
        dtype=[float, float, float, object, int, object],
    ),
    "rad": _parameter_spec(
        primary="rad_value",
        quality="rad_qual",
        columns=[
            "rad_value",
            "rad_err",
            "rad_rel",
            "rad_qual",
            "rad_source_idref",
            "rad_ref",
        ],
        dtype=[float, float, object, object, int, object],
    ),
    "h_link": _parameter_spec(
        primary="membership",
        quality=None,
        columns=[
            "child_object_idref",
            "parent_object_idref",
            "h_link_source_idref",
            "h_link_ref",
            "membership",
        ],
        dtype=[int, int, int, object, int],
    ),
}


for mag in _mag_parameter_names():
    PARAMETER_SPECS[mag] = _parameter_spec(
        primary=f"{mag}_value",
        quality=f"{mag}_qual",
        columns=_mag_columns_for_parameter(mag),
        dtype=[float, float, object, int, object, object],
    )


def _columns_for_parameters(parameters: list[str]) -> list[str]:
    columns = []
    for parameter in parameters:
        columns.extend(PARAMETER_SPECS[parameter]["columns"])
    return columns


def _dtypes_for_parameters(parameters: list[str]) -> list[type]:
    dtypes = []
    for parameter in parameters:
        dtypes.extend(PARAMETER_SPECS[parameter]["dtype"])
    return dtypes


STAR_BASIC_PARAMETERS = [
    "coo",
    "coo_gal",
    *_mag_parameter_names(),
    "plx",
    "dist_st",
    "sptype",
    "class",
    "teff_st",
    "radius_st",
    "mass_st",
    "binary",
    "sep_ang",
]

STAR_BASIC_BEST_PARAMETER_GROUPS = [
    "coo",
    "plx",
    "dist_st",
    "coo_gal",
    *_mag_parameter_names(),
    "class",
    "sptype",
]


TABLE_SCHEMAS = {
    "sources": {
        "parameters": [],
        "names": ["ref", "provider_name", "source_id"],
        "dtype": [object, object, int],
    },
    "objects": {
        "parameters": [],
        "names": ["type", "ids", "main_id"],
        "dtype": [object, object, object],
    },
    "provider": {
        "parameters": [],
        "names": [
            "provider_name",
            "provider_url",
            "provider_bibcode",
            "provider_access",
        ],
        "dtype": [object, object, object, object],
    },
    "ident": {
        "parameters": ["id"],
        "names": ["object_idref", "id", "id_source_idref"],
        "dtype": [int, object, int],
    },
    "h_link": {
        "parameters": ["h_link"],
        "names": PARAMETER_SPECS["h_link"]["columns"],
        "dtype": PARAMETER_SPECS["h_link"]["dtype"],
    },
    "star_basic": {
        "parameters": STAR_BASIC_BEST_PARAMETER_GROUPS,
        "names": ["object_idref", *_columns_for_parameters(STAR_BASIC_PARAMETERS)],
        "dtype": [int, *_dtypes_for_parameters(STAR_BASIC_PARAMETERS)],
    },
    "planet_basic": {
        "parameters": ["mass_pl"],
        "names": ["object_idref", *PARAMETER_SPECS["mass_pl"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["mass_pl"]["dtype"]],
    },
    "disk_basic": {
        "parameters": ["rad"],
        "names": ["object_idref", *PARAMETER_SPECS["rad"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["rad"]["dtype"]],
    },
    "mes_mass_pl": {
        "parameters": ["mass_pl"],
        "names": ["object_idref", *PARAMETER_SPECS["mass_pl"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["mass_pl"]["dtype"]],
    },
    "mes_teff_st": {
        "parameters": ["teff_st"],
        "names": ["object_idref", *PARAMETER_SPECS["teff_st"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["teff_st"]["dtype"]],
    },
    "mes_radius_st": {
        "parameters": ["radius_st"],
        "names": ["object_idref", *PARAMETER_SPECS["radius_st"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["radius_st"]["dtype"]],
    },
    "mes_mass_st": {
        "parameters": ["mass_st"],
        "names": ["object_idref", *PARAMETER_SPECS["mass_st"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["mass_st"]["dtype"]],
    },
    "mes_binary": {
        "parameters": ["binary"],
        "names": ["object_idref", *PARAMETER_SPECS["binary"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["binary"]["dtype"]],
    },
    "mes_sep_ang": {
        "parameters": ["sep_ang"],
        "names": ["object_idref", *PARAMETER_SPECS["sep_ang"]["columns"]],
        "dtype": [int, *PARAMETER_SPECS["sep_ang"]["dtype"]],
    },
    "best_h_link": {
        "parameters": ["h_link"],
        "names": PARAMETER_SPECS["h_link"]["columns"],
        "dtype": PARAMETER_SPECS["h_link"]["dtype"],
    },
}


def _empty_table() -> Table:
    return Table()


def _table_from_schema(table_name: str) -> Table:
    schema = TABLE_SCHEMAS[table_name]
    return Table(names=schema["names"], dtype=schema["dtype"])


def _create_empty_dict() -> dict[str, Table]:
    return {table_name: _empty_table() for table_name in table_names}


def _create_empty_dict_with_columns() -> dict[str, Table]:
    return {table_name: _table_from_schema(table_name) for table_name in table_names}


def _create_provider_tables_dict() -> dict[str, dict[str, Table]]:
    return {provider: _create_empty_dict() for provider in PROVIDERS}


table_names = list(TABLE_SCHEMAS)

filter_abbrev = FILTER_ABBREV

mags = _mag_parameter_names()

mag_columns = [
    [f"{mag}_value" for mag in mags],
    [f"{mag}_ref" for mag in mags],
    [f"{mag}_err" for mag in mags],
    [f"{mag}_qual" for mag in mags],
    [f"{mag}_sys" for mag in mags],
]

paras_dict = {
    table_name: schema["parameters"]
    for table_name, schema in TABLE_SCHEMAS.items()
}

para_columns_dict = {
    "star_basic": {
        parameter: {
            "primary": PARAMETER_SPECS[parameter]["primary"],
            "quality": PARAMETER_SPECS[parameter]["quality"],
            "columns": PARAMETER_SPECS[parameter]["columns"],
        }
        for parameter in STAR_BASIC_BEST_PARAMETER_GROUPS
    },
}

empty_dict = _create_empty_dict()

empty_dict_wit_columns = _create_empty_dict_with_columns()

empty_provider_tables_dict = _create_provider_tables_dict()
