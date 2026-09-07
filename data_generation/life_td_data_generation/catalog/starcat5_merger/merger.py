from collections.abc import Sequence
from typing import Any, Literal

import importlib

import matplotlib.pyplot as plt
import numpy as np
from astropy.io.ascii import read
from astropy.table import MaskedColumn, Table, vstack
from catalog.starcat5_merger.fcts_cat_merge import (
    get_cat2_in_cat1_match_info,
    nearest_neighbor_distances_units,
    add_cat2_to_cat1_by_name_or_coords,
)
from provider.utils import nullvalues
from scipy.optimize import curve_fit
from utils.analysis import catalog_versions, finalplot
from utils.io import load, save, stringtoobject

importlib.reload(catalog_versions)


ADDITIONAL_DATA_LOCATION = "../../../../additional_data/"
HPIC_FILE = "../../../../additional_data/HPICv1.0/full_HPIC.txt"
STARCAT5_PATH = "catalogs/StarCat5"


def model_exp_decay(
    x: np.ndarray | float,
    a: float,
    b: float,
    c: float,
) -> np.ndarray | float:
    """
    Evaluate an exponential decay model.

    :param x: Input value or array of values.
    :type x: numpy.ndarray or float
    :param a: Amplitude of the exponential term.
    :type a: float
    :param b: Decay coefficient.
    :type b: float
    :param c: Constant offset.
    :type c: float
    :returns: Model value ``a * exp(-b * x) + c``.
    :rtype: numpy.ndarray or float
    """
    return a * np.exp(-b * x) + c


def get_catalog(name: Literal["hpic", "starcat5"]) -> Table:
    """
    Load one of the catalogs used for the HPIC-StarCat5 merger.

    :param name: Catalog selector. Supported values are ``"hpic"`` and
        ``"starcat5"``.
    :type name: str
    :returns: Loaded catalog table.
    :rtype: astropy.table.Table
    """
    if name == "hpic":
        catalog = read(HPIC_FILE, delimiter="|")

    if name == "starcat5":
        [catalog] = load(
            [STARCAT5_PATH],
            location=ADDITIONAL_DATA_LOCATION,
        )

    return catalog


def get_radius(catalog_ra: np.ndarray, catalog_dec: np.ndarray) -> float:
    """
    Estimate a positional matching radius from nearest-neighbor distances.

    Nearest-neighbor angular distances are histogrammed and fitted with an
    exponential decay. The matching radius is chosen as the first model-grid
    position where the fitted curve falls below 90 percent of the fitted value
    at the first histogram-bin center.

    :param catalog_ra: Right ascensions of catalog objects.
    :type catalog_ra: numpy.ndarray
    :param catalog_dec: Declinations of catalog objects.
    :type catalog_dec: numpy.ndarray
    :returns: Estimated matching radius in arcseconds.
    :rtype: float
    """
    nearest_neighbor_radius = nearest_neighbor_distances_units(
        catalog_ra,
        catalog_dec,
    )

    plt.figure()

    bin_heights, bin_borders, _ = plt.hist(
        nearest_neighbor_radius,
        bins=100,
        alpha=0.5,
    )
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    popt, _ = curve_fit(
        model_exp_decay,
        bin_centers,
        bin_heights,
        p0=[100000, 0.05, 300],
    )
    a_opt, b_opt, c_opt = popt

    #print(a_opt, b_opt, c_opt)
    #print(np.linalg.cond(pcov)) # is a bit high, could be overparametrized

    # Create fitted curve
    x_model = np.linspace(bin_centers[0], max(bin_borders), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt)

    plt.plot(x_model, y_model, label="fit")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("nearest-neighbor angular distances")
    plt.ylabel("number of objects")

    ymax = model_exp_decay(bin_centers[0], a_opt, b_opt, c_opt)
    x_temp = x_model[np.where(y_model < 0.9 * ymax)]
    radius = min(x_temp)

    print("radius: ", radius)

    plt.plot(
        [radius],
        [y_model[np.where(x_model == radius)]],
        "o",
        color="red",
        label="radius",
    )
    plt.savefig("../../../../plots/radius_fit.png")
    plt.show()


    return radius


def rename_cols(
    catalog: Table,
    colnames: Sequence[str],
    new_colnames: Sequence[str],
) -> Table:
    """
    Copy a catalog, convert string columns to objects, and rename columns.

    :param catalog: Input catalog table.
    :type catalog: astropy.table.Table
    :param colnames: Existing column names.
    :type colnames: sequence[str]
    :param new_colnames: Replacement column names.
    :type new_colnames: sequence[str]
    :returns: Prepared copy with renamed columns.
    :rtype: astropy.table.Table
    """
    pre_merge_cat = catalog.copy()
    pre_merge_cat = stringtoobject(pre_merge_cat)
    pre_merge_cat.rename_columns(colnames, new_colnames)

    return pre_merge_cat


def deal_with_nulls(
    catalog: Table,
    null_columns: Sequence[str],
    null: Any,
) -> Table:
    """
    Mask selected columns where they contain a provider-specific null value.

    :param catalog: Catalog table to update in place.
    :type catalog: astropy.table.Table
    :param null_columns: Columns in which ``null`` should be masked.
    :type null_columns: sequence[str]
    :param null: Provider-specific null sentinel value.
    :type null: object
    :returns: Catalog with selected null values masked.
    :rtype: astropy.table.Table
    """
    for col in null_columns:
        mask = np.where(catalog[col] == null, True, False)
        catalog[col] = MaskedColumn(catalog[col], mask=mask)

    return catalog


def deal_with_resto_of_hpic_cols(
    catalog: Table,
    binary_flag_col: str,
    float_cols: Sequence[str],
) -> Table:
    """
    Normalize remaining HPIC columns before stacking with StarCat5 rows.

    The binary flag is converted from HPIC's numeric representation to bool and
    then to object dtype. Selected numeric columns have ``np.nan`` normalized via
    :func:`provider.utils.nullvalues` and are converted to float.

    :param catalog: HPIC catalog table to update.
    :type catalog: astropy.table.Table
    :param binary_flag_col: Name of the binary-flag column.
    :type binary_flag_col: str
    :param float_cols: Numeric columns to normalize and cast to float.
    :type float_cols: sequence[str]
    :returns: Updated HPIC catalog table.
    :rtype: astropy.table.Table
    """
    catalog[binary_flag_col] = catalog[binary_flag_col].astype(bool)
    catalog[binary_flag_col] = catalog[binary_flag_col].astype(object)

    for col in float_cols:
        catalog = nullvalues(catalog, col, np.nan)
        catalog[col] = catalog[col].astype(float)

    return catalog


def scatter_plot(
    catalogs: Sequence[Table],
    x_values: Sequence[str],
    y_values: Sequence[str],
    ylabel: str,
) -> None:
    """
    Create a scatter plot for one or more catalogs.

    :param catalogs: Catalog tables to plot.
    :type catalogs: sequence[astropy.table.Table]
    :param x_values: Column names to use as x-values for each catalog.
    :type x_values: sequence[str]
    :param y_values: Column names to use as y-values for each catalog.
    :type y_values: sequence[str]
    :param ylabel: Label for the y-axis.
    :type ylabel: str
    :returns: None.
    :rtype: None
    """
    _, ax = plt.subplots(figsize=(9, 6))

    for cat, x_col, y_col in zip(catalogs, x_values, y_values):
        ax.scatter(cat[x_col], cat[y_col], s=2, alpha=0.5)

    ax.set_yscale("log")
    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel(ylabel)

    plt.show()


def merger_analysis(
    pre_merge_hpic_masked: Table,
    starcat5_not_in_hpic: Table,
    catalog: Table,
    float_colnames: Sequence[str],
) -> None:
    """
    Create diagnostic plots for the merged HPIC-StarCat5 catalog.

    :param pre_merge_hpic_masked: Prepared HPIC table before stacking.
    :type pre_merge_hpic_masked: astropy.table.Table
    :param starcat5_not_in_hpic: StarCat5 rows not matched to HPIC.
    :type starcat5_not_in_hpic: astropy.table.Table
    :param catalog: Stacked merger catalog.
    :type catalog: astropy.table.Table
    :param float_colnames: Base names of numeric columns to compare.
    :type float_colnames: sequence[str]
    :returns: None.
    :rtype: None
    """
    scatter_plot(
        [starcat5_not_in_hpic, pre_merge_hpic_masked],
        ["temp_dist_st_value", "temp_dist_st_value"],
        ["temp_teff_st_value", "temp_teff_st_value"],
        "Stellar Temperature [K]",
    )

    catalog_versions.plot_cat_paras(
        ["temp_teff_st_value", "temp_radius_st_value"],
        [starcat5_not_in_hpic, pre_merge_hpic_masked],
        label_list=["starcat5_not_in_hpic", "pre_merge_hpic_masked"],
    )

    scatter_plot(
        [catalog],
        ["temp_dist_st_value"],
        ["temp_teff_st_value"],
        "Stellar Temperature [K]",
    )

    for col in float_colnames:
        hpic_values = pre_merge_hpic_masked["temp_" + col]
        hpic_prepped = hpic_values[~np.isnan(hpic_values)]

        catalog_values = catalog["temp_" + col]
        catalog_prepped = catalog_values[~np.isnan(catalog_values)]

        catalog_versions.threecatboxplot(
            [starcat5_not_in_hpic[col], hpic_prepped, catalog_prepped],
            col,
            ["starcat5_not_in_hpic", "hpic", "merger"],
        )


def plot_para_vs_para(hpic: Table, starcat5_wo_hpic: Table) -> None:
    """
    Plot selected parameter pairs for HPIC and StarCat5-only additions.

    :param hpic: Prepared HPIC table.
    :type hpic: astropy.table.Table
    :param starcat5_wo_hpic: Prepared StarCat5 rows not present in HPIC.
    :type starcat5_wo_hpic: astropy.table.Table
    :returns: None.
    :rtype: None
    """
    labels = ["HPIC", "StarCat5_addition_to_HPIC"]

    catalog_versions.plot_cat_paras(
        ["temp_teff_st_value", "temp_radius_st_value"],
        [hpic, starcat5_wo_hpic],
        label_list=labels,
        min = [2000,0],
        max = [10000,2.5],
    )
    catalog_versions.plot_cat_paras(
        ["temp_teff_st_value", "temp_mass_st_value"],
        [hpic, starcat5_wo_hpic],
        label_list=labels,
        min=[2000, 0],
        max=[10000, 2],
    )
    catalog_versions.plot_cat_paras(
        ["temp_radius_st_value", "temp_mass_st_value"],
        [hpic, starcat5_wo_hpic],
        label_list=labels,
        min=[0, 0],
        max=[2.5, 2],
    )
    catalog_versions.plot_cat_paras(
        ["temp_coo_ra", "temp_coo_dec"],
        [hpic, starcat5_wo_hpic],
        label_list=labels,
    )


def analysis_starcat5_not_in_hpic(catalog: Table) -> None:
    """
    Create diagnostic plots for StarCat5 rows not matched to HPIC.

    :param catalog: StarCat5-only catalog subset.
    :type catalog: astropy.table.Table
    :returns: None.
    :rtype: None
    """
    scatter_plot(
        [catalog],
        ["temp_dist_st_value"],
        ["temp_teff_st_value"],
        "Stellar Temperature [K]",
    )

    scatter_plot(
        [catalog],
        ["temp_dist_st_value"],
        ["temp_mag_j_value"],
        "J Magnitude",
    )

    finalplot.starcat_distribution_plot(
        [catalog["class_temp", "temp_dist_st_value"]],
        ["StarCat5_addition_to_HPIC"],
    )

    for temp_class in ["O", "B", "A", "F", "G"]:
        print(
            catalog["temp_main_id", "class_temp"][
                np.where(catalog["class_temp"] == temp_class)
            ]
        )



def get_merge_column_names() -> tuple[list[str], list[str], list[str]]:
    """
    Define matching StarCat5, HPIC, and temporary merger column names.

    :returns: Tuple with StarCat5 source columns, HPIC source columns, and the
        corresponding temporary target names.
    :rtype: (list[str], list[str], list[str])
    """
    starcat5_merge_colnames = [
        "main_id",
        "coo_ra",
        "coo_dec",
        "sptype_string",
        "plx_value",
        "dist_st_value",
        "teff_st_value",
        "teff_ref",
        "mass_st_value",
        "mass_ref",
        "radius_st_value",
        "radius_ref",
        "mag_i_value",
        "mag_i_err",
        "mag_i_ref",
        "mag_j_value",
        "mag_j_err",
        "mag_j_ref",
        "mag_u_value",
        "mag_u_err",
        "mag_u_ref",
        "mag_g_value",
        "mag_g_err",
        "mag_g_ref",
        "binary_flag",
        "sep_ang_value",
    ]

    hpic_merge_colnames = [
        "star_name",
        "ra",
        "dec",
        "st_spectype",
        "sy_plx",
        "sy_dist",
        "st_teff",
        "st_teff_reflink",
        "st_mass",
        "st_mass_reflink",
        "st_rad",
        "st_rad_reflink",
        "sy_icmag",
        "sy_icmagerr",
        "sy_icmag_reflink",
        "sy_jmag",
        "sy_jmagerr",
        "sy_jmag_reflink",
        "sy_ujmag",
        "sy_ujmagerr",
        "sy_ujmag_reflink",
        "sy_gaiamag",
        "sy_gaiamagerr",
        "sy_gaiamag_reflink",
        "known_binary_fl",
        "wds_sep",
    ]

    new_colnames = [f"temp_{col}" for col in starcat5_merge_colnames]

    return starcat5_merge_colnames, hpic_merge_colnames, new_colnames


def get_null_column_settings() -> tuple[list[str], str, list[str], str]:
    """
    Define provider-specific null columns and sentinel values.

    :returns: StarCat5 null columns, StarCat5 null value, HPIC null columns, and
        HPIC null value.
    :rtype: (list[str], str, list[str], str)
    """
    starcat5_null_columns = [
        "temp_sptype_string",
        "temp_radius_ref",
        "temp_teff_ref",
        "temp_mass_ref",
        "temp_mag_i_ref",
        "temp_mag_j_ref",
        "temp_mag_u_ref",
        "temp_mag_g_ref",
    ]

    hpic_null_columns = [
        "temp_sptype_string",
        "temp_mass_st_value",
        "temp_mass_ref",
        "temp_radius_st_value",
        "temp_radius_ref",
        "temp_teff_st_value",
        "temp_teff_ref",
        "temp_mag_i_value",
        "temp_mag_j_value",
        "temp_mag_u_value",
        "temp_mag_g_value",
        "temp_mag_i_err",
        "temp_mag_j_err",
        "temp_mag_u_err",
        "temp_mag_g_err",
        "temp_mag_i_ref",
        "temp_mag_j_ref",
        "temp_mag_u_ref",
        "temp_mag_g_ref",
        "temp_sep_ang_value",
        "temp_plx_value",
        "temp_dist_st_value",

    ]

    return starcat5_null_columns, "", hpic_null_columns, "null"


def get_float_colnames() -> list[str]:
    """
    Define numeric merger columns that need float conversion for HPIC.

    :returns: Base names of numeric columns.
    :rtype: list[str]
    """
    return [
        "plx_value",
        "dist_st_value",
        "mag_i_value",
        "mag_j_value",
        "mag_u_value",
        "mag_g_value",
        "mag_i_err",
        "mag_j_err",
        "mag_u_err",
        "mag_g_err",
        "teff_st_value",
        "radius_st_value",
        "mass_st_value",
        "sep_ang_value",
    ]


def prepare_pre_merge_catalogs(
    hpic: Table,
    starcat5_not_in_hpic: Table,
) -> tuple[Table, Table, list[str]]:
    """
    Rename, mask, and cast HPIC and StarCat5 tables before stacking.

    :param hpic: Original HPIC catalog.
    :type hpic: astropy.table.Table
    :param starcat5_not_in_hpic: StarCat5 subset not matched to HPIC.
    :type starcat5_not_in_hpic: astropy.table.Table
    :returns: Prepared StarCat5 table, prepared HPIC table, and float-column
        base names.
    :rtype: (astropy.table.Table, astropy.table.Table, list[str])
    """
    (
        starcat5_merge_colnames,
        hpic_merge_colnames,
        new_colnames,
    ) = get_merge_column_names()

    (
        starcat5_null_columns,
        starcat5_null,
        hpic_null_columns,
        hpic_null,
    ) = get_null_column_settings()

    renamed_cols_starcat = rename_cols(
        starcat5_not_in_hpic,
        starcat5_merge_colnames,
        new_colnames,
    )
    renamed_cols_hpic = rename_cols(
        hpic,
        hpic_merge_colnames,
        new_colnames,
    )

    pre_merge_starcat = deal_with_nulls(
        renamed_cols_starcat,
        starcat5_null_columns,
        starcat5_null,
    )
    pre_merge_hpic = deal_with_nulls(
        renamed_cols_hpic,
        hpic_null_columns,
        hpic_null,
    )

    float_colnames = get_float_colnames()
    temp_float_colnames = [f"temp_{col}" for col in float_colnames]

    pre_merge_hpic = deal_with_resto_of_hpic_cols(
        pre_merge_hpic,
        "temp_binary_flag",
        temp_float_colnames,
    )

    return pre_merge_starcat, pre_merge_hpic, float_colnames


def hpic_merger() -> tuple[Table, Table, Table, Table, list[str]]:
    """
    Merge HPIC with StarCat5 rows that are not already present in HPIC.

    The routine loads HPIC and StarCat5, estimates a StarCat5 positional
    matching radius, flags which rows are present in the other catalog, prepares
    both catalogs to shared temporary column names, stacks them, and saves the
    resulting intermediate products.

    :returns: Tuple containing the merged catalog, prepared StarCat5-only table,
        original StarCat5 table with match flag, prepared HPIC table with match
        flag, and float-column base names.
    :rtype: (
        astropy.table.Table,
        astropy.table.Table,
        astropy.table.Table,
        astropy.table.Table,
        list[str]
    )
    """
    hpic = get_catalog("hpic")
    starcat5 = get_catalog("starcat5")

    radius = get_radius(starcat5["coo_ra"], starcat5["coo_dec"])

    match_info = get_cat2_in_cat1_match_info(
        name_cat1=hpic["simbad_name"],
        ra_cat1=hpic["ra"],
        dec_cat1=hpic["dec"],
        name_cat2=starcat5["main_id"],
        ra_cat2=starcat5["coo_ra"],
        dec_cat2=starcat5["coo_dec"],
        r_arcsec=radius,
    )

    mask_cat2_in_cat1 = match_info["mask_cat2_in_cat1"]

    starcat5["in_hpic"] = mask_cat2_in_cat1

    hpic_has_starcat5_match = np.zeros(len(hpic), dtype=bool)
    matched_hpic_indices = match_info["cat1_index_for_cat2"][mask_cat2_in_cat1]
    hpic_has_starcat5_match[matched_hpic_indices] = True

    starcat5_not_in_hpic = starcat5[np.invert(starcat5["in_hpic"])]

    # adding mag_u_columns to hpic
    hpic_mag_u_joined = add_cat2_to_cat1_by_name_or_coords(
        cat1=hpic,
        cat2=starcat5,
        name_cat1_col="simbad_name",
        ra_cat1_col="ra",
        dec_cat1_col="dec",
        name_cat2_col="main_id",
        ra_cat2_col="coo_ra",
        dec_cat2_col="coo_dec",
        r_arcsec=radius,
        cat2_cols=[
            "main_id",
            "coo_ra",
            "coo_dec",
            "mag_u_sdss_value",
            "mag_u_sdss_err",
            "mag_u_sdss_ref",
            "mag_u_sdss_sys",
        ],
        match_info=match_info,
    )

    hpic_mag_u_joined["in_starcat5"] = hpic_has_starcat5_match

    pre_merge_starcat, pre_merge_hpic, float_colnames = (
        prepare_pre_merge_catalogs(hpic_mag_u_joined, starcat5_not_in_hpic)
    )

    pre_merge_starcat["in_hpic"] = False
    pre_merge_starcat["in_starcat5"] = True
    pre_merge_hpic["in_hpic"] = True

    HPIC_StarCat = vstack([pre_merge_hpic, pre_merge_starcat])

    save(
        [
            HPIC_StarCat,
            starcat5,
            pre_merge_hpic,
            pre_merge_starcat,
        ],
        [
            "HPIC_StarCat",
            "StarCat5_with_HPIC_flag",
            "pre_merge_hpic",
            "pre_merge_starcat",
        ],
        location=ADDITIONAL_DATA_LOCATION,
    )

    # merger_analysis(pre_merge_hpic, starcat5, catalog, float_colnames)

    return (
        HPIC_StarCat,
        pre_merge_starcat,
        starcat5,
        pre_merge_hpic,
        float_colnames,
    )
