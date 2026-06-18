from utils.io import load, stringtoobject, save
from astropy.io.ascii import read
from utils.io import Path
from utils.fcts_cat_merge import nearest_neighbor_distances_units, get_mask_cat2_in_cat1
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.table import vstack, MaskedColumn
from provider.utils import nullvalues
from utils.analysis import catalog_versions



def model_exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

def get_catalog(name):
    if name == "hpic":
        catalog = read(Path().additional_data + "/HPICv1.0/full_HPIC.txt",
                    delimiter="|")

    if name == "starcat5":
        [catalog] = load(["catalogs/StarCat5"])

    return catalog

def get_radius(catalog_ra,catalog_dec):
    radius = nearest_neighbor_distances_units(catalog_ra,
                                              catalog_dec)

    plt.figure()

    bin_heights, bin_borders, _ = plt.hist(
        radius,
        bins=100,
        alpha=0.5
    )

    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    popt, pcov = curve_fit(model_exp_decay,
                        bin_centers,
                        bin_heights,
                        p0=[100000, 0.05, 300])
    a_opt, b_opt, c_opt = popt

    #print(a_opt, b_opt, c_opt)
    #print(np.linalg.cond(pcov)) # is a bit high, could be overparametrized

    # Create fitted curve
    x_model = np.linspace(bin_centers[0], max(bin_borders), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt)

    # Plot fitted curve
    plt.plot(x_model, y_model, label=f'fit')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("nearest-neighbor angular distances")
    plt.ylabel("number of objects")


    # get x of fit where y is 0.1*ymax
    # ymax = y (x=0)
    ymax = model_exp_decay(bin_centers[0], a_opt, b_opt, c_opt)
    #print("ymax: ",ymax)
    x_temp = x_model[np.where(y_model < 0.9*ymax)]
    radius = min(x_temp)
    print("radius: ",radius)

    plt.plot([radius],[y_model[np.where(x_model == radius)]],'o',color='red',label='radius')

    plt.show()
    return radius

def rename_cols(catalog,
                  colnames,
                  new_colnames):
    pre_merge_cat = catalog.copy()
    pre_merge_cat = stringtoobject(pre_merge_cat)

    pre_merge_cat.rename_columns(colnames, new_colnames)

    return pre_merge_cat

def deal_with_nulls(catalog,
                    null_columns,
                    null):
    for col in null_columns:
        mask = np.where(catalog[col] == null, True, False)
        catalog[col] = MaskedColumn(catalog[col], mask=mask)
    return catalog

def deal_with_resto_of_hpic_cols(catalog,binary_flag_col,float_cols):
    # change binary column (0-False, 1-True)
    catalog[binary_flag_col] = catalog[
        binary_flag_col].astype(bool)
    catalog[binary_flag_col] = catalog[
        binary_flag_col].astype(object)

    for col in float_cols:
        catalog = nullvalues(catalog, col, np.nan)
        catalog[col] = catalog[col].astype(float)
    # null value below mask can't get transformed into float... do I use np.nan?

    return catalog

def analysis(pre_merge_hpic_masked,starcat5,catalog,float_colnames,):
    ylabel = "Stellar Temperature [K]"

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )  # subplots so that I can overplot old version?

    ax.scatter(starcat5["dist_st_value"], starcat5["teff_st_value"], s=2,
               alpha=0.5)
    ax.scatter(pre_merge_hpic_masked["temp_dist_st_value"],
               pre_merge_hpic_masked["temp_teff_st_value"], s=2, alpha=0.5)
    ax.set_yscale("log")

    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel(ylabel)
    plt.show()

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )  # subplots so that I can overplot old version?

    ax.scatter(catalog["temp_dist_st_value"],
               catalog["temp_teff_st_value"], s=2, alpha=0.5)
    ax.set_yscale("log")

    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel(ylabel)
    plt.show()

    for col in float_colnames:
        x = pre_merge_hpic_masked['temp_' + col]
        hpic_prepped = x[~np.isnan(x)]
        y = catalog['temp_' + col]
        catalog_prepped = y[~np.isnan(y)]
        catalog_versions.threecatboxplot([starcat5[col],
                                          hpic_prepped, catalog_prepped],
                                         col, ["starcat5", "hpic", "merger"])
    return

def hpic_merger():
    hpic = get_catalog("hpic")
    starcat5 = get_catalog("starcat5")
    radius = get_radius(starcat5["coo_ra"], starcat5["coo_dec"])

    mask_cat2_in_cat1 = get_mask_cat2_in_cat1(name_cat1=hpic["simbad_name"],
                                              ra_cat1=hpic["ra"],
                                              dec_cat1=hpic["dec"],
                                              name_cat2=starcat5["main_id"],
                                              ra_cat2=starcat5["coo_ra"],
                                              dec_cat2=starcat5["coo_dec"],
                                              r_arcsec=radius)

    masked_starcat5 = starcat5[mask_cat2_in_cat1]

    starcat5_merge_colnames = ["main_id", "coo_ra", "coo_dec", "sptype_string",
                               "plx_value", "dist_st_value", "teff_st_value",
                               "teff_ref", "mass_st_value", "mass_ref",
                               "radius_st_value", "radius_ref", "mag_i_value",
                               "mag_j_value", "binary_flag", "sep_ang_value"]

    hpic_merge_colnames = ["star_name", "ra", "dec", "st_spectype", "sy_plx",
                           "sy_dist", "st_teff", "st_teff_reflink", "st_mass",
                           "st_mass_reflink", "st_rad", "st_rad_reflink",
                           "sy_icmag", "sy_jmag", "known_binary_fl", "wds_sep"]

    new_colnames = ["temp_" + col for col in starcat5_merge_colnames]

    starcat5_null_columns = ["temp_sptype_string", "temp_radius_ref",
                             "temp_teff_ref", "temp_mass_ref"]
    starcat5_null = ""

    hpic_null_columns = ["temp_sptype_string", "temp_mass_st_value",
                         "temp_radius_st_value", "temp_radius_ref",
                         "temp_mag_i_value", "temp_mag_j_value",
                         "temp_sep_ang_value", "temp_mass_ref",
                         "temp_plx_value", "temp_dist_st_value",
                         "temp_teff_st_value"]
    hpic_null = "null"

    renamed_cols_starcat = rename_cols(masked_starcat5,
                                       starcat5_merge_colnames,
                                       new_colnames)

    renamed_cols_hpic = rename_cols(hpic,
                                    hpic_merge_colnames,
                                    new_colnames)

    pre_merge_starcat = deal_with_nulls(renamed_cols_starcat,
                                        starcat5_null_columns,
                                        starcat5_null)

    pre_merge_hpic = deal_with_nulls(renamed_cols_hpic,
                                     hpic_null_columns,
                                     hpic_null)

    # float columns null values
    float_colnames = ["plx_value", "mag_i_value", "mag_j_value",
                      "dist_st_value",
                      "teff_st_value", "radius_st_value", "mass_st_value",
                      "sep_ang_value"]
    temp_float_colnames = ['temp_' + float_colnames[j] for j in
                           range(len(float_colnames))]

    pre_merge_hpic = deal_with_resto_of_hpic_cols(pre_merge_hpic,
                                                  "temp_binary_flag",
                                                  temp_float_colnames)

    catalog = vstack([pre_merge_hpic, pre_merge_starcat])
    save([catalog], ["HPIC_StarCat"])

    analysis(pre_merge_hpic, starcat5,catalog,float_colnames)

    return catalog








