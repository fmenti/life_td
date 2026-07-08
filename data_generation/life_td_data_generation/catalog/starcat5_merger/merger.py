from utils.io import load, stringtoobject, save
from astropy.io.ascii import read
from utils.io import Path
from catalog.starcat5_merger.fcts_cat_merge import nearest_neighbor_distances_units, get_mask_cat2_in_cat1
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.table import vstack, MaskedColumn
from provider.utils import nullvalues
from utils.analysis import catalog_versions, finalplot
import importlib
importlib.reload(catalog_versions)



def model_exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

def get_catalog(name):
    if name == "hpic":
        catalog = read("../../../../additional_data/HPICv1.0/full_HPIC.txt",
                    delimiter="|")

    if name == "starcat5":
        [catalog] = load(["catalogs/StarCat5"],
                         location = "../../../../additional_data/")

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

def scatter_plot(catalogs,x_values,y_values,ylabel):

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )  # subplots so that I can overplot old version?

    for cat,x,y in zip(catalogs,x_values,y_values):
        ax.scatter(cat[x], cat[y], s=2,
               alpha=0.5)

    ax.set_yscale("log")

    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel(ylabel)
    plt.show()

def merger_analysis(pre_merge_hpic_masked,starcat5_not_in_hpic,catalog,float_colnames,):
    scatter_plot([starcat5_not_in_hpic,pre_merge_hpic_masked],
                 ["temp_dist_st_value","temp_dist_st_value"],
                 ["temp_teff_st_value","temp_teff_st_value"],
                 "Stellar Temperature [K]")

    catalog_versions.plot_cat_paras(
        ["temp_teff_st_value", "temp_radius_st_value"],
        [starcat5_not_in_hpic, pre_merge_hpic_masked],
        label_list=["starcat5_not_in_hpic", "pre_merge_hpic_masked"])

    scatter_plot([catalog],
                 ["temp_dist_st_value"],
                 [ "temp_teff_st_value"],
                 "Stellar Temperature [K]")

    for col in float_colnames:
        x = pre_merge_hpic_masked['temp_' + col]
        hpic_prepped = x[~np.isnan(x)]
        y = catalog['temp_' + col]
        catalog_prepped = y[~np.isnan(y)]
        catalog_versions.threecatboxplot([starcat5_not_in_hpic[col],
                                          hpic_prepped, catalog_prepped],
                                         col, ["starcat5_not_in_hpic", "hpic", "merger"])
    return


def plot_para_vs_para(hpic,starcat5_wo_hpic):

    catalog_versions.plot_cat_paras(["temp_teff_st_value","temp_radius_st_value"],
                   [hpic,starcat5_wo_hpic],
                   label_list = ["HPIC","StarCat5_addition_to_HPIC"])
    catalog_versions.plot_cat_paras(["temp_teff_st_value", "temp_mass_st_value"],
                   [hpic, starcat5_wo_hpic],
                   label_list = ["HPIC","StarCat5_addition_to_HPIC"])
    catalog_versions.plot_cat_paras(["temp_radius_st_value", "temp_mass_st_value"],
                   [hpic, starcat5_wo_hpic],
                   label_list = ["HPIC","StarCat5_addition_to_HPIC"])

    catalog_versions.plot_cat_paras(["temp_coo_ra", "temp_coo_dec"],
                   [hpic, starcat5_wo_hpic],
                   label_list = ["HPIC","StarCat5_addition_to_HPIC"])

    return

def analysis_starcat5_not_in_hpic(catalog):
    scatter_plot([catalog],
                 ["temp_dist_st_value"],
                 ["temp_teff_st_value"],
                 "Stellar Temperature [K]")

    scatter_plot([catalog],
                 ["temp_dist_st_value"],
                 ["temp_mag_j_value"],
                 "J Magnitude")

    finalplot.starcat_distribution_plot(
        [catalog["class_temp", "temp_dist_st_value"]], ["StarCat5_addition_to_HPIC"]
    )

    # print(catalog)

    for temp_class in ["O","B","A","F","G"]:
        print(catalog["temp_main_id","class_temp"][np.where(
            catalog["class_temp"] == temp_class)])

    return

def spec_dist_plot(catalogs,x):
    plt.figure()

    for cat,label in zip(catalogs,["HPIC","StarCat5_addition_to_HPIC"]):
        spectype = np.array(cat["temp_sptype_string"]).astype(str)

        # Reduce spectral types to first two characters, e.g. "M3.4" -> "M3"
        spectype = np.array([s[:2] for s in spectype])

        # Keep only spectral types listed in x
        spectype = spectype[np.isin(spectype, x)]


        plt.hist(spectype, bins=np.arange(len(x) + 1) - 0.5, edgecolor="black",label = label, alpha = 0.5)
    plt.xticks(range(len(x)), x)
    plt.xlabel("Spectral Subclass")
    plt.ylabel("Number of stars")
    plt.legend()

    plt.show()

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

    starcat5_in_hpic = starcat5[mask_cat2_in_cat1].copy()

    starcat5_not_in_hpic = starcat5[np.invert(mask_cat2_in_cat1)]

    starcat5_merge_colnames = ["main_id", "coo_ra", "coo_dec", "sptype_string",
                               "plx_value", "dist_st_value", "teff_st_value",
                               "teff_ref", "mass_st_value", "mass_ref",
                               "radius_st_value", "radius_ref", "mag_i_value",
                               "mag_j_value", "mag_u_value", "binary_flag",
                               "sep_ang_value"]

    hpic_merge_colnames = ["star_name", "ra", "dec", "st_spectype", "sy_plx",
                           "sy_dist", "st_teff", "st_teff_reflink", "st_mass",
                           "st_mass_reflink", "st_rad", "st_rad_reflink",
                           "sy_icmag", "sy_jmag", "sy_ujmag", "known_binary_fl",
                           "wds_sep"]

    new_colnames = ["temp_" + col for col in starcat5_merge_colnames]

    starcat5_null_columns = ["temp_sptype_string", "temp_radius_ref",
                             "temp_teff_ref", "temp_mass_ref"]
    starcat5_null = ""

    hpic_null_columns = ["temp_sptype_string", "temp_mass_st_value",
                         "temp_radius_st_value", "temp_radius_ref",
                         "temp_mag_i_value", "temp_mag_j_value",
                         "temp_mag_u_value",
                         "temp_sep_ang_value", "temp_mass_ref",
                         "temp_plx_value", "temp_dist_st_value",
                         "temp_teff_st_value"]
    hpic_null = "null"

    renamed_cols_starcat = rename_cols(starcat5_not_in_hpic,
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
    float_colnames = ["plx_value", "mag_i_value", "mag_j_value", "mag_u_value",
                      "dist_st_value",
                      "teff_st_value", "radius_st_value", "mass_st_value",
                      "sep_ang_value"]
    temp_float_colnames = ['temp_' + float_colnames[j] for j in
                           range(len(float_colnames))]

    pre_merge_hpic = deal_with_resto_of_hpic_cols(pre_merge_hpic,
                                                  "temp_binary_flag",
                                                  temp_float_colnames)

    catalog = vstack([pre_merge_hpic, pre_merge_starcat])
    print(catalog)

    save([catalog,starcat5_not_in_hpic,pre_merge_hpic,pre_merge_starcat],
         ["HPIC_StarCat","starcat5_not_in_hpic","pre_merge_hpic","pre_merge_starcat"],
         location="../../../../additional_data/"
         )

    #merger_analysis(pre_merge_hpic, starcat5,catalog,float_colnames)

    return catalog, pre_merge_starcat, starcat5, pre_merge_hpic, float_colnames, starcat5_in_hpic








