import importlib  # reloading external functions after modification
import sys
from typing import Any

import astropy as ap  # votables
import numpy as np  # arrays
import seaborn as sns
from matplotlib import pyplot as plt

sys.path.append(
    "../../implementation/life/data_generation/life_td_data_generation"
)

# self created modules
import utils.analysis.analysis as la
from utils.analysis.catalog_comparison import create_common
from utils.io import Path, load, save

importlib.reload(la)  # reload module after changing it
import provider as p

importlib.reload(p)  # reload module after changing it

para3 = [
    "sim_ra",
    "sim_dec",
    "sim_plx",
    "distance",
    "gal_coord_l",
    "gal_coord_b",
    "mod_Teff",
    "mod_R",
    "mod_M",
    "sep_phys",
    "sim_i",
    "sim_j",
]
para4 = [
    "coo_ra",
    "coo_dec",
    "plx_value",
    "dist_st_value",
    "coo_gal_l",
    "coo_gal_b",
    "teff_st_value",
    "radius_st_value",
    "mass_st_value",
    "sep_phys_value",
    "mag_i_value",
    "mag_j_value",
]

para_hwo = ["ra", "dec", "sy_plx", "sy_dist", "st_teff", "st_rad", "st_mass"]
para_starcat4 = [
    "coo_ra",
    "coo_dec",
    "plx_value",
    "dist_st_value",
    "teff_st_value",
    "radius_st_value",
    "mass_st_value",
]


def analysis(StarCat4, distance_cut_in_pc):
    """
    This Function analyses the StarCat4.

    :param StarCat4:
    :type StarCat4:
    :param float distance_cut_in_pc:
    """

    la.spechistplot(
        [
            StarCat4["class_temp"][
                np.where(StarCat4["binary_flag"] == "False")
            ],
            StarCat4["class_temp"][np.where(StarCat4["binary_flag"] == "True")],
            StarCat4["class_temp"],
        ],
        ["Singles", "Binaries", "Total"],
        "spechist",
    )

    la.final_plot(
        [StarCat4["class_temp", "dist_st_value"]],
        ["LIFE-StarCat4"],
        distance_cut_in_pc,
    )

    within45deg = StarCat4[np.where(StarCat4["ecliptic_pm45deg"] == "True")]

    print(
        "if we cut for architecture reason everything above and below +45 deg, -45 deg respectively"
    )
    print(
        f"we are left with {np.round(100 * len(within45deg) / len(StarCat4), 1)}% of our stars"
    )
    print(
        f"this corresponds to a factor of {np.round(len(StarCat4) / len(within45deg), 1)} fewer"
    )
    remaining_singles = len(
        within45deg[np.where(within45deg["binary_flag"] == "False")]
    ) / len(within45deg)
    print(
        "from the remaining stars almost all(",
        np.round(100 * remaining_singles, 1),
        "%) are singles",
    )

    return


def prepare_dataframe_for_sns(table, paras, labels):
    """
    This function prepares the dataframe for seaborn boxplot.

    :param table: Table containing two data sets and the column
        names given in paras as well as class_temp column.
    :type table: astropy.table.table.Table
    :param paras: Column names of table with first list corresponding to first lables entry.
    :type paras: list(list(str))
    :param labels: Touple with lable for the two data sets.
    :type labels: list(str)
    :returns: Dataframe ready to be ingested by seaborn boxplot.
    :rtype: pandas.core.frame.DataFrame
    """
    seaborn_join = table[["class_temp"] + paras[0]].copy()
    seaborn_join["catalog"] = [labels[0] for j in range(len(seaborn_join))]
    temp = table[["class_temp"] + paras[1]].copy()
    temp["catalog"] = [labels[1] for j in range(len(temp))]
    seaborn_join.rename_columns(paras[0], paras[1])
    seaborn_join = ap.table.vstack([seaborn_join, temp])
    df = seaborn_join.to_pandas()
    return df


def snsplot(df, y, path):
    """
    This function creates a boxplot.

    :param df:
    :type df:
    :param y:
    :type y:
    :param str path:
    """

    plt.figure()
    sns.boxplot(
        data=df, x="class_temp", y=y, hue="catalog", order=["F", "G", "K", "M"]
    )
    #plt.savefig(path + ".png")
    plt.show()
    plt.figure(figsize=(2, 3))
    sns.boxplot(data=df, y=y, hue="catalog", legend=False)
    #plt.savefig(path + "total.png", bbox_inches="tight")
    plt.show()
    return


def threecatboxplot(data, para, labels):
    """
    This function creates a boxplot.

    :param data:
    :type data:
    :param para:
    :type para:
    """

    fig = plt.figure()

    # Creating axes instance
    ax = fig.add_axes([0, 0, 1, 1])
    plt.ylabel(para)
    # Creating plot
    bp = ax.boxplot(data, labels=labels)
    #plt.savefig(Path().plot + para + "_allcat_total.png", bbox_inches="tight")
    # show plot
    plt.show()
    return


def distance_boxplot_catalog_comparison(ax2):
    # don't have a join thingy so might need to adapt stuff
    # df=prepare_dataframe_for_sns(join_34,[para3,para4],labels)
    # I need class_temp for spectral type sorting. can do that too with the specsomething function

    seaborn_join = table[["class_temp"] + paras[0]].copy()
    seaborn_join["catalog"] = [labels[0] for j in range(len(seaborn_join))]
    temp = table[["class_temp"] + paras[1]].copy()
    temp["catalog"] = [labels[1] for j in range(len(temp))]
    seaborn_join.rename_columns(paras[0], paras[1])
    seaborn_join = ap.table.vstack([seaborn_join, temp])
    df = seaborn_join.to_pandas()

    snsplot(df, "dist_st_value")
    return


def boxplot_all_objects(catalogs: list[Any], labels, paras):
    y_axis = [
        "Stellar_Distance",
        "Stellar_Mass",
        "Stellar_Radius",
        "Stellar_Effective_Temperature",
    ]
    for i in range(len(y_axis)):
        # comparing all objects in the catalogs -> wait, I want only two cats -> add very old
        data = [catalogs[0][paras[0][i]], catalogs[1][paras[1][i]]]
        threecatboxplot(data, y_axis[i], labels)


def boxplot_common_objects(common_0_1, common_1_0, labels, paras):
    # comparing only objects that are in both catalogs

    # preparing the data frame
    if paras[0] != paras[1]:
        join_0_1 = ap.table.join(common_0_1, common_1_0, keys="main_id")
        save(
            [join_0_1],
            [f"join_{labels[0]}_{labels[1]}"],
            location="../" + Path().additional_data,
        )
        df_0_1 = prepare_dataframe_for_sns(
            join_0_1, [paras[0], paras[1]], [labels[0], labels[1]]
        )
    else:
        common_0_1["catalog"] = [labels[0] for j in range(len(common_0_1))]
        common_1_0["catalog"] = [labels[1] for j in range(len(common_1_0))]
        df_0_1 = ap.table.vstack([common_0_1, common_1_0])
        df_0_1 = df_0_1.to_pandas()

    for y in paras[0]:
        snsplot(df_0_1, y, Path().plot + f"sns_{labels[0]}_{labels[1]}{y}")


def ltc_compare(labels, paras, catalogs=[], paths=[]):
    """
    This function compares StarCat4 to StarCat3 and HWO catalog SPORES.

    :param bool load: Defaults to True, determines if tables are created or loaded.
    :param float distance_cut_in_pc: Defaults to 30., determines distance cut of StarCat4 if load=False.
    """

    if catalogs == []:
        for path in paths:
            catalogs.append(load([path], location=Path().additional_data)[0])

    common_0_1, common_1_0 = create_common(catalogs[0], catalogs[1])

    boxplot_common_objects(common_0_1, common_1_0, labels, paras)

    boxplot_all_objects(catalogs, labels, paras)

    la.final_plot(
        [
            catalogs[0]["class_temp", "dist_st_value"],
            catalogs[1]["class_temp", "dist_st_value"],
        ],
        labels,
        path=Path().plot + f"final_plot_{labels[0]}_{labels[1]}.png",
    )

    la.final_plot(
        [
            common_0_1["class_temp", "dist_st_value"],
            common_1_0["class_temp", "dist_st_value"],
        ],
        ["common" + labels[0], "common" + labels[1]],
        path=Path().plot + f"final_plot_common_{labels[0]}_{labels[1]}.png",
    )

    return


def compare():
    para4 = [
        "coo_ra",
        "coo_dec",
        "plx_value",
        "dist_st_value",
        "coo_gal_l",
        "coo_gal_b",
        "teff_st_value",
        "radius_st_value",
        "mass_st_value",
        "sep_phys_value",
        "mag_i_value",
        "mag_j_value",
    ]
    paras = [para4, para4]
    labels = ["StarCat5", "StarCat4"]
    paths = ["catalogs/StarCat5", "StarCat4"]
    ltc_compare(labels, paras, paths=paths)

def plot_cat_paras(paras,catalog_list,
                   label_list=["cat1","cat2"]):
    fig, ax = plt.subplots(
        figsize=(9, 6)
    )  # subplots so that I can overplot old version?
    for catalog,lab in zip(catalog_list,label_list):
        arr = catalog[paras[0], paras[1]]
        arr2 = arr[np.where(arr[paras[0]] != 1e20)]
        data = arr2[np.where(arr2[paras[1]] != 1e20)]

        ax.scatter(data[paras[0]], data[paras[1]], s=2, label = lab, alpha=0.5)
    # ax.set_yscale("log")

    ax.set_xlabel(paras[0])
    ax.set_ylabel(paras[1])
    plt.legend()
    plt.show()
