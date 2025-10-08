import matplotlib.pyplot as plt
import numpy as np  # Used for arrays
from astropy import io
from astropy.table import Table
from utils.io import Path, load
from utils.starcat4 import starcat_creation


def test_star_cat_looks_fine():
    starcat4 = starcat_creation(30)
    # why do I have difficulties running
    # this function and why do I not get the print I want? because was called integration test instead of test...
    # print(starcat4)
    K_in_30 = starcat4[np.where(starcat4["class_temp"] == "K")]
    K_in_20 = K_in_30[np.where(K_in_30["dist_st_value"] < 20.0)]
    print(
        K_in_20
    )  # 126, felix reports 154 which is concerning because in ltc3 we had 245
    # find out which ones are missing and why, maybe start with G stars as there are fewer
    G_in_30 = starcat4[np.where(starcat4["class_temp"] == "G")]
    G_in_20 = G_in_30[np.where(G_in_30["dist_st_value"] < 20.0)]
    print(len(G_in_20))  # 34 felix reports 35 and in ltcv 65

    assert type(starcat4) == type(Table())


def test_star_cat_contains_specific_objects():
    starcat4 = starcat_creation(30)

    assert len(starcat4[np.where(starcat4["main_id"] == "TRAPPIST-1")]) == 1


def teff_vs_dist_plot(dist_data, teff_data):
    # plt.figure()
    fig, ax = plt.subplots(figsize=(9, 6))
    # TBD take two data sets and use different colors -> tell AI to do so
    ax.scatter(dist_data, teff_data, s=2)
    ax.set_yscale("log")

    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel("Temperature [K]")
    # plt.savefig(Path().plot+'/teff_vs_dist', dpi=300)
    plt.show()
    return


def teff_and_dist_db_data():
    # loading the correct table
    [table] = load(["star_basic"], location=Path().data)
    # exctracting the correct columns
    arr = table["dist_st_value", "teff_st_value"]
    arr2 = arr[np.where(arr["dist_st_value"] != 1e20)]
    data = arr2[np.where(arr2["teff_st_value"] != 1e20)]
    return data


def teff_and_dist_cat_data(number, dist_colname, temp_colname):
    # loading the correct table
    table = io.ascii.read(
        Path().additional_data + "catalogs/LIFE-StarCat" + number + ".csv"
    )
    # exctracting the correct columns
    arr = table[dist_colname, temp_colname]
    arr2 = arr[np.where(arr[dist_colname] != 1e20)]
    data = arr2[np.where(arr2[temp_colname] != 1e20)]
    return data


def test_teff_vs_dist():
    # what changes with the different versions
    # (old version one color of dots, new version other color of dots,
    # visual comparison of areas where only one color)

    # data
    db_data = teff_and_dist_db_data()
    # now do the same with cat3 and cat4
    ltc3_dist_colname = "distance"
    ltc3_temp_colname = "mod_Teff"
    ltc3_data = teff_and_dist_cat_data(
        "3", ltc3_dist_colname, ltc3_temp_colname
    )

    ltc4_dist_colname = "dist_st_value"
    ltc4_temp_colname = "teff_st_value"
    ltc4_data = teff_and_dist_cat_data(
        "4", ltc4_dist_colname, ltc4_temp_colname
    )

    # function
    teff_vs_dist_plot(db_data["dist_st_value"], db_data["teff_st_value"])

    teff_vs_dist_plot(
        ltc3_data[ltc3_dist_colname], ltc3_data[ltc3_temp_colname]
    )

    teff_vs_dist_plot(
        ltc4_data[ltc4_dist_colname], ltc4_data[ltc4_temp_colname]
    )
