import matplotlib.pyplot as plt
import numpy as np  # Used for arrays
from astropy import io
from astropy.table import Table
from utils.io import Path, load
from catalog.starcat4 import starcat_creation
from catalog.starcat5 import (
    query_stars,
    query_children,
    query_systems,
)

path = Path().additional_data + "catalogs/"




def test_star_cat_looks_fine():
    starcat4 = starcat_creation(30,path=path)

    K_in_30 = starcat4[np.where(starcat4["class_temp"] == "K")]
    K_in_20 = K_in_30[np.where(K_in_30["dist_st_value"] < 20.0)]
    print(
        K_in_20
    )  # 126, felix reports 154 which is concerning because in ltc3 we had 245
    # find out which ones are missing and why, maybe start with G stars as there are fewer
    G_in_30 = starcat4[np.where(starcat4["class_temp"] == "G")]
    G_in_20 = G_in_30[np.where(G_in_30["dist_st_value"] < 20.0)]
    print(len(G_in_20))  # 34 felix reports 35 and in ltcv 65

    # formatting (csv file is correctly readable)
    # number of stars
    assert type(starcat4) == type(Table())


def test_star_cat_contains_specific_objects():
    starcat4 = starcat_creation(30,path=path)

    # list of important stars
    l_golden_targets = ['* alf Cen A', '* alf Cen B', '* alf Cen C',
                        'NAME Teegarden\'s Star', '* tau Cet', '* eps Eri',
                        '* eps Ind', 'HD  88230', '* omi02 Eri', '* sig Dra',
                        'HD 131977', '* eta Cas A', '* eta Cas B', 'TRAPPIST-1']
    # list of reasons why excluded
    # tbd read it out directly from code detail_criteria instead of copy pasting from hand
    dropout_reasons = ['has more than one sibling',
                       'has more than one sibling',
                       'has more than one sibling',
                       'system object but with found child object',
                       'companion does not fit spectral type criteria',
                       'multiple parents',
                       'companion does not fit spectral type criteria',
                       'system without child',
                       'multiple parents',
                       'multiple parents',
                       'multiple parents',
                       'system without child',
                       'system without child',
                       '']
    # check remaining stars present

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


#starcat5
def test_query_stars():
    #data
    colnames = ['main_id', 'coo_ra', 'coo_dec', 'sptype_string', 'plx_value',
                'dist_st_value', 'coo_gal_l', 'coo_gal_b', 'teff_st_value',
                'teff_ref', 'mass_st_value', 'mass_ref', 'radius_st_value',
                'radius_ref', 'binary_flag', 'binary_ref', 'mag_i_value',
                'mag_j_value', 'class_lum', 'class_temp', 'parent_main_id',
                'sep_ang_value']
    #execute
    query = query_stars("http://localhost:8080/tap",30)
    #assert
    assert len(query)>9000
    assert query.colnames == colnames

def test_query_children():
    #data
    colnames = ['child_main_id', 'child_type', 'parent_object_idref']
    #execute
    query = query_children("http://localhost:8080/tap")
    #assert
    assert len(query)>6000
    assert query.colnames == colnames

def test_query_systems():
    #data
    colnames = ['object_id', 'main_id', 'coo_ra', 'coo_dec', 'sptype_string',
                'plx_value', 'dist_st_value', 'coo_gal_l', 'coo_gal_b',
                'teff_st_value', 'teff_ref', 'mass_st_value', 'mass_ref',
                'radius_st_value', 'radius_ref', 'binary_flag', 'binary_ref',
                'mag_i_value', 'mag_j_value', 'class_lum', 'class_temp',
                'parent_main_id', 'sep_ang_value']
    #execute
    query = query_systems("http://localhost:8080/tap",30)
    #assert
    assert len(query)>2000
    assert query.colnames == colnames


import subprocess
import sys
from pathlib import Path as pp


def test_runs_as_script():
    here = pp(__file__).resolve().parent  # .../integration_tests/post_db_creation
    script = here.parents[
                 1] / "catalog" / "starcat5.py"  # up 2 levels, then catalog/starcat5.py

    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=here,  # run "as if" you launched it from post_db_creation
        check=False,
    )
    assert result.returncode == 0

def test_processing_doesnt_loose_objects():
    # rewrite this into plots instead of doubling code from main
    # run starcat5 then plot tables with different flags and assert lengths

    #StarCat5 = how to run main?


    #assert len(flag_non_ms)==len(stars_with_ub)
    #assert len(singles) > 6000
    #assert len(multiples) > 4000
    assert False
