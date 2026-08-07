import numpy as np
from astropy.table import Table, MaskedColumn
from catalog.starcat5_merger.fcts_cat_merge import (
    get_mask_cat2_in_cat1,
    get_cat2_in_cat1_match_info,
    add_cat2_to_cat1_by_name_or_coords,
    _invert_cat2_matches_to_cat1,
)




def test_get_mask_cat2_in_cat1():
    # data
    hpic = Table(
        (
            np.array(
                ["* eta UMa", "HD 127024", "test1","test2","LP  137-54"]
            ),
            np.array([206.88435898896, 217.95318702116, 10.0,1.,252.89855630184]),
            np.array([49.31320798716, -64.28385140399, 3.0,-1,51.63245588655]),
        ),
        names=(
            "simbad_name",
            "ra",
            "dec",
        ),
        dtype=[object, float, float],
    )

    starcat5 = Table(
        (
            np.array(
                ["test3","* eta UMa", "HD 127024","LP  137-54"]
            ),
            np.array([5.0, 206.88515734206297, 217.95318811753998,252.89855562033]),
            np.array([1.0, 49.31326672942533, -64.28385151882,51.6324558137]),
        ),
        names=(
            "main_id",
            "coo_ra",
            "coo_dec",
        ),
        dtype=[object, float, float],
    )

    radius = 50.

    # execute
    mask_cat2_in_cat1 = get_mask_cat2_in_cat1(name_cat1=hpic["simbad_name"],
                          ra_cat1=hpic["ra"],
                          dec_cat1=hpic["dec"],
                          name_cat2=starcat5["main_id"],
                          ra_cat2=starcat5["coo_ra"],
                          dec_cat2=starcat5["coo_dec"],
                          r_arcsec=radius)

    # assert
    assert starcat5[np.invert(mask_cat2_in_cat1)]["main_id"][0] == "test3"
    assert len(starcat5[np.invert(mask_cat2_in_cat1)]) == 1

def test_get_mask_cat2_in_cat1_match_info():
    # data
    hpic = Table(
        (
            np.array(
                ["* eta UMa", "HD 127024", "test1", "test2", "LP  137-54"]
            ),
            np.array(
                [206.88435898896, 217.95318702116, 10.0, 1., 252.89855630184]),
            np.array(
                [49.31320798716, -64.28385140399, 3.0, -1, 51.63245588655]),
        ),
        names=(
            "simbad_name",
            "ra",
            "dec",
        ),
        dtype=[object, float, float],
    )

    starcat5 = Table(
        (
            np.array(
                ["test3", "* eta UMa", "HD 127024", "LP wrong name"]
            ),
            np.array(
                [5.0, 206.88515734206297, 217.95318811753998, 252.89855562033]),
            np.array([1.0, 49.31326672942533, -64.28385151882, 51.6324558137]),
        ),
        names=(
            "main_id",
            "coo_ra",
            "coo_dec",
        ),
        dtype=[object, float, float],
    )

    radius = 50.

    # execute
    match_info = get_cat2_in_cat1_match_info(name_cat1=hpic["simbad_name"],
                                              ra_cat1=hpic["ra"],
                                              dec_cat1=hpic["dec"],
                                              name_cat2=starcat5["main_id"],
                                              ra_cat2=starcat5["coo_ra"],
                                              dec_cat2=starcat5["coo_dec"],
                                              r_arcsec=radius)

    # assert
    assert list(match_info["match_by_name"]) == [False,  True,  True,  False]
    assert list(match_info["match_by_coords"]) == [False, False, False, True]

def test_add_cat2_to_cat1_by_name_or_coords():
    # data
    hpic = Table(
        (
            np.array(
                ["* eta UMa", "HD 127024", "test1", "test2", "LP  137-54"]
            ),
            np.array(
                [206.88435898896, 217.95318702116, 10.0, 1., 252.89855630184]),
            np.array(
                [49.31320798716, -64.28385140399, 3.0, -1, 51.63245588655]),
        ),
        names=(
            "simbad_name",
            "ra",
            "dec",
        ),
        dtype=[object, float, float],
    )

    starcat5 = Table(
        (
            np.array(
                ["test3", "* eta UMa", "HD 127024", "LP wrong name"]
            ),
            np.array(
                [5.0, 206.88515734206297, 217.95318811753998, 252.89855562033]),
            np.array([1.0, 49.31326672942533, -64.28385151882, 51.6324558137]),
        ),
        names=(
            "main_id",
            "coo_ra",
            "coo_dec",
        ),
        dtype=[object, float, float],
    )

    starcat5["mag_u_sdss_value"] = MaskedColumn(dtype=float, length=len(starcat5), mask = True)
    starcat5["mag_u_sdss_value"][2] = 2.2
    starcat5["mag_u_sdss_value"][3] = 3.3

    radius = 50.

    hpic_mag_u_joined = add_cat2_to_cat1_by_name_or_coords(
        cat1 = hpic,
        cat2 = starcat5,
        name_cat1_col = "simbad_name",
        ra_cat1_col = "ra",
        dec_cat1_col = "dec",
        name_cat2_col = "main_id",
        ra_cat2_col = "coo_ra",
        dec_cat2_col = "coo_dec",
        r_arcsec = radius,
        cat2_cols = [
            "main_id",
            "coo_ra",
            "coo_dec",
            "mag_u_sdss_value",
        ],

    )

    new_columns = [
        "main_id",
        "coo_ra",
        "coo_dec",
        "mag_u_sdss_value",
        "cat2_match_by_name",
        "cat2_match_by_coords",
        "cat2_match_sep_arcsec",
        "cat2_match_index",
        ]


    assert len(hpic_mag_u_joined) == len(hpic)
    assert list(hpic_mag_u_joined.columns) == list(hpic.columns) + new_columns

def test_get_cat2_in_cat1_match_info():
    # data
    hpic = Table(
        (
            np.array(
                ["* eta UMa", "HD 127024", "test1", "test2", "LP  137-54"]
            ),
            np.array(
                [206.88435898896, 217.95318702116, 10.0, 1., 252.89855630184]),
            np.array(
                [49.31320798716, -64.28385140399, 3.0, -1, 51.63245588655]),
        ),
        names=(
            "simbad_name",
            "ra",
            "dec",
        ),
        dtype=[object, float, float],
    )

    starcat5 = Table(
        (
            np.array(
                ["test3", "* eta UMa", "HD 127024", "LP wrong name"]
            ),
            np.array(
                [5.0, 206.88515734206297, 217.95318811753998, 252.89855562033]),
            np.array([1.0, 49.31326672942533, -64.28385151882, 51.6324558137]),
        ),
        names=(
            "main_id",
            "coo_ra",
            "coo_dec",
        ),
        dtype=[object, float, float],
    )

    radius = 50.

    cat2_match_info = get_cat2_in_cat1_match_info(
        name_cat1=hpic["simbad_name"],
        ra_cat1=hpic["ra"],
        dec_cat1=hpic["dec"],
        name_cat2=starcat5["main_id"],
        ra_cat2=starcat5["coo_ra"],
        dec_cat2=starcat5["coo_dec"],
        r_arcsec=radius,
    )

    print(cat2_match_info)

    hpic_index_for_starcat = [-1,  0,  1,  4]
    name_match = [False,  True,  True, False]
    coord_match = [False, False, False,  True]
    mask = [False,  True,  True,  True]

    cat2_match_info["mask_cat2_in_cat1"] == mask
    cat2_match_info["cat1_index_for_cat2"] == hpic_index_for_starcat
    cat2_match_info["match_by_name"] == name_match
    cat2_match_info["match_by_coords"] == coord_match


def test__invert_cat2_matches_to_cat1():
    # data
    hpic = Table(
        (
            np.array(
                ["* eta UMa", "HD 127024", "test1", "test2", "LP  137-54"]
            ),
            np.array(
                [206.88435898896, 217.95318702116, 10.0, 1., 252.89855630184]),
            np.array(
                [49.31320798716, -64.28385140399, 3.0, -1, 51.63245588655]),
        ),
        names=(
            "simbad_name",
            "ra",
            "dec",
        ),
        dtype=[object, float, float],
    )

    starcat5 = Table(
        (
            np.array(
                ["test3", "* eta UMa", "HD 127024", "LP wrong name"]
            ),
            np.array(
                [5.0, 206.88515734206297, 217.95318811753998, 252.89855562033]),
            np.array([1.0, 49.31326672942533, -64.28385151882, 51.6324558137]),
        ),
        names=(
            "main_id",
            "coo_ra",
            "coo_dec",
        ),
        dtype=[object, float, float],
    )


    radius = 50.

    match_info = get_cat2_in_cat1_match_info(
        name_cat1=hpic["simbad_name"],
        ra_cat1=hpic["ra"],
        dec_cat1=hpic["dec"],
        name_cat2=starcat5["main_id"],
        ra_cat2=starcat5["coo_ra"],
        dec_cat2=starcat5["coo_dec"],
        r_arcsec=radius,
    )

    cat1_match_info = _invert_cat2_matches_to_cat1(
        match_info,
        len(hpic),
    )

    starcat5_index_for_hpic = [ 1,  2, -1, -1,  3]
    name_match = [ True,  True, False, False, False]
    coord_match = [False, False, False, False,  True]
    mask = [False, False,  True,  True, False]

    cat1_match_info["cat2_match_index"] == starcat5_index_for_hpic
    cat1_match_info["cat2_match_by_name"] == name_match
    cat1_match_info["cat2_match_by_coords"] == coord_match
    cat1_match_info["cat2_match_sep_arcsec"].mask == mask

