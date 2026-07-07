import numpy as np
from astropy.table import MaskedColumn, Table
from catalog.starcat5_merger.fcts_cat_merge import (
    get_mask_cat2_in_cat1,
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
    print(mask_cat2_in_cat1) #false true true
    print(starcat5[mask_cat2_in_cat1])
    print(starcat5[np.invert(mask_cat2_in_cat1)])

    # I want to get starcat5[mask_cat2_in_cat1] is only test3 object -> inverse?


