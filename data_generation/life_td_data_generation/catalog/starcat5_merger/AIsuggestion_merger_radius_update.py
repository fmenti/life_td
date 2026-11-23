from collections.abc import Sequence
from typing import Any, Literal

import importlib

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io.ascii import read
from astropy.table import MaskedColumn, Table, vstack
from astropy import units as u
from catalog.starcat5_merger.fcts_cat_merge import (
    get_cat2_in_cat1_match_info,
    add_cat2_to_cat1_by_name_or_coords,
)

# ... existing code ...


def get_radius(hpic: Table, starcat5: Table) -> float:
    """
    Estimate the default positional matching radius from ID-matched stars.

    The radius is derived by:
    1. selecting stars whose IDs occur in both HPIC and StarCat5,
    2. finding the nearest HPIC coordinate neighbor for each selected StarCat5
       star,
    3. checking that the nearest HPIC neighbor is the same row as the ID match,
    4. returning the largest validated nearest-neighbor distance.

    Returning the largest validated distance gives the smallest radius that still
    includes all verified ID matches.

    :param hpic: HPIC catalog.
    :type hpic: astropy.table.Table
    :param starcat5: StarCat5 catalog.
    :type starcat5: astropy.table.Table
    :returns: Default matching radius in arcseconds.
    :rtype: float
    :raises ValueError: If no ID matches are found or if any ID match is not the
        nearest coordinate match.
    """
    hpic_name_to_index = {}
    for hpic_index, name in enumerate(hpic["simbad_name"]):
        if name not in hpic_name_to_index:
            hpic_name_to_index[name] = hpic_index

    matched_starcat_indices = []
    matched_hpic_indices = []

    for starcat_index, name in enumerate(starcat5["main_id"]):
        if name in hpic_name_to_index:
            matched_starcat_indices.append(starcat_index)
            matched_hpic_indices.append(hpic_name_to_index[name])

    if len(matched_starcat_indices) == 0:
        raise ValueError("Cannot estimate radius: no ID-matched stars found.")

    matched_starcat_indices = np.asarray(matched_starcat_indices, dtype=int)
    matched_hpic_indices = np.asarray(matched_hpic_indices, dtype=int)

    hpic_coords = SkyCoord(
        ra=np.asarray(hpic["ra"]) * u.deg,
        dec=np.asarray(hpic["dec"]) * u.deg,
    )
    matched_starcat_coords = SkyCoord(
        ra=np.asarray(starcat5["coo_ra"][matched_starcat_indices]) * u.deg,
        dec=np.asarray(starcat5["coo_dec"][matched_starcat_indices]) * u.deg,
    )

    nearest_hpic_indices, sep2d, _ = matched_starcat_coords.match_to_catalog_sky(
        hpic_coords,
    )

    nearest_is_id_match = nearest_hpic_indices == matched_hpic_indices

    if not np.all(nearest_is_id_match):
        failed_indices = np.where(~nearest_is_id_match)[0]
        failed_starcat_indices = matched_starcat_indices[failed_indices]

        raise ValueError(
            "Cannot estimate radius: not all ID-matched StarCat5 stars have "
            "their ID-matched HPIC star as nearest coordinate neighbor. "
            f"Number of failed rows: {len(failed_starcat_indices)}."
        )

    radius = float(np.max(sep2d.arcsec))

    print("radius: ", radius)

    return radius


# ... existing code ...


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

    radius = get_radius(hpic, starcat5)

    match_info = get_cat2_in_cat1_match_info(
        name_cat1=hpic["simbad_name"],
        ra_cat1=hpic["ra"],
        dec_cat1=hpic["dec"],
        name_cat2=starcat5["main_id"],
        ra_cat2=starcat5["coo_ra"],
        dec_cat2=starcat5["coo_dec"],
        r_arcsec=radius,
    )

# ... existing code ...
