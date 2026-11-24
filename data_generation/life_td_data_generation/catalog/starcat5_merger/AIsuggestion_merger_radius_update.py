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

def get_radius_mismatch_table(
    hpic: Table,
    starcat5: Table,
    matched_starcat_indices: np.ndarray,
    matched_hpic_indices: np.ndarray,
    nearest_hpic_indices: np.ndarray,
    nearest_distances: u.Quantity,
) -> Table:
    """
    Build a diagnostics table for ID-matched stars whose nearest HPIC coordinate
    match is not their ID-matched HPIC row.

    :param hpic: HPIC catalog.
    :type hpic: astropy.table.Table
    :param starcat5: StarCat5 catalog.
    :type starcat5: astropy.table.Table
    :param matched_starcat_indices: StarCat5 indices for rows with ID matches.
    :type matched_starcat_indices: numpy.ndarray
    :param matched_hpic_indices: HPIC indices of the corresponding ID matches.
    :type matched_hpic_indices: numpy.ndarray
    :param nearest_hpic_indices: HPIC indices of nearest coordinate matches.
    :type nearest_hpic_indices: numpy.ndarray
    :param nearest_distances: Angular distances to nearest HPIC matches.
    :type nearest_distances: astropy.units.Quantity
    :returns: Diagnostics table for mismatched nearest-neighbor cases.
    :rtype: astropy.table.Table
    """
    nearest_is_id_match = nearest_hpic_indices == matched_hpic_indices
    failed_indices = np.where(~nearest_is_id_match)[0]

    failed_starcat_indices = matched_starcat_indices[failed_indices]
    failed_id_matched_hpic_indices = matched_hpic_indices[failed_indices]
    failed_nearest_hpic_indices = nearest_hpic_indices[failed_indices]

    failed_starcat_coords = SkyCoord(
        ra=np.asarray(starcat5["coo_ra"][failed_starcat_indices]) * u.deg,
        dec=np.asarray(starcat5["coo_dec"][failed_starcat_indices]) * u.deg,
    )
    failed_id_matched_hpic_coords = SkyCoord(
        ra=np.asarray(hpic["ra"][failed_id_matched_hpic_indices]) * u.deg,
        dec=np.asarray(hpic["dec"][failed_id_matched_hpic_indices]) * u.deg,
    )

    id_matched_distances = failed_starcat_coords.separation(
        failed_id_matched_hpic_coords,
    )

    return Table(
        {
            "starcat_name": starcat5["main_id"][failed_starcat_indices],
            "id_matched_hpic_name": hpic["simbad_name"][
                failed_id_matched_hpic_indices
            ],
            "nearest_hpic_name": hpic["simbad_name"][failed_nearest_hpic_indices],
            "id_matched_dist_arcsec": id_matched_distances.arcsec,
            "nearest_hpic_dist_arcsec": nearest_distances[
                failed_indices
            ].arcsec,
            "starcat_index": failed_starcat_indices,
            "id_matched_hpic_index": failed_id_matched_hpic_indices,
            "nearest_hpic_index": failed_nearest_hpic_indices,
        }
    )


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
    hpic_name_to_index = {} # dictionary [name: index]
    for hpic_index, name in enumerate(hpic["simbad_name"]):
        if name not in hpic_name_to_index:
            hpic_name_to_index[name] = hpic_index

    matched_starcat_indices = []
    matched_hpic_indices = []
    # if name in hpic and starcat then add index to matched lists each

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
        radius_mismatch_table = get_radius_mismatch_table(
            hpic=hpic,
            starcat5=starcat5,
            matched_starcat_indices=matched_starcat_indices,
            matched_hpic_indices=matched_hpic_indices,
            nearest_hpic_indices=nearest_hpic_indices,
            nearest_distances=sep2d,
        )

        print(radius_mismatch_table)
        print(radius_mismatch_table["starcat_name","nearest_hpic_name",
              "nearest_hpic_dist_arcsec","id_matched_dist_arcsec"])

        raise ValueError(
            "Cannot estimate radius: not all ID-matched StarCat5 stars have "
            "their ID-matched HPIC star as nearest coordinate neighbor. "
            f"Number of failed rows: {len(radius_mismatch_table)}."
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
