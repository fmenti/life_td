from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from astropy.table import hstack, Table

def _as_degree_quantity(values):
    """
    Return values as an Astropy Quantity in degrees.

    If values already has angular units, convert to degrees.
    If values has no unit, assume it is given in degrees.
    """
    unit = getattr(values, "unit", None)

    if unit is None:
        return np.asarray(values) * u.deg

    return values.to(u.deg)


def get_mask_cat2_in_cat1(name_cat1, ra_cat1, dec_cat1, name_cat2, ra_cat2, dec_cat2, r_arcsec):
    """
    For each target in cat2, check if in cat1 (first by name, then by coords)


    Parameters
    ----------
    name_cat1, name_cat1: array_lie
    	Star name as string
    ra_cat1, ra_cat2: array_like
        Right ascension values in degrees.
    dec_cat1, dec_cat2: array_like
        Declination values in degrees.
    r_arcsec: float
    	Matching radius in arcsec

    Returns
    -------
    mask_cat2_in_cat1 : numpy.ndarray, same lenght as cat2
        Boolean mask. True if a target in cat2 is also in cat1
    """

    # 1) Names
    mask_cat2_in_cat1_name = np.isin(name_cat2, name_cat1)

    # 2) Coordinates
    coords_cat1 = SkyCoord(ra=_as_degree_quantity(ra_cat1),
                           dec=_as_degree_quantity(dec_cat1))
    coords_cat2 = SkyCoord(ra=_as_degree_quantity(ra_cat2),
                           dec=_as_degree_quantity(dec_cat2))

    mask_cat2_in_cat1_coords = np.zeros(len(name_cat2), dtype=bool)

    for ii, c in enumerate(coords_cat2):
        if not mask_cat2_in_cat1_name[ii]:
            sep = c.separation(coords_cat1)
            if np.min(sep.arcsec) <= r_arcsec :
                mask_cat2_in_cat1_coords[ii] = True

    mask_cat2_in_cat1 = mask_cat2_in_cat1_name | mask_cat2_in_cat1_coords

    return mask_cat2_in_cat1

def get_cat2_match_indices_in_cat1(
    name_cat1,
    ra_cat1,
    dec_cat1,
    name_cat2,
    ra_cat2,
    dec_cat2,
    r_arcsec,
):
    """
    For each row in cat1, find the matching row in cat2.

    Matching is done first by name, then by coordinates for rows that
    did not get a name match.

    Returns
    -------
    cat2_index_for_cat1 : numpy.ndarray
        Integer array with same length as cat1.
        Value is the matching row index in cat2, or -1 if no match exists.

    match_by_name : numpy.ndarray
        Boolean array with same length as cat1. True where match was by name.

    match_by_coords : numpy.ndarray
        Boolean array with same length as cat1. True where match was by coordinates.
    """

    cat2_index_for_cat1 = np.full(len(name_cat1), -1, dtype=int)
    match_by_name = np.zeros(len(name_cat1), dtype=bool)
    match_by_coords = np.zeros(len(name_cat1), dtype=bool)

    # 1) Name matching
    cat2_name_to_index = {}
    for ii, name in enumerate(name_cat2):
        if name not in cat2_name_to_index:
            cat2_name_to_index[name] = ii

    for ii, name in enumerate(name_cat1):
        if name in cat2_name_to_index:
            cat2_index_for_cat1[ii] = cat2_name_to_index[name]
            match_by_name[ii] = True

    # 2) Coordinate matching for rows without name match
    coords_cat1 = SkyCoord(
        ra=_as_degree_quantity(ra_cat1),
        dec=_as_degree_quantity(dec_cat1),
    )
    coords_cat2 = SkyCoord(
        ra=_as_degree_quantity(ra_cat2),
        dec=_as_degree_quantity(dec_cat2),
    )

    idx_cat2, sep2d, _ = coords_cat1.match_to_catalog_sky(coords_cat2)

    unmatched_by_name = ~match_by_name
    close_enough = sep2d.arcsec <= r_arcsec
    coord_matches = unmatched_by_name & close_enough

    cat2_index_for_cat1[coord_matches] = idx_cat2[coord_matches]
    match_by_coords[coord_matches] = True

    return cat2_index_for_cat1, match_by_name, match_by_coords


def add_cat2_to_cat1_by_name_or_coords(
    cat1,
    cat2,
    name_cat1_col,
    ra_cat1_col,
    dec_cat1_col,
    name_cat2_col,
    ra_cat2_col,
    dec_cat2_col,
    r_arcsec,
    cat2_cols=None,
    cat2_suffix="_cat2",
):
    """
    Add columns from cat2 to cat1 using a left-join-like match.

    Rows are matched first by name, then by nearest sky coordinate within
    r_arcsec. All rows from cat1 are preserved.

    Parameters
    ----------
    cat1, cat2 : astropy.table.Table
        Catalogs to combine.
    name_cat1_col, name_cat2_col : str
        Name/id columns used for exact name matching.
    ra_cat1_col, dec_cat1_col, ra_cat2_col, dec_cat2_col : str
        Coordinate columns in degrees or with Astropy angular units.
    r_arcsec : float
        Maximum coordinate matching radius in arcseconds.
    cat2_cols : list[str] or None
        Columns from cat2 to append. If None, all cat2 columns are appended.
    cat2_suffix : str
        Suffix added to appended cat2 columns to avoid name collisions.

    Returns
    -------
    joined : astropy.table.Table
        cat1 with selected cat2 columns appended.
    """

    if cat2_cols is None:
        cat2_cols = list(cat2.colnames)

    cat2_index_for_cat1, match_by_name, match_by_coords = get_cat2_match_indices_in_cat1(
        name_cat1=cat1[name_cat1_col],
        ra_cat1=cat1[ra_cat1_col],
        dec_cat1=cat1[dec_cat1_col],
        name_cat2=cat2[name_cat2_col],
        ra_cat2=cat2[ra_cat2_col],
        dec_cat2=cat2[dec_cat2_col],
        r_arcsec=r_arcsec,
    )

    cat2_matched = Table(masked=True)

    for col in cat2_cols:
        new_col_name = col
        if new_col_name in cat1.colnames:
            new_col_name = f"{col}{cat2_suffix}"

        cat2_matched[new_col_name] = np.ma.masked_all(
            len(cat1),
            dtype=cat2[col].dtype,
        )

        has_match = cat2_index_for_cat1 >= 0
        cat2_matched[new_col_name][has_match] = cat2[col][cat2_index_for_cat1[has_match]]

    cat2_matched["cat2_match_by_name"] = match_by_name
    cat2_matched["cat2_match_by_coords"] = match_by_coords
    cat2_matched["cat2_match_index"] = cat2_index_for_cat1

    joined = hstack([cat1, cat2_matched])

    return joined



def nearest_neighbor_distances_units(ra, dec):
    """
    Compute nearest-neighbor angular distances for a catalog of stars.

    Parameters
    ----------
    ra : array_like
        Right ascension values in degrees, or Astropy quantities/columns with angular units.
    dec : array_like
        Declination values in degrees, or Astropy quantities/columns with angular units.

    Returns
    -------
    distances_arcsec : numpy.ndarray
        Array of nearest-neighbor distances in arcseconds,
        same length as input arrays.
    """
    stars = SkyCoord(
        ra=_as_degree_quantity(ra),
        dec=_as_degree_quantity(dec),
    )

    if len(stars) < 2:
        return np.full(len(stars), np.nan)

    # nthneighbor=2 because the closest match in the same catalog is the object itself.
    _, sep2d, _ = stars.match_to_catalog_sky(stars, nthneighbor=2)

    return sep2d.arcsec

