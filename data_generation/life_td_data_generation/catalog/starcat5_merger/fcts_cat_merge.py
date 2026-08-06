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

def get_cat2_in_cat1_match_info(
    name_cat1,
    ra_cat1,
    dec_cat1,
    name_cat2,
    ra_cat2,
    dec_cat2,
    r_arcsec,
):
    """
    Match cat2 rows to cat1 rows, first by name and then by sky coordinates.

    Parameters
    ----------
    name_cat1, name_cat2 : array_like
        Object names.
    ra_cat1, ra_cat2 : array_like
        Right ascension values in degrees, or Astropy quantities/columns
        with angular units.
    dec_cat1, dec_cat2 : array_like
        Declination values in degrees, or Astropy quantities/columns
        with angular units.
    r_arcsec : float
        Matching radius in arcsec.

    Returns
    -------
    match_info : dict
        Dictionary with one entry per cat2 row:
        - ``mask_cat2_in_cat1``: True if cat2 row matched cat1.
        - ``cat1_index_for_cat2``: matched cat1 row index, or -1.
        - ``match_by_name``: True if matched by name.
        - ``match_by_coords``: True if matched by coordinates.
        - ``sep_arcsec``: coordinate separation in arcsec for coordinate
          nearest-neighbor match. Name-only matches can have non-zero values
          because this is computed independently from the name matching.
    """

    cat1_index_for_cat2 = np.full(len(name_cat2), -1, dtype=int)
    match_by_name = np.zeros(len(name_cat2), dtype=bool)
    match_by_coords = np.zeros(len(name_cat2), dtype=bool)

    # 1) Name matching: cat2 -> cat1
    cat1_name_to_index = {}
    for ii, name in enumerate(name_cat1):
        if name not in cat1_name_to_index:
            cat1_name_to_index[name] = ii

    for ii, name in enumerate(name_cat2):
        if name in cat1_name_to_index:
            cat1_index_for_cat2[ii] = cat1_name_to_index[name]
            match_by_name[ii] = True

    # 2) Coordinate matching for cat2 rows without name match
    coords_cat1 = SkyCoord(
        ra=_as_degree_quantity(ra_cat1),
        dec=_as_degree_quantity(dec_cat1),
    )
    coords_cat2 = SkyCoord(
        ra=_as_degree_quantity(ra_cat2),
        dec=_as_degree_quantity(dec_cat2),
    )

    # astropy SkyCoord matching
    idx_cat1, sep2d, _ = coords_cat2.match_to_catalog_sky(coords_cat1)

    unmatched_by_name = ~match_by_name
    close_enough = sep2d.arcsec <= r_arcsec
    coord_matches = unmatched_by_name & close_enough

    cat1_index_for_cat2[coord_matches] = idx_cat1[coord_matches]
    match_by_coords[coord_matches] = True

    mask_cat2_in_cat1 = match_by_name | match_by_coords

    return {
        "mask_cat2_in_cat1": mask_cat2_in_cat1,
        "cat1_index_for_cat2": cat1_index_for_cat2,
        "match_by_name": match_by_name,
        "match_by_coords": match_by_coords,
        "sep_arcsec": sep2d.arcsec,
    }


def get_mask_cat2_in_cat1(
    name_cat1,
    ra_cat1,
    dec_cat1,
    name_cat2,
    ra_cat2,
    dec_cat2,
    r_arcsec,
):
    """
    For each target in cat2, check if it is in cat1.

    Matching is done first by name, then by coordinates.

    Returns
    -------
    mask_cat2_in_cat1 : numpy.ndarray
        Boolean mask with same length as cat2. True if a target in cat2 is
        also in cat1.
    """

    match_info = get_cat2_in_cat1_match_info(
        name_cat1=name_cat1,
        ra_cat1=ra_cat1,
        dec_cat1=dec_cat1,
        name_cat2=name_cat2,
        ra_cat2=ra_cat2,
        dec_cat2=dec_cat2,
        r_arcsec=r_arcsec,
    )

    return match_info["mask_cat2_in_cat1"]

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
    match_info=None,
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
    match_info : dict or None
        Optional precomputed output from get_cat2_in_cat1_match_info. If given,
        matching is not recomputed.

    Returns
    -------
    joined : astropy.table.Table
        cat1 with selected cat2 columns appended.
    """

    if cat2_cols is None:
        cat2_cols = list(cat2.colnames)

    if match_info is None:
        match_info = get_cat2_in_cat1_match_info(
            name_cat1=cat1[name_cat1_col],
            ra_cat1=cat1[ra_cat1_col],
            dec_cat1=cat1[dec_cat1_col],
            name_cat2=cat2[name_cat2_col],
            ra_cat2=cat2[ra_cat2_col],
            dec_cat2=cat2[dec_cat2_col],
            r_arcsec=r_arcsec,
        )

    mask_cat2_in_cat1 = match_info["mask_cat2_in_cat1"]
    cat1_index_for_cat2 = match_info["cat1_index_for_cat2"]

    cat2_index_for_cat1 = np.full(len(cat1), -1, dtype=int)
    cat1_has_match = np.zeros(len(cat1), dtype=bool)
    cat1_match_by_name = np.zeros(len(cat1), dtype=bool)
    cat1_match_by_coords = np.zeros(len(cat1), dtype=bool)
    cat1_match_sep_arcsec = np.ma.masked_all(len(cat1), dtype=float)

    for cat2_index in np.where(mask_cat2_in_cat1)[0]:
        cat1_index = cat1_index_for_cat2[cat2_index]

        if cat1_index < 0:
            continue

        # Keep the first match if multiple cat2 rows point to the same cat1 row.
        if cat1_has_match[cat1_index]:
            # maybe change to keep the first non null match of the value
            continue

        cat2_index_for_cat1[cat1_index] = cat2_index
        cat1_has_match[cat1_index] = True
        cat1_match_by_name[cat1_index] = match_info["match_by_name"][cat2_index]
        cat1_match_by_coords[cat1_index] = match_info["match_by_coords"][cat2_index]
        cat1_match_sep_arcsec[cat1_index] = match_info["sep_arcsec"][cat2_index]

    cat2_matched = Table(masked=True)

    for col in cat2_cols:
        new_col_name = col
        if new_col_name in cat1.colnames:
            new_col_name = f"{col}{cat2_suffix}"

        cat2_matched[new_col_name] = np.ma.masked_all(
            len(cat1),
            dtype=cat2[col].dtype,
        )

        cat1_rows_with_match = cat2_index_for_cat1 >= 0
        cat2_matched[new_col_name][cat1_rows_with_match] = cat2[col][
            cat2_index_for_cat1[cat1_rows_with_match]
        ]

    cat2_matched["cat2_match_by_name"] = cat1_match_by_name
    cat2_matched["cat2_match_by_coords"] = cat1_match_by_coords
    cat2_matched["cat2_match_sep_arcsec"] = cat1_match_sep_arcsec
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

