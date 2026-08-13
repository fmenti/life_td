from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack


MatchInfo = dict[str, np.ndarray]


def _as_degree_quantity(values: Any) -> u.Quantity:
    """
    Return angular values as an Astropy quantity in degrees.

    If ``values`` already has angular units, it is converted to degrees. If it
    has no unit, it is assumed to be given in degrees.

    :param values: Angular values, either unitless or with Astropy units.
    :type values: object
    :returns: Angular values as degree quantity.
    :rtype: astropy.units.Quantity
    """
    unit = getattr(values, "unit", None)

    if unit is None:
        return np.asarray(values) * u.deg

    return values.to(u.deg)


def _make_skycoord(ra: Any, dec: Any) -> SkyCoord:
    """
    Build a sky coordinate object from right ascension and declination.

    This function does not seem necessary but also not too bad if kept.

    :param ra: Right ascension values in degrees or with angular units.
    :type ra: object
    :param dec: Declination values in degrees or with angular units.
    :type dec: object
    :returns: Sky coordinates.
    :rtype: astropy.coordinates.SkyCoord
    """
    return SkyCoord(
        ra=_as_degree_quantity(ra),
        dec=_as_degree_quantity(dec),
    )


def _match_cat2_to_cat1_by_name(
    name_cat1: Any,
    name_cat2: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Match rows from catalog 2 to catalog 1 by exact object name.

    If a name occurs multiple times in catalog 1, the first occurrence is used.
    This preserves the previous first-match behavior.
    A bit difficult to understand with the loops and dictionaries.

    :param name_cat1: Object names in catalog 1.
    :type name_cat1: object
    :param name_cat2: Object names in catalog 2.
    :type name_cat2: object
    :returns: Catalog-1 indices for catalog-2 rows and name-match mask.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    cat1_index_for_cat2 = np.full(len(name_cat2), -1, dtype=int)
    match_by_name = np.zeros(len(name_cat2), dtype=bool)

    cat1_name_to_index: dict[Any, int] = {}
    for cat1_index, name in enumerate(name_cat1):
        if name not in cat1_name_to_index:
            cat1_name_to_index[name] = cat1_index

    for cat2_index, name in enumerate(name_cat2):
        if name in cat1_name_to_index:
            cat1_index_for_cat2[cat2_index] = cat1_name_to_index[name]
            match_by_name[cat2_index] = True

    return cat1_index_for_cat2, match_by_name


def _match_cat2_to_cat1_by_coords(
    ra_cat1: Any,
    dec_cat1: Any,
    ra_cat2: Any,
    dec_cat2: Any,
    r_arcsec: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Match rows from catalog 2 to nearest catalog-1 coordinates.

    :param ra_cat1: Right ascension values of catalog 1.
    :type ra_cat1: object
    :param dec_cat1: Declination values of catalog 1.
    :type dec_cat1: object
    :param ra_cat2: Right ascension values of catalog 2.
    :type ra_cat2: object
    :param dec_cat2: Declination values of catalog 2.
    :type dec_cat2: object
    :param r_arcsec: Maximum accepted matching distance in arcseconds.
    :type r_arcsec: float
    :returns: Nearest catalog-1 indices, coordinate-match mask, separations.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    coords_cat1 = _make_skycoord(ra_cat1, dec_cat1)
    coords_cat2 = _make_skycoord(ra_cat2, dec_cat2)

    cat1_indices, sep2d, _ = coords_cat2.match_to_catalog_sky(coords_cat1)
    sep_arcsec = sep2d.arcsec
    match_by_coords = sep_arcsec <= r_arcsec

    return cat1_indices, match_by_coords, sep_arcsec


def get_cat2_in_cat1_match_info(
    name_cat1: Any,
    ra_cat1: Any,
    dec_cat1: Any,
    name_cat2: Any,
    ra_cat2: Any,
    dec_cat2: Any,
    r_arcsec: float,
) -> MatchInfo:
    """
    Match catalog-2 rows to catalog-1 rows by name and sky coordinates.

    Matching is done in two steps:

    1. exact name matching,
    2. nearest-neighbor coordinate matching for rows not matched by name.

    The returned arrays all have the same length as catalog 2.
    Is something like cat 2 (in length and indices) with info where in cat 1
    one can find matched stuff.

    :param name_cat1: Object names in catalog 1.
    :type name_cat1: object
    :param ra_cat1: Right ascension values of catalog 1.
    :type ra_cat1: object
    :param dec_cat1: Declination values of catalog 1.
    :type dec_cat1: object
    :param name_cat2: Object names in catalog 2.
    :type name_cat2: object
    :param ra_cat2: Right ascension values of catalog 2.
    :type ra_cat2: object
    :param dec_cat2: Declination values of catalog 2.
    :type dec_cat2: object
    :param r_arcsec: Maximum coordinate matching distance in arcseconds.
    :type r_arcsec: float
    :returns: Dictionary containing match masks, indices, and separations.
    :rtype: dict[str, numpy.ndarray]
    """
    cat1_index_for_cat2, match_by_name = _match_cat2_to_cat1_by_name(
        name_cat1,
        name_cat2,
    )
    coord_indices, coord_mask, sep_arcsec = _match_cat2_to_cat1_by_coords(
        ra_cat1,
        dec_cat1,
        ra_cat2,
        dec_cat2,
        r_arcsec,
    )

    unmatched_by_name = ~match_by_name
    match_by_coords = unmatched_by_name & coord_mask
    cat1_index_for_cat2[match_by_coords] = coord_indices[match_by_coords]

    mask_cat2_in_cat1 = match_by_name | match_by_coords

    return {
        "mask_cat2_in_cat1": mask_cat2_in_cat1,
        "cat1_index_for_cat2": cat1_index_for_cat2,
        "match_by_name": match_by_name,
        "match_by_coords": match_by_coords,
        "sep_arcsec": sep_arcsec,
    }


def get_mask_cat2_in_cat1(
    name_cat1: Any,
    ra_cat1: Any,
    dec_cat1: Any,
    name_cat2: Any,
    ra_cat2: Any,
    dec_cat2: Any,
    r_arcsec: float,
) -> np.ndarray:
    """
    Return a mask showing which catalog-2 rows are present in catalog 1.

    Matching is done first by exact object name and then by sky coordinates.

    :param name_cat1: Object names in catalog 1.
    :type name_cat1: object
    :param ra_cat1: Right ascension values of catalog 1.
    :type ra_cat1: object
    :param dec_cat1: Declination values of catalog 1.
    :type dec_cat1: object
    :param name_cat2: Object names in catalog 2.
    :type name_cat2: object
    :param ra_cat2: Right ascension values of catalog 2.
    :type ra_cat2: object
    :param dec_cat2: Declination values of catalog 2.
    :type dec_cat2: object
    :param r_arcsec: Maximum coordinate matching distance in arcseconds.
    :type r_arcsec: float
    :returns: Boolean mask with one entry per catalog-2 row.
    :rtype: numpy.ndarray
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


def _invert_cat2_matches_to_cat1(
    match_info: MatchInfo,
    len_cat1: int,
) -> dict[str, np.ndarray]:
    """
    Convert catalog-2-based match information to catalog-1-based arrays.

    ``get_cat2_in_cat1_match_info`` stores one entry per catalog-2 row. For a
    left-join-like table, one entry per catalog-1 row is needed instead. If
    multiple catalog-2 rows match the same catalog-1 row, the first match is
    kept.
    Difficult to understand but okay, is something like cat 1 (in length and
    indices) with info where in cat 2 one can find matched stuff.

    :param match_info: Match dictionary from ``get_cat2_in_cat1_match_info``.
    :type match_info: dict[str, numpy.ndarray]
    :param len_cat1: Number of rows in catalog 1.
    :type len_cat1: int
    :returns: Dictionary with one entry per catalog-1 row.
    :rtype: dict[str, numpy.ndarray]
    """
    cat2_index_for_cat1 = np.full(len_cat1, -1, dtype=int)
    cat1_has_match = np.zeros(len_cat1, dtype=bool)
    cat1_match_by_name = np.zeros(len_cat1, dtype=bool)
    cat1_match_by_coords = np.zeros(len_cat1, dtype=bool)
    cat1_match_sep_arcsec = np.ma.masked_all(len_cat1, dtype=float)

    matched_cat2_indices = np.where(match_info["mask_cat2_in_cat1"])[0]

    for cat2_index in matched_cat2_indices:
        # get index in match info corresponding to cat 1
        cat1_index = match_info["cat1_index_for_cat2"][cat2_index]

        # exit if no match
        if cat1_index < 0:
            continue

        # only first row
        if cat1_has_match[cat1_index]:
            continue

        # inverse index
        cat2_index_for_cat1[cat1_index] = cat2_index

        # fill in rest of match info
        cat1_has_match[cat1_index] = True
        cat1_match_by_name[cat1_index] = match_info["match_by_name"][
            cat2_index
        ]
        cat1_match_by_coords[cat1_index] = match_info["match_by_coords"][
            cat2_index
        ]
        cat1_match_sep_arcsec[cat1_index] = match_info["sep_arcsec"][
            cat2_index
        ]

    return {
        "cat2_match_index": cat2_index_for_cat1,
        "cat2_match_by_name": cat1_match_by_name,
        "cat2_match_by_coords": cat1_match_by_coords,
        "cat2_match_sep_arcsec": cat1_match_sep_arcsec,
    }


def _get_output_column_name(
    column_name: str,
    existing_colnames: list[str],
    suffix: str,
) -> str:
    """
    Return an output column name that avoids name collisions.

    :param column_name: Original column name.
    :type column_name: str
    :param existing_colnames: Names already present in the output table.
    :type existing_colnames: list[str]
    :param suffix: Suffix to add if the name already exists.
    :type suffix: str
    :returns: Safe output column name.
    :rtype: str
    """
    if column_name in existing_colnames:
        return f"{column_name}{suffix}"

    return column_name


def _build_masked_cat2_columns(
    cat1: Table,
    cat2: Table,
    cat2_cols: list[str],
    cat2_index_for_cat1: np.ndarray,
    cat2_suffix: str,
) -> Table:
    """
    Build masked catalog-2 table with columns aligned to catalog-1 rows.

    :param cat1: Left-side catalog whose row count defines the output length.
    :type cat1: astropy.table.Table
    :param cat2: Right-side catalog from which values are copied.
    :type cat2: astropy.table.Table
    :param cat2_cols: Catalog-2 columns to append.
    :type cat2_cols: list[str]
    :param cat2_index_for_cat1: Catalog-2 row index for each catalog-1 row.
    :type cat2_index_for_cat1: numpy.ndarray
    :param cat2_suffix: Suffix for catalog-2 columns that collide with cat1.
    :type cat2_suffix: str
    :returns: Masked table with catalog-2 columns aligned to catalog 1.
    :rtype: astropy.table.Table
    """
    cat2_matched = Table(masked=True)
    cat1_rows_with_match = cat2_index_for_cat1 >= 0

    for col in cat2_cols:
        output_col = _get_output_column_name(
            col,
            cat1.colnames,
            cat2_suffix,
        )
        cat2_matched[output_col] = np.ma.masked_all(
            len(cat1),
            dtype=cat2[col].dtype,
        )
        cat2_matched[output_col][cat1_rows_with_match] = cat2[col][
            cat2_index_for_cat1[cat1_rows_with_match]
        ]

    return cat2_matched


def add_cat2_to_cat1_by_name_or_coords(
    cat1: Table,
    cat2: Table,
    name_cat1_col: str,
    ra_cat1_col: str,
    dec_cat1_col: str,
    name_cat2_col: str,
    ra_cat2_col: str,
    dec_cat2_col: str,
    r_arcsec: float,
    cat2_cols: list[str] | None = None,
    cat2_suffix: str = "_cat2",
    match_info: MatchInfo | None = None,
) -> Table:
    """
    Add columns from catalog 2 to catalog 1 by name or coordinate match.

    This behaves like a left join: all catalog-1 rows are preserved. Selected
    catalog-2 columns are appended where a match is found. Matching is done by
    exact name first and then by nearest sky coordinate within ``r_arcsec``.

    If ``match_info`` is provided, matching is not recomputed.

    :param cat1: Left-side catalog.
    :type cat1: astropy.table.Table
    :param cat2: Right-side catalog.
    :type cat2: astropy.table.Table
    :param name_cat1_col: Name column in catalog 1.
    :type name_cat1_col: str
    :param ra_cat1_col: Right ascension column in catalog 1.
    :type ra_cat1_col: str
    :param dec_cat1_col: Declination column in catalog 1.
    :type dec_cat1_col: str
    :param name_cat2_col: Name column in catalog 2.
    :type name_cat2_col: str
    :param ra_cat2_col: Right ascension column in catalog 2.
    :type ra_cat2_col: str
    :param dec_cat2_col: Declination column in catalog 2.
    :type dec_cat2_col: str
    :param r_arcsec: Maximum coordinate matching distance in arcseconds.
    :type r_arcsec: float
    :param cat2_cols: Catalog-2 columns to append. Appends all if ``None``.
    :type cat2_cols: list[str] or None
    :param cat2_suffix: Suffix for catalog-2 columns colliding with cat1.
    :type cat2_suffix: str
    :param match_info: Optional precomputed catalog-2-to-catalog-1 matches.
    :type match_info: dict[str, numpy.ndarray] or None
    :returns: Catalog 1 with matched catalog-2 columns appended.
    :rtype: astropy.table.Table
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

    cat1_match_info = _invert_cat2_matches_to_cat1(
        match_info,
        len(cat1),
    )
    cat2_matched = _build_masked_cat2_columns(
        cat1=cat1,
        cat2=cat2,
        cat2_cols=cat2_cols,
        cat2_index_for_cat1=cat1_match_info["cat2_match_index"],
        cat2_suffix=cat2_suffix,
    )

    cat2_matched["cat2_match_by_name"] = cat1_match_info[
        "cat2_match_by_name"
    ]
    cat2_matched["cat2_match_by_coords"] = cat1_match_info[
        "cat2_match_by_coords"
    ]
    cat2_matched["cat2_match_sep_arcsec"] = cat1_match_info[
        "cat2_match_sep_arcsec"
    ]
    cat2_matched["cat2_match_index"] = cat1_match_info["cat2_match_index"]

    return hstack([cat1, cat2_matched])


def nearest_neighbor_distances_units(ra: Any, dec: Any) -> np.ndarray:
    """
    Compute nearest-neighbor angular distances for a catalog of stars.

    The nearest neighbor is computed on the sky. For catalogs with fewer than
    two rows, no nearest neighbor exists and ``np.nan`` is returned for each
    row.

    :param ra: Right ascension values in degrees or with angular units.
    :type ra: object
    :param dec: Declination values in degrees or with angular units.
    :type dec: object
    :returns: Nearest-neighbor distances in arcseconds.
    :rtype: numpy.ndarray
    """
    stars = _make_skycoord(ra, dec)

    if len(stars) < 2:
        return np.full(len(stars), np.nan)

    # nthneighbor=2 because the closest match in the same catalog is the object itself.
    _, sep2d, _ = stars.match_to_catalog_sky(stars, nthneighbor=2)

    return sep2d.arcsec
