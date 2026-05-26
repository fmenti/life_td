from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

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

