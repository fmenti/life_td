from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

import astropy as ap  # used for votables / tables
import numpy as np  # used for arrays
from astropy.table import Table

# Self created modules
from provider.utils import query
from utils.io import objecttostring, save


TAP_URL_HEID = "http://dc.zah.uni-heidelberg.de/tap"
TAP_URL_GVO = "http://dc.g-vo.org/tap"
TAP_URL_DEV = "http://localhost:8080/tap"


def flag_ecliptic(ang: float, ra: np.ndarray, dec: np.ndarray) -> list[str]:
    """
    Compute flags for objects close to the ecliptic.

    The ecliptic is approximated as::

        dec_ecl = 23.4 * sin(2*pi*ra/360)

    An object is flagged "True" if its declination lies within +/- ``ang``
    degrees of this curve (using RA in degrees).

    :param ang: Half-width of the band around the ecliptic, in degrees.
    :type ang: float
    :param ra: Right ascension values in degrees.
    :type ra: numpy.ndarray
    :param dec: Declination values in degrees.
    :type dec: numpy.ndarray
    :returns: List of string flags ("True"/"False") aligned with ``ra``.
    :rtype: list[str]
    """
    ecliptic = 23.4 * np.sin(2 * np.pi * ra / 360.0)
    return [
        "True" if (-ang + ecliptic[j]) < dec[j] < (ang + ecliptic[j])
        else "False"
        for j in range(len(ra))
    ]

def apply_stability_constraint(hz_stability: Table, a_max: float) -> Table:
    """
    Filter to binaries with stable  S-type planet orbits inside ``a_max``.

    The input is expected to contain (at least) the columns:
    ``parent_main_id`` and ``a_crit_s`` and to be grouped in pairs per
    ``parent_main_id`` (as produced by :func:`assign_critical_separation`).

    A pair is kept if::

        a_max < min(a_crit_s_of_pair)

    :param hz_stability: Table of candidate binaries suitable for stability
        checks.
    :type hz_stability: astropy.table.Table
    :param a_max: Maximum planet semi-major axis of interest (AU).
    :type a_max: float
    :returns: Table containing only the rows that pass the stability cut.
    :rtype: astropy.table.Table
    """
    final = hz_stability[:0].copy()  # keep only columns / metadata
    indices = hz_stability.group_by("parent_main_id").groups.indices

    for i in range(len(indices) - 1):
        left = indices[i]
        right = indices[i] + 1
        if a_max < min(hz_stability["a_crit_s"][left],
                       hz_stability["a_crit_s"][right]):
            final.add_row(hz_stability[left])
            final.add_row(hz_stability[right])

    return final

def assign_critical_separation(multiples: Table) -> tuple[Table, Table]:
    """
    Compute critical semi-major axis for stability for suitable binaries.

    This function:
    1) sorts by ``parent_main_id``,
    2) selects rows where ``suitable_companions`` is True,
    3) adds/overwrites the column ``a_crit_s`` based on :func:`crit_sep`.

    Notes:
    - The routine assumes the suitable companions are ordered in pairs
      (primary/secondary alternating). This matches the original behavior.
    - The returned first element is the (sorted) input table.

    :param multiples: Table of objects flagged as multiples.
    :type multiples: astropy.table.Table
    :returns: Tuple of (sorted multiples table, hz_stability subset table).
    :rtype: (astropy.table.Table, astropy.table.Table)
    """
    multiples.sort("parent_main_id")
    hz_stability = multiples[multiples["suitable_companions"] == True]

    # Initialize a_crit_s with same type/unit as sep_phys_value.
    hz_stability["a_crit_s"] = hz_stability["sep_phys_value"]

    for i in range(len(hz_stability)):
        m_p = hz_stability["mass_st_value"][i]
        m_s = (
            hz_stability["mass_st_value"][i + 1]
            if i % 2 == 0
            else hz_stability["mass_st_value"][i - 1]
        )
        mu = m_s / (m_p + m_s)
        hz_stability["a_crit_s"][i] = crit_sep(
            eps=0.0,
            mu=mu,
            a_bin=hz_stability["sep_phys_value"][i],
        )[0]

    return multiples, hz_stability

def crit_sep(eps: float, mu: float, a_bin: Any) -> tuple[Any, Any]:
    """
    Compute critical semi-major axes for planetary orbit stability.

    Implements the Holman & Wiegert (1999) fits for:
    - S-type (circumstellar) stability limit ``a_crit_s``
    - P-type (circumbinary) stability limit ``a_crit_p``

    ``a_bin`` is treated as a scalar that supports multiplication (often an
    Astropy Quantity, hence typed as ``Any`` here).

    :param eps: Binary orbit eccentricity.
    :type eps: float
    :param mu: Mass fraction ``m_s / (m_p + m_s)``.
    :type mu: float
    :param a_bin: Binary semi-major axis (e.g. AU quantity).
    :type a_bin: Any
    :returns: ``(a_crit_s, a_crit_p)`` stability limits.
    :rtype: tuple[Any, Any]
    """
    a_crit_s = (
        0.464
        - 0.38 * mu
        - 0.631 * eps
        + 0.586 * mu * eps
        + 0.15 * eps**2
        - 0.198 * mu * eps**2
    ) * a_bin
    a_crit_p = (
        1.6
        + 5.1 * eps
        - 2.22 * eps**2
        + 4.12 * mu
        - 4.27 * eps * mu
        - 5.09 * mu**2
        + 4.61 * eps**2 * mu**2
    ) * a_bin
    return a_crit_s, a_crit_p

def deal_with_separation(multiples: Table) -> Table:
    """
    Convert angular separations into physical separations.

    Adds/overwrites:
    - ``sep_flag``: True where ``sep_ang_value`` is not masked
    - ``sep_phys_value``: ``sep_ang_value * dist_st_value`` rounded to 0.1 AU

    Units are set to AU.

    :param multiples: Multiples table with ``sep_ang_value`` and
        ``dist_st_value``.
    :type multiples: astropy.table.Table
    :returns: The updated table (same object, modified in-place).
    :rtype: astropy.table.Table
    """
    multiples["sep_flag"] = np.invert(multiples["sep_ang_value"].mask)

    # Initiate new column with same properties.
    multiples["sep_phys_value"] = multiples["sep_ang_value"]
    multiples["sep_phys_value"].unit = ap.units.AU

    for i in range(len(multiples)):
        if multiples["sep_flag"][i] == True:
            multiples["sep_phys_value"][i] = np.round(
                multiples["sep_ang_value"][i] * multiples["dist_st_value"][i], 1
            )

    return multiples

def sorting_number_of_id(
    input_column: Iterable[Any],
    occurences: int,
    match_column: Iterable[Any],
) -> np.ndarray:
    """
    Flag rows whose IDs appear a specific number of times.

    This is used to detect "exactly N occurrences" patterns, e.g.:
    - parent ids that occur exactly twice (binary pairs)
    - ids that occur exactly once (unique parent)

    :param input_column: Sequence of ids to count occurrences in.
    :type input_column: Iterable[Any]
    :param occurences: Required count.
    :type occurences: int
    :param match_column: Sequence to flag; returned boolean array aligns with
        this sequence.
    :type match_column: Iterable[Any]
    :returns: Boolean NumPy array: True where element in ``match_column`` is in
        the subset of ids that occur exactly ``occurences`` times.
    :rtype: numpy.ndarray
    """
    unique_id, number_of_repetitions = np.unique(input_column,
                                                 return_counts=True)
    subset = unique_id[number_of_repetitions == occurences]
    return np.isin(match_column, subset)

def flag_hz_orbit_stability(multiples: Table) -> Table:
    """
    Flag multiples for habitable-zone orbit stability (simple binary model).

    Pipeline (behavior preserved):
    - compute physical separations
    - build a ``requirement_flag`` mask based on multiple boolean criteria
    - mark ``suitable_companions`` where each parent has exactly 2 suitable rows
    - compute critical separations ``a_crit_s``
    - keep only those pairs stable for ``a_max=10.0`` AU
    - set ``stableHZ`` to "True"/"False" based on membership in the final set

    :param multiples: Multiples table.
    :type multiples: astropy.table.Table
    :returns: Updated table with added/overwritten columns.
    :rtype: astropy.table.Table
    """
    multiples = deal_with_separation(multiples)

    multiples["requirement_flag"] = (
        multiples["sep_flag"]
        & (multiples["ms_temp_class"] == "True")
        & (multiples["ms_lum_class"] == "True")
        & multiples["mass_flag"]
        & multiples["trivial_binaries"]
    )

    multiples["suitable_companions"] = sorting_number_of_id(
        multiples["parent_main_id"][multiples["requirement_flag"]],
        2,
        multiples["parent_main_id"],
    )

    multiples, hz_stability = assign_critical_separation(multiples)

    final = apply_stability_constraint(hz_stability, a_max=10.0)

    multiples["stableHZ"] = np.where(
        np.isin(multiples["main_id"], final["main_id"]), "True", "False"
    )

    return multiples


def flag_trivial_binaries(
    catalog: Table,
    children: Table
) -> tuple[Table, Table]:
    """
    Split catalog into singles and multiples, and flag "trivial" binaries.

    A "trivial" binary here means:
    - not a higher-order multiple (parent is not also a child)
    - each object has a single parent
    - each parent has exactly 2 entries in the multiples table

    :param catalog: Input catalog table containing a ``binary_flag`` column.
    :type catalog: astropy.table.Table
    :param children: Child table with at least ``child_main_id``.
    :type children: astropy.table.Table
    :returns: ``(singles, multiples)`` tables.
    :rtype: (astropy.table.Table, astropy.table.Table)
    """
    singles = catalog[np.where(catalog["binary_flag"] == "False")]
    multiples = catalog[np.where(catalog["binary_flag"] == "True")]

    # Parent object is also a child -> higher order multiple.
    multiples["higher_order_multiples"] = np.isin(
        multiples["parent_main_id"], children["child_main_id"]
    )

    # Objects that have multiple parents.
    multiples["single_parent"] = sorting_number_of_id(
        multiples["main_id"], 1, multiples["main_id"]
    )

    # Parents with exactly two entries in the multiples table.
    binaries_in_multiples_table = sorting_number_of_id(
        multiples["parent_main_id"], 2, multiples["parent_main_id"]
    )

    multiples["trivial_binaries"] = (
        (multiples["higher_order_multiples"] == False)
        & multiples["single_parent"]
        & binaries_in_multiples_table
    )
    return singles, multiples


def flag_non_main_sequence_stars(catalog: Table) -> Table:
    """
    Add simple main-sequence classification flags.

    Adds:
    - ``ms_temp_class``: "True" if ``class_temp`` in OBAFGKM
    - ``ms_lum_class``: "True" if ``class_lum`` is V
    - ``mass_flag``: True where ``mass_st_value`` is not masked

    :param catalog: Input table with ``class_temp``, ``class_lum``,
        ``mass_st_value``.
    :type catalog: astropy.table.Table
    :returns: Updated table.
    :rtype: astropy.table.Table
    """
    ms_tempclass = np.array(["O", "B", "A", "F", "G", "K", "M"])
    catalog["ms_temp_class"] = np.where(
        np.isin(catalog["class_temp"], ms_tempclass), "True", "False"
    )

    ms_lumclass = np.array(["V"])
    catalog["ms_lum_class"] = np.where(
        np.isin(catalog["class_lum"], ms_lumclass), "True", "False"
    )

    catalog["mass_flag"] = np.invert(catalog["mass_st_value"].mask)
    return catalog

def add_unresolved_binaries(systems: Table,
                            children: Table,
                            stars: Table
                            ) -> Table:
    """
    Add unresolved binary systems to the star table.

    A system is considered "unresolved" if it has no child entry after joining
    with the children table (left join on system object_id -> parent_object_idref).

    The result is a vstack of:
    - ``stars`` with ``unresolved_binaries="False"``
    - unresolved systems with ``unresolved_binaries="True"`` and a few columns
      removed.

    :param systems: Systems table from TAP (must include ``object_id``).
    :type systems: astropy.table.Table
    :param children: Children table (must include ``parent_object_idref`` and
        ``child_main_id``).
    :type children: astropy.table.Table
    :param stars: Stars table to extend.
    :type stars: astropy.table.Table
    :returns: Combined catalog table.
    :rtype: astropy.table.Table
    """
    systems_with_child_info = ap.table.join(
        systems,
        children,
        keys_left="object_id",
        keys_right="parent_object_idref",
        join_type="left",
    )

    unresolved_binaries = systems_with_child_info[
        systems_with_child_info["child_main_id"].mask
    ]

    stars["unresolved_binaries"] = ["False" for _ in range(len(stars))]
    unresolved_binaries["unresolved_binaries"] = [
        "True" for _ in range(len(unresolved_binaries))
    ]

    unresolved_binaries.remove_columns(
        ["object_id", "child_main_id", "child_type", "parent_object_idref"]
    )
    return ap.table.vstack([stars, unresolved_binaries])


def _query_star_like(
    *,
    service: str,
    distance_cut: float,
    object_type: str,
    include_object_id: bool,
) -> Table:
    """
    Query star_basic joined to object/h_link for a specific object type.

    This consolidates the shared SQL between :func:`query_stars` and
    :func:`query_systems` while preserving the returned columns.

    :param service: TAP base URL.
    :type service: str
    :param distance_cut: Maximum distance in pc.
    :type distance_cut: float
    :param object_type: Object type filter ('st' or 'sy').
    :type object_type: str
    :param include_object_id: If True, include ``o.object_id`` in SELECT.
    :type include_object_id: bool
    :returns: Query result table.
    :rtype: astropy.table.Table
    """
    object_id_select = "o.object_id,\n        " if include_object_id else ""
    adql_query = f"""
    SELECT {object_id_select}o.main_id, sb.coo_ra, sb.coo_dec, sb.sptype_string,
        sb.plx_value, sb.dist_st_value, sb.coo_gal_l, sb.coo_gal_b,
        sb.teff_st_value, teff_source.ref AS teff_ref,
        sb.mass_st_value, mass_source.ref AS mass_ref,
        sb.radius_st_value, radius_source.ref AS radius_ref,
        sb.binary_flag, binary_source.ref AS binary_ref,
        sb.mag_i_value, sb.mag_j_value, sb.class_lum, sb.class_temp,
        o_parent.main_id AS parent_main_id, sb_parent.sep_ang_value
    FROM life_td.star_basic AS sb
    JOIN life_td.object AS o ON sb.object_idref=o.object_id
    LEFT JOIN life_td.h_link AS h ON o.object_id=h.child_object_idref
    LEFT JOIN life_td.object AS o_parent ON
        h.parent_object_idref=o_parent.object_id
    LEFT JOIN life_td.star_basic AS sb_parent ON
        o_parent.object_id=sb_parent.object_idref
    LEFT JOIN life_td.source AS radius_source ON
        sb.radius_st_source_idref=radius_source.source_id
    LEFT JOIN life_td.source AS mass_source ON
        sb.mass_st_source_idref=mass_source.source_id
    LEFT JOIN life_td.source AS teff_source ON
        sb.teff_st_source_idref=teff_source.source_id
    LEFT JOIN life_td.source AS binary_source ON
        sb.binary_source_idref=binary_source.source_id
    WHERE o.type = '{object_type}' AND sb.dist_st_value < {distance_cut}
    """
    return query(service, adql_query)


def query_systems(service: str, distance_cut: float) -> Table:
    """
    Query all systems within a distance cut.

    Includes the ``object_id`` column (needed later for joining on children).

    :param service: TAP base URL.
    :type service: str
    :param distance_cut: Maximum distance in pc.
    :type distance_cut: float
    :returns: Systems table.
    :rtype: astropy.table.Table
    """
    return _query_star_like(
        service=service,
        distance_cut=distance_cut,
        object_type="sy",
        include_object_id=True,
    )


def query_children(service: str) -> Table:
    """
    Query child objects (excluding planets and disks) with their parent ids.

    :param service: TAP base URL.
    :type service: str
    :returns: Table with child_main_id, child_type, parent_object_idref.
    :rtype: astropy.table.Table
    """
    adql_query = """
    SELECT o.main_id as child_main_id, o.type as child_type,
           h.parent_object_idref
    FROM life_td.h_link AS h
    JOIN life_td.object AS o on o.object_id=h.child_object_idref
    WHERE o.type not in ('pl','di')
    """
    return query(service, adql_query)

def query_stars(service: str, distance_cut: float) -> Table:
    """
    Query all stars within a distance cut.

    :param service: TAP base URL.
    :type service: str
    :param distance_cut: Maximum distance in pc.
    :type distance_cut: float
    :returns: Stars table.
    :rtype: astropy.table.Table
    """
    return _query_star_like(
        service=service,
        distance_cut=distance_cut,
        object_type="st",
        include_object_id=False,
    )

def choose_service(service: str) -> str:
    """
    Map a short service name to a TAP URL.

    Supported values:
    - ``"heid"`` -> Heidelberg service
    - ``"gvo"`` -> GAVO service
    - anything else -> local development service

    :param service: Short service selector.
    :type service: str
    :returns: TAP base URL.
    :rtype: str
    """
    if service == "heid":
        return TAP_URL_HEID
    if service == "gvo":
        return TAP_URL_GVO
    return TAP_URL_DEV



def save_catalog(starcat5: Table) -> None:
    """
    Persist StarCat5 both as ECSV and via the project save helper.

    The output location is::

        <repo_root>/additional_data/catalogs/StarCat5.ecsv

    ``starcat5`` is converted via :func:`utils.io.objecttostring` first, to
    ensure variable-length strings/objects are serialized consistently.

    :param starcat5: Final catalog table to save.
    :type starcat5: astropy.table.Table
    :returns: None
    :rtype: None
    """
    starcat5 = objecttostring(starcat5)

    project_root = Path(__file__).resolve().parents[3]  # .../life_td
    catalogs_dir = project_root / "additional_data" / "catalogs"
    catalogs_dir.mkdir(parents=True, exist_ok=True)

    out_path = catalogs_dir / "StarCat5.ecsv"
    starcat5.write(str(out_path), delimiter=",", overwrite=True)

    save([starcat5], ["StarCat5"], location=str(catalogs_dir) + "/")


def main(distance_cut = 30.0, service_type = "") -> int:
    """
    Build and save the StarCat5 catalog.

    Steps (high level):
    - choose TAP service
    - query stars, systems, and children
    - add unresolved binaries
    - flag main-sequence membership and trivial binaries
    - flag HZ stability for multiples
    - stack singles + multiples into StarCat5
    - add an ecliptic proximity flag
    - save results

    :returns: Exit code (0 for success).
    :rtype: int
    """
    service = choose_service(service_type)

    queried_stars = query_stars(service, distance_cut)
    queried_children = query_children(service)
    queried_systems = query_systems(service, distance_cut)

    stars_with_ub = add_unresolved_binaries(
        queried_systems,
        queried_children,
        queried_stars,
    )

    # Process the result.
    flagged_non_ms = flag_non_main_sequence_stars(stars_with_ub)
    singles, multiples = flag_trivial_binaries(flagged_non_ms, queried_children)
    multiples = flag_hz_orbit_stability(multiples)

    starcat5 = ap.table.vstack([singles, multiples])

    angle = 45
    starcat5[f"ecliptic_pm{angle}deg"] = flag_ecliptic(
        angle, starcat5["coo_ra"], starcat5["coo_dec"]
    )

    save_catalog(starcat5)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



