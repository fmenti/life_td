import astropy as ap  # Used for votables
import numpy as np  # Used for arrays
import pyvo as vo  # Used for catalog query

# Self created modules
from provider.utils import query
from utils.io import save, objecttostring

def apply_stability_constraint(HZstability,a_max):
    final = HZstability[:0].copy()  # only columns
    # wait, didn't I already define this? -> was before removing some
    ind = HZstability.group_by("parent_main_id").groups.indices

    for i in range(len(ind) - 1):
        if a_max < min(HZstability["a_crit_s"][ind[i]],
                       HZstability["a_crit_s"][ind[i] + 1]):
            final.add_row(HZstability[ind[i]])
            final.add_row(HZstability[ind[i] + 1])
    return final

def assign_critical_separation(multiples):
    multiples.sort("parent_main_id")
    HZstability = multiples[multiples["suitable_companions"] == True]

    # initializing column a... like sep...
    HZstability["a_crit_s"] = HZstability["sep_phys_value"]

    for i in range(len(HZstability)):
        m_p = HZstability["mass_st_value"][i]
        if i % 2 == 0:
            m_s = HZstability["mass_st_value"][i + 1]
        else:
            m_s = HZstability["mass_st_value"][i - 1]
        mu = m_s / (m_p + m_s)
        HZstability["a_crit_s"][i] = \
            crit_sep(0, mu, HZstability["sep_phys_value"][i])[0]
        # assumed circular orbit and sep_phys = a_bin
    return multiples, HZstability

def crit_sep(eps, mu, a_bin):
    """
    Computes critical semimajor-axis for planet orbit stability.

    For binary system as described in Holman and Wiegert 1999.

    :param eps: Binary orbit excentricity.
    :type eps:
    :param mu: mass fraction with mu=m_s/(m_p+m_s), with m_s the mass
        of the star considered as perturbing binary companion and m_p
        the mass of the star the planet is orbiting.
    :type mu:
    :param a_bin: semimajor-axis of the binary stars.
    :type a_bin:
    :returns: Critical separation beyond which a planet on a S-type
        orbit (circumstellar) and on a P-type orbit (circumbinary) is
        not stable any more.
    :rtype:
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

def deal_with_separation(multiples):
    multiples["sep_flag"] = np.invert(multiples["sep_ang_value"].mask)

    # just initiating new column with same properties
    multiples["sep_phys_value"] = multiples["sep_ang_value"]

    multiples["sep_phys_value"].unit = ap.units.AU
    for i in range(len(multiples)):
        if multiples["sep_flag"][i] == True:
            multiples["sep_phys_value"][i] = np.round(
                multiples["sep_ang_value"][i] * multiples["dist_st_value"][
                    i], 1
            )
    return multiples

def sorting_number_of_id(input_column,occurences,match_column):
    """

    :param input_column:
    :param occurences:
    :param match_column: Column to with the return flag_array will belong
    """

    unique_id,number_of_repetitions=np.unique(
        input_column,return_counts=True)
    subset=unique_id[number_of_repetitions==occurences]
    flag_array=np.isin(match_column,subset)
    return flag_array

def flag_hz_orbit_stability(multiples):
    # do I need to refacture this function?

    #there is some print statement in this that I want to get rid off, look at it in jupyter notebook

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
        2, multiples["parent_main_id"])

    multiples,HZstability = assign_critical_separation(multiples)

    final = apply_stability_constraint(HZstability,a_max=10.0)

    # final["stableHZ"]
    # need to rewrite into flag not cut
    # do I use catalog or multiples?
    multiples["stableHZ"] = np.where(
        np.isin(multiples["main_id"], final["main_id"]), "True", "False")

    return multiples


def flag_trivial_binaries(catalog,children):
    # do I need to refacture this function?
    singles = catalog[np.where(catalog["binary_flag"] == "False")]
    multiples = catalog[np.where(catalog["binary_flag"] == "True")]

    multiples["higher_order_multiples"] = np.isin(
        multiples["parent_main_id"], children["child_main_id"]
    )

    multiples["single_parent"] = sorting_number_of_id(
        multiples["main_id"], 1, multiples["main_id"])

    # flag stuff like *  16 Cyg A , B, and C who all have parent *  16 Cyg
    binaries_in_multiples_table = sorting_number_of_id(
        multiples["parent_main_id"], 2, multiples["parent_main_id"])

    multiples["trivial_binaries"] = (
        (multiples["higher_order_multiples"] == False)
        & multiples["single_parent"]
        & binaries_in_multiples_table
    )
    return singles,multiples


def flag_non_main_sequence_stars(catalog):
    ms_tempclass = np.array(["O", "B", "A", "F", "G", "K", "M"])

    catalog["ms_temp_class"] = np.where(
        np.isin(catalog["class_temp"], ms_tempclass), "True", "False")

    ms_lumclass = np.array(["V"])

    catalog["ms_lum_class"] = np.where(
        np.isin(catalog["class_lum"], ms_lumclass), "True", "False")

    catalog["mass_flag"] = np.invert(catalog["mass_st_value"].mask)

    return catalog

def add_unresolved_binaries(systems,children,stars):
    systems_with_child_info = ap.table.join(systems, children,
                                            keys_left="object_id",
                                            keys_right="parent_object_idref",
                                            join_type="left")
    unresolved_binaries = systems_with_child_info[
        systems_with_child_info["child_main_id"].mask]
    stars["unresolved_binaries"] = ["False" for _ in
                                            range(len(stars))]
    unresolved_binaries["unresolved_binaries"] = ["True" for _ in range(
        len(unresolved_binaries))]
    unresolved_binaries.remove_columns(["object_id",
                                        "child_main_id",
                                        "child_type",
                                        "parent_object_idref"])
    catalog = ap.table.vstack([stars, unresolved_binaries])
    return catalog


def query_systems(service, distance_cut):
    adql_query = """
    SELECT o.object_id,
        o.main_id, sb.coo_ra, sb.coo_dec, sb.sptype_string,
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
    WHERE o.type = 'sy' AND sb.dist_st_value < """ + str(distance_cut)

    return query(service,adql_query)


def query_children(service):
    adql_query = """
    SELECT o.main_id as child_main_id, o.type as child_type, h.parent_object_idref
    FROM life_td.h_link AS h
    JOIN life_td.object AS o on o.object_id=h.child_object_idref
    WHERE o.type not in ('pl','di')
    """
    return query(service,adql_query)

def query_stars(service, distance_cut):
    adql_query = """
                 SELECT o.main_id, \
                        sb.coo_ra, \
                        sb.coo_dec, \
                        sb.sptype_string,
                        sb.plx_value, \
                        sb.dist_st_value, \
                        sb.coo_gal_l, \
                        sb.coo_gal_b,
                        sb.teff_st_value, \
                        teff_source.ref   AS teff_ref,
                        sb.mass_st_value, \
                        mass_source.ref   AS mass_ref,
                        sb.radius_st_value, \
                        radius_source.ref AS radius_ref,
                        sb.binary_flag, \
                        binary_source.ref AS binary_ref,
                        sb.mag_i_value, \
                        sb.mag_j_value, \
                        sb.class_lum, \
                        sb.class_temp,
                        o_parent.main_id  AS parent_main_id, \
                        sb_parent.sep_ang_value
                 FROM life_td.star_basic AS sb
                          JOIN life_td.object AS o \
                               ON sb.object_idref = o.object_id
                          LEFT JOIN life_td.h_link AS h \
                                    ON o.object_id = h.child_object_idref
                          LEFT JOIN life_td.object AS o_parent ON
                     h.parent_object_idref = o_parent.object_id
                          LEFT JOIN life_td.star_basic AS sb_parent ON
                     o_parent.object_id = sb_parent.object_idref
                          LEFT JOIN life_td.source AS radius_source ON
                     sb.radius_st_source_idref = radius_source.source_id
                          LEFT JOIN life_td.source AS mass_source ON
                     sb.mass_st_source_idref = mass_source.source_id
                          LEFT JOIN life_td.source AS teff_source ON
                     sb.teff_st_source_idref = teff_source.source_id
                          LEFT JOIN life_td.source AS binary_source ON
                     sb.binary_source_idref = binary_source.source_id
                 WHERE o.type = 'st' \
                   AND sb.dist_st_value < """ + str(distance_cut)
    return query(service,adql_query)

def choose_service(service):
    if service == "heid":
        return "http://dc.zah.uni-heidelberg.de/tap"
    elif service == "gvo":
        return "http://dc.g-vo.org/tap"
    else:
        #development_service
        return "http://localhost:8080/tap"

if __name__ == "__main__":
    service = choose_service("")

    queried_stars = query_stars(service, 30.0)
    queried_children = query_children(service)
    queried_systems = query_systems(service, 30.0)

    stars_with_ub = add_unresolved_binaries(queried_systems,
                                      queried_children,
                                      queried_stars)

    #process the result
    flag_non_ms = flag_non_main_sequence_stars(stars_with_ub)
    singles,multiples = flag_trivial_binaries(flag_non_ms,
                                                    queried_children)
    multiples = flag_hz_orbit_stability(multiples)

    StarCat5 = ap.table.vstack([singles, multiples])

    #ecliptic
    #plots



