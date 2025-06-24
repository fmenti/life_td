from utils.io import save, Path
from provider.utils import query
import numpy as np  # Used for arrays
import astropy as ap  # Used for votables

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
    a_crit_s = (0.464 - 0.38 * mu - 0.631 * eps + 0.586 * mu * eps + 0.15 * eps ** 2 \
                - 0.198 * mu * eps ** 2) * a_bin
    a_crit_p = (1.6 + 5.1 * eps - 2.22 * eps ** 2 + 4.12 * mu - 4.27 * eps * mu - 5.09 * mu ** 2 \
                + 4.61 * eps ** 2 * mu ** 2) * a_bin
    return a_crit_s, a_crit_p


def ecliptic(ang, ra, dec):
    """
    Computes if position is within angle from the ecliptic.

    :param ang: Angle in degrees.
    :type ang:
    :param ra: Right ascention in degrees.
    :type ra: np.array
    :param dec: Array of declination in degrees.
    :type dec: np.array
    :returns: Flags.
    :rtype: np.array
    """
    ecliptic = (23.4) * np.sin(2 * np.pi * ra / 360)
    flag = ['True' if dec[j] > -ang + ecliptic[j] and dec[j] < ang + ecliptic[j] \
                else 'False' for j in range(len(ra))]
    return flag


def starcat_creation(distance_cut):
    """
    LIFE-StarCat4 creation

    :param distance_cut:
    :type distance_cut:
    :param
    :type
    :param
    :type
    :param
    :type
    :returns:
    :rtype:
    """

    # version 3 parameters were: ra, dec, plx, distance, name, sptype,
    # coo_gal_l, coo_gal_b, Teff, R, M, sep_phys, binary_flag, mag_i,
    # mag_j
    adql_query = """ \
                 SELECT o.main_id, \
                        sb.coo_ra, \
                        sb.coo_dec, \
                        sb.plx_value, \
                        sb.dist_st_value, \
                        sb.sptype_string, \
                        sb.coo_gal_l, \
                        sb.coo_gal_b, \
                        sb.teff_st_value, \
                        sb.mass_st_value, \
                        sb.radius_st_value, \
                        sb.binary_flag, \
                        sb.mag_i_value, \
                        sb.mag_j_value, \
                        sb.class_lum, \
                        sb.class_temp, \
                        o_parent.main_id AS parent_main_id, \
                        sb_parent.sep_ang_value \
                 FROM life_td.star_basic AS sb \
                          JOIN life_td.object AS o ON sb.object_idref = o.object_id \
                          LEFT JOIN life_td.h_link AS h ON o.object_id = h.child_object_idref \
                          LEFT JOIN life_td.object AS o_parent ON \
                     h.parent_object_idref = o_parent.object_id \
                          LEFT JOIN life_td.star_basic AS sb_parent ON \
                     o_parent.object_id = sb_parent.object_idref \
                 WHERE o.type = 'st' \
                   AND sb.dist_st_value < """ + str(distance_cut)
    # we are only interested in object type stars, up to a distance cut
    # and well defined luminocity class (to sort out objects not around
    # main sequence)
    service = 'http://localhost:8080/tap'

    catalog = query(service, adql_query)

    ms_tempclass = np.array(['O', 'B', 'A', 'F', 'G', 'K', 'M'])
    cat_ms_tempclass = catalog[np.where(np.in1d(catalog['class_temp'], ms_tempclass))]

    ms_lumclass = np.array(['V'])
    cat_ms_lumclass = cat_ms_tempclass[np.where(np.in1d(cat_ms_tempclass['class_lum'], ms_lumclass))]

    cat_ms_lumclass.remove_rows(cat_ms_lumclass['mass_st_value'].mask.nonzero()[0])

    singles = cat_ms_lumclass[np.where(cat_ms_lumclass['binary_flag'] == 'False')]
    multiples = cat_ms_lumclass[np.where(cat_ms_lumclass['binary_flag'] == 'True')]

    adql_query2 = """ \
                  SELECT o.main_id as child_main_id, o.object_id \
                  FROM life_td.object AS o \
                           JOIN life_td.h_link AS h on o.object_id = h.child_object_idref \
                """
    h_link = query(service, adql_query2)

    higher_order_multiples = np.in1d(multiples['parent_main_id'], \
                                     h_link['child_main_id'])
    multiples.remove_rows(higher_order_multiples)


    multi_parent = []
    grouped = multiples.group_by('main_id')
    ind = grouped.groups.indices
    for i in range(len(ind) - 1):
        if ind[i + 1] - ind[i] != 1:
            multi_parent.append(grouped['main_id'][ind[i]])

    single_parent_multiples = grouped[np.where(np.invert(np.in1d(grouped['main_id'], multi_parent)))]


    sep_multiples = single_parent_multiples[np.where(single_parent_multiples['sep_ang_value'].mask == False)].copy()

    sep_multiples['sep_phys_value'] = sep_multiples['sep_ang_value']  # just initiating new column with same properties
    sep_multiples['sep_phys_value'].unit = ap.units.AU
    for i in range(len(sep_multiples)):
        sep_multiples['sep_phys_value'][i] = np.round(
            sep_multiples['sep_ang_value'][i] * sep_multiples['dist_st_value'][i], 1)

    grouped_multiples = sep_multiples.group_by('parent_main_id')
    ind = grouped_multiples.groups.indices

    result = grouped_multiples[:0].copy()

    for i in range(len(ind) - 1):
        l = ind[i + 1] - ind[i]
        if l == 2:
            result.add_row(grouped_multiples[ind[i]])
            result.add_row(grouped_multiples[ind[i] + 1])


    result['a_crit_s'] = result['sep_phys_value']  # initializing column a... like sep...

    for i in range(len(result)):
        m_p = result['mass_st_value'][i]
        if i % 2 == 0:
            m_s = result['mass_st_value'][i + 1]
        else:
            m_s = result['mass_st_value'][i - 1]
        mu = m_s / (m_p + m_s)
        result['a_crit_s'][i] = crit_sep(0, mu, result['sep_phys_value'][i])[0]
        # assumed circular orbit and sep_phys = a_bin

    final = result[:0].copy()
    # wait, didn't I already define this? -> was before removing some
    ind = result.group_by('parent_main_id').groups.indices
    a_max = 10.

    for i in range(len(ind) - 1):
        if a_max < min(result['a_crit_s'][ind[i]], result['a_crit_s'][ind[i] + 1]):
            final.add_row(result[ind[i]])
            final.add_row(result[ind[i] + 1])


    StarCat4 = ap.table.vstack([singles, final])

    # flag any object whose declination is contained within the region
    # between -(23.4+45)*sin(RA) and +(23.4+45)*sin(RA) with the
    # object's RA in degrees.
    StarCat4['ecliptic_pm45deg'] = ecliptic(45, StarCat4['coo_ra'], \
                                            StarCat4['coo_dec'])

    save([StarCat4], ['integration_test_StarCat4'], location=Path().additional_data+"/catalogs")

    return StarCat4