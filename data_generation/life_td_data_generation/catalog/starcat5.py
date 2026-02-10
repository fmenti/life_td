import astropy as ap  # Used for votables
import numpy as np  # Used for arrays
import pyvo as vo  # Used for catalog query

# Self created modules
from provider.utils import query
from utils.io import save, objecttostring

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
    catalog = flag_non_main_sequence_stars(stars_with_ub)


