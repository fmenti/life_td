from astropy.table import Table
from provider.utils import query

# db needs to be ingested by GAVO DaCHS and hosted on localhost for this to be meaningful


def test_life_db_filter_objects_by_type():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT TOP 10 object_id, main_id FROM life_td.object
                WHERE type='st'"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert len(life_db) == 10
    assert life_db.colnames == ["object_id", "main_id"]


def test_life_db_all_children_of_an_object():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT DISTINCT main_id as Child_main_id, object_id as
                child_object_id
                FROM life_td.h_link
                JOIN life_td.ident as p on p.object_idref=parent_object_idref
                JOIN life_td.object on object_id=child_object_idref
                WHERE p.id = '* alf Cen'"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert len(life_db) == 3
    assert life_db.colnames == ["child_main_id", "child_object_id"]


def test_life_db_all_parents_of_an_object():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT DISTINCT main_id as parent_main_id, object_id as 
                parent_object_id
                FROM life_td.h_link
                JOIN life_td.ident as c on c.object_idref=child_object_idref
                JOIN life_td.object on object_id=parent_object_idref
                WHERE c.id =  '* alf Cen A'"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert len(life_db) == 1
    assert life_db.colnames == ["parent_main_id", "parent_object_id"]


def test_life_db_all_specific_measurements_of_an_object():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT *
                FROM life_td.mes_teff_st
                JOIN life_td.ident USING(object_idref)
                WHERE id = '* alf Cen A'"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert len(life_db) == 1
    assert life_db.colnames == [
        "object_idref",
        "teff_st_value",
        "teff_st_err",
        "teff_st_qual",
        "teff_st_source_idref",
        "id",
        "id_source_idref",
    ]


def test_life_db_all_basic_stellar_data_from_an_object():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT  *
            FROM life_td.star_basic
            JOIN life_td.ident USING(object_idref)
            WHERE id = '* alf Cen'"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert len(life_db) == 1


def test_life_db_all_basic_disk_data_from_host_name():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT DISTINCT main_id disk_main_id, object_id as 
                disk_object_id, db.*
                FROM life_td.h_link
                JOIN life_td.disk_basic as db ON 
                 db.object_idref=child_object_idref
                JOIN life_td.ident as p on p.object_idref=parent_object_idref
                JOIN life_td.object on object_id=child_object_idref
                WHERE p.id = '* bet Pic' and type='di'"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert len(life_db) == 2
    assert life_db.colnames == [
        "disk_main_id",
        "disk_object_id",
        "object_idref",
        "rad_value",
        "rad_err",
        "rad_qual",
        "rad_rel",
        "rad_source_idref",
    ]


def test_life_db_missing_reliable_measurements():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT star_ob.main_id as star_name, plx_value, plx_err,
                plx_qual, plx_source_idref
                FROM life_td.star_basic as s
                JOIN life_td.object as star_ob on
                (s.object_idref=star_ob.object_id)
                WHERE plx_value is Null or plx_qual in ('D','E') or
                plx_qual is Null"""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert life_db.colnames == [
        "star_name",
        "plx_value",
        "plx_err",
        "plx_qual",
        "plx_source_idref",
    ]


def test_life_db_LIFE_StarCat_candidates():
    life_db_provider = Table()
    life_db_provider["provider_url"] = ["http://localhost:8080/tap"]
    adql_query = [
        """SELECT o.main_id, sb.coo_ra, sb.coo_dec, sb.plx_value, 
                 sb.dist_st_value, sb.sptype_string, sb.coo_gal_l, 
                 sb.coo_gal_b, sb.teff_st_value, sb.mass_st_value, 
                 sb.radius_st_value, sb.binary_flag, sb.mag_i_value, 
                 sb.mag_j_value,  sb.class_lum, sb.class_temp, 
                 o_parent.main_id AS parent_main_id, 
                 sb_parent.sep_ang_value
                FROM life_td.star_basic AS sb
                JOIN life_td.object AS o ON sb.object_idref=o.object_id
                LEFT JOIN life_td.h_link AS h ON 
                 o.object_id=h.child_object_idref
                LEFT JOIN life_td.object AS o_parent ON 
                 h.parent_object_idref=o_parent.object_id
                LEFT JOIN life_td.star_basic AS sb_parent ON 
                 o_parent.object_id=sb_parent.object_idref
                WHERE o.type = 'st' AND sb.dist_st_value < 30."""
    ]
    life_db = query(life_db_provider["provider_url"][0], adql_query[0])

    assert type(life_db) == type(Table())
    assert life_db.colnames == [
        "main_id",
        "coo_ra",
        "coo_dec",
        "plx_value",
        "dist_st_value",
        "sptype_string",
        "coo_gal_l",
        "coo_gal_b",
        "teff_st_value",
        "mass_st_value",
        "radius_st_value",
        "binary_flag",
        "mag_i_value",
        "mag_j_value",
        "class_lum",
        "class_temp",
        "parent_main_id",
        "sep_ang_value",
    ]
