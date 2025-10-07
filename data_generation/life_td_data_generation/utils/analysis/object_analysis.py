import numpy as np
from astropy.io.ascii import read
from provider.utils import IdentifierCreator, fetch_main_id
from utils.io import Path
from utils.starcat4 import starcat_creation


def detail_info_object(database_tables, starname):
    # find star in db
    table = database_tables["ident"]
    main_id = table["main_id"][np.where(table["id"] == starname)][0]
    o_idref = table["object_idref"][np.where(table["id"] == starname)][0]
    print("\nIdentifier ")
    print(table[np.where(table["id"] == starname)])

    source_ref = table["id_source_idref"][np.where(table["id"] == starname)]
    source = database_tables["sources"]
    # print(source[np.where(source['source_id']==source_ref[0])])
    # print(source[np.where(source['source_id']==source_ref[1])])
    objects = database_tables["objects"]
    print("\nObject type")
    print(objects["type"][np.where(objects["main_id"] == main_id)])
    # find stars properties
    star_basic = database_tables["star_basic"]
    cols = ["main_id", "class_lum", "class_temp", "class_temp_nr"]
    print("\n Basic stellar data")
    print(star_basic[cols][np.where(star_basic["main_id"] == main_id)])

    # find parent and sibling
    print("\nParents")
    bhlink = database_tables["best_h_link"]
    print(bhlink[np.where(bhlink["child_object_idref"] == o_idref)])
    print("\nSiblings")
    parent_o_idref = bhlink["parent_object_idref"][
        np.where(bhlink["child_object_idref"] == o_idref)
    ][0]
    print(bhlink[np.where(bhlink["parent_object_idref"] == parent_o_idref)])

    print("\nSibling's properties")
    print(
        star_basic[cols][np.where(star_basic["main_id"] == "WDS J10114+4927B")]
    )


def starcat_versions_lost_stars():
    ltc3_g_stars, ltc3_k_stars = get_g_and_k_ltc3_stars()
    # update their main id
    ltc3_g_stars_update = fetch_main_id(
        ltc3_g_stars,
        id_creator=IdentifierCreator(name="main_id", colname="sim_name"),
    )
    # currently not running, but simbad query isn't either so maybe issue with internet

    ltc4 = starcat_creation(30)

    lost_g_stars = ltc3_g_stars_update["main_id"][
        np.where(
            np.invert(np.isin(ltc3_g_stars_update["main_id"], ltc4["main_id"]))
        )
    ]
    return lost_g_stars


def get_g_and_k_ltc3_stars():
    starcat3 = read(Path().additional_data + "catalogs/LIFE-StarCat3.csv")
    starcat3["class_temp"] = starcat3["sim_sptype"]
    starcat3["class_temp"] = starcat3["class_temp"].astype("S1")

    ltc3_g_stars = starcat3[np.where(starcat3["class_temp"] == "G")]
    ltc3_k_stars = starcat3[np.where(starcat3["class_temp"] == "K")]
    return ltc3_g_stars, ltc3_k_stars
