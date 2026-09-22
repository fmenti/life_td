from astropy.table import Table
from provider.wds import create_star_basic, assign_names
from utils.io import load
import numpy as np
from sdata import empty_dict

def test_create_star_basic():
    # data
    wds: dict[str, Table] = empty_dict.copy()
    [wds_helptab,wds_ident] = load(["wds_helptab","wds_ident"])
    wds["ident"] = wds_ident

    wds_star_basic = create_star_basic(wds_helptab,wds)

    # test if there are empty main_id entries
    assert len(wds_star_basic["main_id"].mask.nonzero()[0]) == 0
    assert len(wds_star_basic["main_id"][np.where(wds_star_basic["main_id"]=="")]) == 0


    # test if all main_id entries are also in ident table
    mask = np.isin(wds_star_basic["main_id"], wds_ident["main_id"])
    main_id_not_in_ident = wds_star_basic["main_id"][~mask]
    assert len(main_id_not_in_ident) == 0


def test_assign_names():
    # data
    [wds_querried] = load(["wds_querried"])

    wds_helptab = assign_names(wds_querried)

    # how comes that I have entries of primary and secondary that are empty?
    # can wds_name which is basic for those columns be empty? -> no...
    print(wds_helptab["wds_name","primary","secondary"])

    assert len(wds_helptab["wds_name"].mask.nonzero()[0]) == 0
    assert len(
        wds_helptab["wds_name"][np.where(wds_helptab["wds_name"]=="")]) == 0
    assert len(wds_helptab["primary"].mask.nonzero()[0]) == 0
    assert len(
        wds_helptab["primary"][np.where(wds_helptab["primary"] == "")]) == 0
    assert len(wds_helptab["secondary"].mask.nonzero()[0]) == 0
    assert len(
        wds_helptab["secondary"][np.where(wds_helptab["secondary"] == "")]) == 0

    # so why where there so many empty / masked entries in the star_basic creation one?
    # this is before the distance cut on system, primary and secondary
    # because I can't use main id from helptab but need to use ident table
