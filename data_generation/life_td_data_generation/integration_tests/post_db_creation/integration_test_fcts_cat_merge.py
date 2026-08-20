from catalog.starcat5_merger.fcts_cat_merge import (
    _match_cat2_to_cat1_by_name,
    _match_cat2_to_cat1_by_coords,
    get_cat2_in_cat1_match_info,
)
from astropy.io.ascii import read
from utils.io import load
from astropy.table import setdiff
import numpy as np

def test_nearest_neighbor_distances_units_return_works_for_coord_match():
    # data
    hpic = read("../../additional_data/HPICv1.0/full_HPIC.txt",
                delimiter="|")
    #hpic = get_catalog("hpic")
    [starcat5] = load(["catalogs/StarCat5"],
                      location="../../additional_data/")

    radius =  103.43301245378981

    # only get name matches
    cat1_index_for_cat2, match_by_name = _match_cat2_to_cat1_by_name(
        hpic["simbad_name"],
        starcat5["main_id"],
    )

    name_matched = starcat5[match_by_name]


    # do coordinate matching only on the ones that were previously name matched
    coord_indices, coord_mask, sep_arcsec = _match_cat2_to_cat1_by_coords(
        hpic["ra"],
        hpic["dec"],
        name_matched["coo_ra"],
        name_matched["coo_dec"],
        r_arcsec = radius,
    )
    coord_matched = name_matched[coord_mask]

    matched_objects_by_name = cat1_index_for_cat2[match_by_name]
    matched_objects_by_coord = coord_indices

    name_matched["ind_matched_objects_by_name"] = cat1_index_for_cat2[match_by_name]
    coord_matched["ind_matched_objects_by_coord"] = coord_indices

    name_matched["ind"] = name_matched["ind_matched_objects_by_name"]
    coord_matched["ind"] = coord_matched["ind_matched_objects_by_coord"]

    sdiff_name = setdiff(name_matched["main_id","ind"],
                    coord_matched["main_id","ind"],
                    keys=['ind'])
    sdiff_coord = setdiff(
                         coord_matched["main_id", "ind"],
                         name_matched["main_id", "ind"],
                         keys=['ind'])
    print(sdiff_name)
    # * alf Cen A    7
    #   HD 352860  326
    # *  70 Oph A  599
    # BD+66    34B 7093
    print("\n")
    print(sdiff_coord)
    # BD+66    34B 7092 -> wrongly matched
    #    HD 352860 8535

    # since in sdiff_coord two stars not shown see if those are correctly matched -> no
    print("\n")
    print(coord_matched["main_id","ind"][np.where(coord_matched["main_id"] == "* alf Cen A")])
    # 5
    print(hpic["simbad_name"][5])
    # * alf Cen B -> wrongly matched
    print(coord_matched["main_id", "ind"][
              np.where(coord_matched["main_id"] == "*  70 Oph A")])
    # 928
    print(hpic["simbad_name"][928])
    # *  70 Oph B -> wrongly matched

    # check all the matches if correct
    print(hpic["simbad_name"][7])
    # correct match
    print(hpic["simbad_name"][326])
    # correct match
    print(hpic["simbad_name"][599])
    # correct match
    print(hpic["simbad_name"][7093])
    # correct match
    print(hpic["simbad_name"][7092])
    # wrong match BD+66    34A instead of B
    print(hpic["simbad_name"][8535])
    # correct match

    # how can both hd 352860 be correctly matched but different indices?
    # check again -> same result
    print(hpic["simbad_name","star_name","ra","dec"][np.where(hpic["simbad_name"] == "HD 352860")])
    # different TIC objects with very similar coordinates

    # fazit: if leaving name away I get 4 wrongly matched objects, all binaries
    # (* alf Cen A, HD 352860, *  70 Oph A, BD+66     34B)

    # now see how that changes with different radii
    # 50 instead of 103.43301245378981 -> all the same 4 mismatched objects
    # 200 instead of 103.43301245378981 -> still same 4 mismatched objects
    # 10 -> error message
    # 1000 -> still same 4 mismatched objects
    # does this outcome make sense? I would have expected it to be dependent on r
    # but maybe pm missing is reason that stuff independent of r?
    # is something wrong with my coord fitting function?

    # so what does radius actually mean?
    # nearest neighbor separations of a whole catalog computed
    # pick a treshold where 90% exponential decay.
    # so most objects have nearest neighbors distances below this radius
    # since the catalog I inputed has multiple objects that are the same (sy, st)
    # ...maybe I can use report on exoplanet false positive test function here
    # look in exact detail what happens to alf cen pair
    # find its nearest neighbor distances -> should be in match info
    # using it in next test below so try continuing work there

    # run into issue of no coords / wrong ident there
    print(coord_matched["main_id", "ind", "coo_ra", "coo_dec"][
              np.where(coord_matched["main_id"] == "BD+66    34B")])
    # but here they have the coords... so why not below?
    # now just copy pasted again this id and worked, but only for B not A

    # match by coord matches them, if skycoord distance < radius


    # find out if I miss objects or get false positives
    assert len(name_matched) == len(coord_matched)
    # 13247 != 15850 when using starcat5 instead of name_matched as base
    # working fine else -> confirmes that all names get recovered
    # well at least number of objects is same, could be matched up with wrong objects still
    # does it also rule out, that false positives are created?
    # maybe do coord match for all and confirm that is same length as with coord and names match

    # assert list(cat1_index_for_cat2[match_by_name]) == list(coord_indices)
    # At index 7747 diff: np.int64(326) != np.int64(8535)
    # full diff...
    # so how many different objects do I have? maybe do via setdiff from ap
    # 4 objects are different
    # which are they wrongly paired to?
    # also in end try out with different radiusd
    assert False

def test_false_positives_coord_match():
    # data
    hpic = read("../../additional_data/HPICv1.0/full_HPIC.txt",
                delimiter="|")
    # hpic = get_catalog("hpic")
    [starcat5] = load(["catalogs/StarCat5"],
                      location="../../additional_data/")

    radius = 103.43301245378981

    # do only coordinate matching
    coord_indices, coord_mask, sep_arcsec = _match_cat2_to_cat1_by_coords(
        hpic["ra"],
        hpic["dec"],
        starcat5["coo_ra"],
        starcat5["coo_dec"],
        r_arcsec=radius,
    )
    coord_matched = starcat5[coord_mask]

    # do normal matching
    match_info = get_cat2_in_cat1_match_info(
        name_cat1=hpic["simbad_name"],
        ra_cat1=hpic["ra"],
        dec_cat1=hpic["dec"],
        name_cat2=starcat5["main_id"],
        ra_cat2=starcat5["coo_ra"],
        dec_cat2=starcat5["coo_dec"],
        r_arcsec=radius,
    )

    print(match_info["sep_arcsec"][np.where(starcat5["main_id"] == "BD+66    34A")])
    # 11.88040573
    print(match_info["sep_arcsec"][np.where(starcat5["main_id"] == "BD+66    34B")])
    # 13.89590673

    # ok so radius would need to be very small (11, 13) to not match the alf cens
    # test for others that went wrong above
    # 11 & 13 for alf cen A and B
    # 0.01746714 HD 352860
    # 5 and 0.01 for 70 oph A and B
    # ? and 0.01 no nearest neighbor distances for BD+66     34A and B -> ??
    print(starcat5["coo_ra","coo_dec"][np.where(starcat5["main_id"] == "BD+66    34A")])
    print(starcat5["coo_ra", "coo_dec"][
              np.where(starcat5["main_id"] == "BD+66    34B")])

    # why do I not have coordinates for BD+66      34A and B?
    # simbad does have coords for them
    # are also in db with coord but with slightly different name syntax: BD+66    34A

    # even with that name sntax not found

    all_matched_len = (len(starcat5[match_info["match_by_name"]]) +
                       len(starcat5[match_info["match_by_coords"]]))

    assert len(coord_matched) == all_matched_len
    assert False
