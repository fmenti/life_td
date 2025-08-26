import provider.simbad as simbad_module
import pytest
from astropy.table import MaskedColumn, Table, setdiff
from provider.simbad import (
    create_simbad_helpertable,
    creating_helpertable_stars,
    create_ident_table,
    create_h_link_table,
    stars_in_multiple_system,
    expanding_helpertable_stars,
    create_objects_table,
    create_sim_sources_table
)
import numpy as np  # arrays
from provider.utils import nullvalues

@pytest.fixture
def query_returns(monkeypatch):
    # Shared setup for tests that need the helper table
    distance_cut_in_pc = 5
    test_objects = []

    t0 = Table(
        {
            "main_id": np.array(
                ["star1", "star2", "star3", "star4"], dtype=object
            ),
            "coo_ra": MaskedColumn([160.0, 165.0, 100.0, 60.0]),
            "coo_dec": MaskedColumn([35.0, 45.0, 30.0, 20.0]),
            "coo_err_angle": MaskedColumn(
                [90.0, 17.0, np.ma.masked, np.ma.masked]
            ),
            "coo_err_maj": MaskedColumn(
                [0.02, 180.0, np.ma.masked, np.ma.masked]
            ),
            "coo_err_min": MaskedColumn(
                [0.01, 160.0, np.ma.masked, np.ma.masked]
            ),
            "oid": [1, 2, 3, 6],
            "coo_ref": np.array(
                ["ref_coo", "ref_coo", "ref_coo2", "ref_coo"], dtype=object
            ),
            "coo_qual": ["A", "A", "C", "C"],
            "sptype_string": np.array(
                ["M7.0V", "", "M5.0V+M9", ""], dtype=object
            ),
            "sptype_qual": ["C", "", "B", ""],
            "sptype_ref": np.array(
                ["ref_sptype", "", "ref_sptype2", ""], dtype=object
            ),
            "plx_err": [0.05, 0.9, 3, 2],
            "plx_value": MaskedColumn([190.0, 230.0, 375.0, 300.0]),
            "plx_ref": np.array(
                ["ref_plx", "ref_plx", "ref_plx2", "ref_plx"], dtype=object
            ),
            "plx_qual": MaskedColumn(data=["A", "A", "C", "A"], fill_value="N"),
            "membership": MaskedColumn(
                [np.ma.masked, np.ma.masked, 100, np.ma.masked]
            ),
            "parent_oid": MaskedColumn([4, np.ma.masked, 4, np.ma.masked]),
            "h_link_ref": np.array(
                ["ref_hlink", "", "ref_hlink", ""], dtype=object
            ),
            "otypes": np.array(["*|**", "BD?", "*|**", "*"], dtype=object),
            "ids": np.array(
                ["star1|id1", "star2|id2", "star3|id3", "star4|id4"],
                dtype=object,
            ),
            "mag_i_value": MaskedColumn(
                [np.ma.masked, np.ma.masked, 7.0, np.ma.masked]
            ),
            "mag_j_value": MaskedColumn([8.0, np.ma.masked, 6.0, np.ma.masked]),
            "mag_k_value": MaskedColumn([9.0, np.ma.masked, 8.0, np.ma.masked]),
        }
    )
    parents_without_plx = t0[:0].copy()
    parents_without_plx.add_row(
        [
            "system1",
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
            4,
            "",
            "",
            "",
            "",
            "",
            np.ma.masked,
            np.ma.masked,
            "",
            "",
            np.ma.masked,
            np.ma.masked,
            "",
            "**",
            "system1",
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
        ]
    )
    children_without_plx = t0[:0].copy()
    children_without_plx.add_row(
        [
            "planet1",
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
            5,
            "",
            "",
            "",
            "",
            "",
            np.ma.masked,
            np.ma.masked,
            "",
            "",
            np.ma.masked,
            3,
            "ref_hlink2",
            "Pl",
            "planet1",
            np.ma.masked,
            np.ma.masked,
            np.ma.masked,
        ]
    )

    ident = Table(
        {
            "main_id": np.array(
                [
                    "star1",
                    "star1",
                    "star3",
                    "star3",
                    "system1",
                    "planet1",
                    "star4",
                    "star4",
                ],
                dtype=object,
            ),
            "oid": [1, 1, 3, 3, 4, 5, 6, 6],
            "id": np.array(
                [
                    "star1",
                    "id1",
                    "star3",
                    "id3",
                    "system1",
                    "planet1",
                    "star4",
                    "id4",
                ],
                dtype=object,
            ),
        }
    )

    tables_iter = iter([t0, parents_without_plx, children_without_plx, ident])

    def fake_query(url, adql, uploads=None):
        return next(tables_iter)

    monkeypatch.setattr(simbad_module, "query", fake_query)

    result_sim_helptab, result_sim = create_simbad_helpertable(
        distance_cut_in_pc, test_objects
    )

    result_sim["ident"] = create_ident_table(result_sim_helptab, result_sim)

    return result_sim_helptab, result_sim, t0, result_sim["ident"]


def test_create_simbad_helpertable(query_returns):
    result_sim_helptab, result_sim, t0, ident = query_returns

    # Verify
    # all helptab columns are present
    assert set(result_sim_helptab.colnames) == set(
        t0.colnames + ["type", "binary_flag"]
    )
    assert set(result_sim_helptab["type"]) == set(
        ["sy", "sy", "sy", "pl", "st"]
    )
    # objects got correctly added and removed from helptab
    assert (
        len(result_sim_helptab) == 5
    )  # stars, systems and planets but removed star2
    # provider table is correctly filled in
    assert len(result_sim["provider"]) == 1


def test_creating_helpertable_stars(query_returns):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns

    # Execute
    result_stars = creating_helpertable_stars(result_sim_helptab)

    # Verify
    expected_stars = result_sim_helptab
    expected_stars.remove_row(-1)
    for null, colname in zip(
        [
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
        ],
        [
            "membership",
            "coo_err_angle",
            "coo_err_maj",
            "coo_err_min",
            "parent_oid",
            "mag_i_value",
            "mag_j_value",
            "mag_k_value",
            "coo_ra",
            "coo_dec",
            "plx_value",
        ],
    ):
        result_stars = nullvalues(result_stars, colname, null)
        expected_stars = nullvalues(expected_stars, colname, null)

    assert len(result_stars) == 4  # stars and systems but removed planets
    assert len(setdiff(result_stars, expected_stars)) == 0
    assert len(setdiff(expected_stars, result_stars)) == 0


def test_create_ident_table(query_returns):
    # Execute
    result_sim_helptab, result_sim, t0, return_ident = query_returns

    # Verify
    ref = result_sim["provider"]["provider_bibcode"][0]
    expected_ident = Table(
        {
            "main_id": np.array(
                [
                    "star1",
                    "star1",
                    "star3",
                    "star3",
                    "system1",
                    "planet1",
                    "star4",
                    "star4",
                ],
                dtype=object,
            ),
            "id_ref": np.array(
                [ref, ref, ref, ref, ref, ref, ref, ref], dtype=object
            ),
            "id": np.array(
                [
                    "id1",
                    "star1",
                    "id3",
                    "star3",
                    "system1",
                    "planet1",
                    "star4",
                    "id4",
                ],
                dtype=object,
            ),
        }
    )

    assert len(setdiff(return_ident, expected_ident)) == 0
    assert len(setdiff(expected_ident, return_ident)) == 0


def test_create_h_link_table(query_returns, monkeypatch):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    result_stars = creating_helpertable_stars(result_sim_helptab)
    main_id_table = Table(
        {
            "main_id": np.array(["star1", "star3", "planet1"], dtype=object),
            "parent_main_id": np.array(
                ["system1", "system1", "star3"], dtype=object
            ),
            "h_link_ref": np.array(
                ["ref_hlink", "ref_hlink", "ref_hlink2"], dtype=object
            ),
            "membership": [np.ma.masked, 100, np.ma.masked],
            "parent_oid": [1, 4, 3],
        }
    )

    def fake_fetch_main_id(sim_h_link, id_creator):
        return main_id_table

    monkeypatch.setattr(simbad_module, "fetch_main_id", fake_fetch_main_id)

    # Execute
    result_h_link = create_h_link_table(
        result_sim_helptab, result_sim, result_stars
    )
    # issue: wrong null value in membership column. find out what simbad returns, str, need to have them as object.
    # sim returns int and masked entries

    expected_h_link = Table(
        {
            "main_id": np.array(["star1", "star3", "planet1"], dtype=object),
            "parent_main_id": np.array(
                ["system1", "system1", "star3"], dtype=object
            ),
            "h_link_ref": np.array(
                ["ref_hlink", "ref_hlink", "ref_hlink2"], dtype=object
            ),
            "membership": [999999, 100, 999999],
        }
    )

    # Verify
    assert len(setdiff(result_h_link, expected_h_link)) == 0
    assert len(setdiff(expected_h_link, result_h_link)) == 0


def test_stars_in_multiple_system():
    # Data
    cat = Table(
        {
            "main_id": np.array(["system1", "star1", "star2"], dtype=object),
            "type": np.array(["sy", "sy", "sy"], dtype=object),
            "sptype_string": np.array(
                ["M7.0V+M5", "M7.0V", "M5.0V+M9"], dtype=object
            ),
        }
    )
    sim_h_link = Table(
        {
            "main_id": np.array(["star1", "star2", "planet1"], dtype=object),
            "parent_main_id": np.array(
                ["system1", "system1", "star1"], dtype=object
            ),
            "h_link_ref": np.array(["ref1", "ref1", "ref2"], dtype=object),
        }
    )
    all_objects = Table(
        {
            "type": np.array(["sy", "sy", "sy", "pl"], dtype=object),
            "main_id": np.array(
                ["system1", "star1", "star2", "planet1"], dtype=object
            ),
        }
    )
    # Execute
    return_cat = stars_in_multiple_system(cat, sim_h_link, all_objects)

    # Verify
    expected_cat = Table(
        {
            "main_id": np.array(["system1", "star1", "star2"], dtype=object),
            "type": np.array(["sy", "st", "sy"], dtype=object),
            "sptype_string": np.array(
                ["M7.0V+M5", "M7.0V", "M5.0V+M9"], dtype=object
            ),
        }
    )
    assert len(setdiff(return_cat, expected_cat)) == 0
    assert len(setdiff(expected_cat, return_cat)) == 0


def expanded_helptab_stars(expected_stars, ref):
    expected_stars.remove_row(-1)  # removing planet object row
    # doing stars in multiple systems function
    expected_stars["binary_flag"] = np.array(
        ["True" for j in range(len(expected_stars))], dtype=object
    )
    expected_stars["binary_flag"][
        np.where(expected_stars["main_id"] == "star4")
    ] = "False"
    # expected_stars['plx_qual'][np.where(expected_stars['main_id']=='system1')] = '?'
    expected_stars["plx_qual"] = expected_stars["plx_qual"].astype(object)

    expected_stars["mag_i_ref"] = MaskedColumn(
        data=np.array(["", ref, "", ""], dtype=object),
        mask=[True, False, True, True],
        fill_value="N",
    )
    expected_stars["mag_j_ref"] = MaskedColumn(
        data=np.array([ref, ref, "", ""], dtype=object),
        mask=[False, False, True, True],
        fill_value="N",
    )
    expected_stars["mag_k_ref"] = MaskedColumn(
        data=np.array([ref, ref, "", ""], dtype=object),
        mask=[False, False, True, True],
        fill_value="N",
    )
    for col in ["sptype_ref", "coo_ref", "plx_ref"]:
        expected_stars[col] = MaskedColumn(expected_stars[col])
        expected_stars[col].mask[
            np.where(expected_stars["main_id"] == "system1")
        ] = True
        expected_stars[col].mask[np.where(expected_stars[col] == "")] = True

    expected_stars["binary_ref"] = [ref, ref, ref, ref]
    expected_stars["binary_qual"] = ["D", "D", "D", "D"]
    expected_stars["type"][np.where(expected_stars["main_id"] == "star1")] = (
        "st"  # stars in multiple systems function
    )
    return expected_stars


def test_expanding_helpertable_stars(query_returns):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    stars = result_sim_helptab
    expected_stars = expanded_helptab_stars(
        stars.copy(), result_sim["provider"]["provider_bibcode"][0]
    )

    stars.remove_row(-1)  # removing planet object row
    result_sim["h_link"] = Table(
        {
            "main_id": np.array(["star1", "star3", "planet1"], dtype=object),
            "parent_main_id": np.array(
                ["system1", "system1", "star3"], dtype=object
            ),
            "h_link_ref": np.array(
                ["ref_hlink", "ref_hlink", "ref_hlink2"], dtype=object
            ),
            "membership": [999999, 100, 999999],
        }
    )

    # Execute
    stars["plx_qual"].mask = [False, False, False, True]
    result_stars = expanding_helpertable_stars(
        result_sim_helptab, result_sim, stars.copy()
    )

    # Verify
    for null, colname in zip(
        [
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            999999,
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            999999,
            999999,
            999999,
            "",
            "",
            "",
        ],
        [
            "membership",
            "coo_err_angle",
            "coo_err_maj",
            "coo_err_min",
            "parent_oid",
            "mag_i_value",
            "mag_j_value",
            "mag_k_value",
            "mag_i_ref",
            "mag_j_ref",
            "mag_k_ref",
            "plx_qual",
            "coo_qual",
            "sptype_string",
            "sptype_qual",
            "h_link_ref",
            "coo_ra",
            "coo_dec",
            "plx_value",
            "coo_ref",
            "sptype_ref",
            "plx_ref",
        ],
    ):
        result_stars = nullvalues(result_stars, colname, null)
        expected_stars = nullvalues(expected_stars, colname, null)

    assert len(setdiff(result_stars, expected_stars)) == 0
    assert len(setdiff(expected_stars, result_stars)) == 0


def test_create_object_table(query_returns):
    # Execute
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    stars = result_sim_helptab.copy()
    stars.remove_row(-1)  # removing planet object row

    stars["type"][np.where(stars["main_id"] == "star1")] = (
        "st"  # stars in multiple systems function
    )

    return_object = create_objects_table(result_sim_helptab, stars)

    # Verify
    expected_object = Table(
        {
            "main_id": np.array(
                ["star1", "star3", "system1", "planet1", "star4"], dtype=object
            ),
            "type": np.array(["st", "sy", "sy", "pl", "st"], dtype=object),
            "ids": np.array(
                ["star1|id1", "star3|id3", "system1", "planet1", "star4|id4"],
                dtype=object,
            ),
        }
    )

    assert len(setdiff(return_object, expected_object)) == 0
    assert len(setdiff(expected_object, return_object)) == 0


def test_create_sim_sources_table(query_returns):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    stars = result_sim_helptab.copy()
    expected_stars = expanded_helptab_stars(
        stars.copy(), result_sim["provider"]["provider_bibcode"][0]
    )
    result_sim["h_link"] = Table(
        {
            "main_id": np.array(["star1", "star3", "planet1"], dtype=object),
            "parent_main_id": np.array(
                ["system1", "system1", "star3"], dtype=object
            ),
            "h_link_ref": np.array(
                ["ref_hlink", "ref_hlink", "ref_hlink2"], dtype=object
            ),
            "membership": [999999, 100, 999999],
        }
    )

    # Execute
    return_sources = create_sim_sources_table(expected_stars, result_sim)

    # Verify
    expected_sources = Table(
        {
            "ref": [
                "ref_coo",
                "ref_coo2",
                "ref_sptype",
                "ref_sptype2",
                "ref_plx",
                "ref_plx2",
                "ref_hlink",
                "ref_hlink2",
                result_sim["provider"]["provider_bibcode"][0],
            ],
            "provider_name": ["SIMBAD" for j in range(9)],
        }
    )

    assert len(setdiff(return_sources, expected_sources)) == 0
    assert len(setdiff(expected_sources, return_sources)) == 0

    # to do next:
    # clean up provider_sim while making sure that tests still work
    # clean up tests for provider sim afterwards
