from provider.utils import (
    fill_sources_table,
    create_sources_table,
    create_provider_table,
    fetch_main_id,
    lower_quality,
    IdentifierCreator,
    OidCreator,
    distance_cut
)
import pytest
import provider.utils as utils_module
from datetime import datetime
import numpy as np  # arrays
from astropy.table import (
    Table,
    setdiff
)

def test_fill_sources_table():
    # Data
    cat = Table({"para1_ref": ["ref1", "ref2"], "para2_ref": ["ref2", "ref3"]})
    ref_columns = ["para1_ref", "para2_ref"]
    provider = "SIMBAD"

    # Execute
    return_sources = fill_sources_table(cat, ref_columns, provider)

    # Verify
    expected_sources = Table(
        {
            "ref": ["ref1", "ref2", "ref3"],
            "provider_name": ["SIMBAD" for j in range(3)],
        }
    )
    assert len(setdiff(return_sources, expected_sources)) == 0


def test_create_sources_table():
    # Data
    cat = Table({"para1_ref": ["ref1", "ref2"], "para2_ref": ["ref2", "ref3"]})
    cat2 = Table({"para3_ref": ["ref4", "ref2"], "para4_ref": ["ref2", "ref5"]})
    tables = [cat, cat2]
    ref_columns = [["para1_ref", "para2_ref"], ["para3_ref", "para4_ref"]]
    provider_name = "SIMBAD"

    # Execute
    return_sources = create_sources_table(tables, ref_columns, provider_name)

    # Verify
    expected_sources = Table(
        {
            "ref": ["ref1", "ref2", "ref3", "ref4", "ref5"],
            "provider_name": ["SIMBAD" for j in range(5)],
        }
    )
    assert len(setdiff(return_sources, expected_sources)) == 0


def test_create_provider_table_date_given():
    gk_provider = create_provider_table(
        "Grant Kennedy Disks",
        "http://drgmk.com/sdb/",
        "priv. comm.",
        "2024-02-09",
    )
    assert gk_provider["provider_name"] == "Grant Kennedy Disks"
    assert gk_provider["provider_url"] == "http://drgmk.com/sdb/"
    assert gk_provider["provider_bibcode"] == "priv. comm."
    assert gk_provider["provider_access"] == "2024-02-09"


def test_create_provider_table_no_date_given():
    gk_provider = create_provider_table(
        "Grant Kennedy Disks", "http://drgmk.com/sdb/", "priv. comm."
    )
    assert gk_provider["provider_access"] == datetime.now().strftime("%Y-%m-%d")


# @pytest.mark.dependency()
def test_fetch_main_id():
    # data
    id_cat = Table(
        data=[["HD 166", "6 Cet"]], names=["sdb_host_main_id"], dtype=[object]
    )
    oid_cat = Table(data=[[1800872, 508940]], names=["parent_oid"], dtype=[int])

    # function
    id_result = fetch_main_id(
        id_cat,
        id_creator=IdentifierCreator(
            name="main_id", colname="sdb_host_main_id"
        ),
    )
    oid_result = fetch_main_id(
        oid_cat,
        id_creator=OidCreator(name="parent_main_id", colname="parent_oid"),
    )

    # assert
    assert (
        id_result["main_id"][
            np.where(id_result["sdb_host_main_id"] == "HD 166")
        ]
        == "HD    166"
    )
    assert (
        id_result["main_id"][np.where(id_result["sdb_host_main_id"] == "6 Cet")]
        == "*   6 Cet"
    )
    assert (
        oid_result["parent_main_id"][
            np.where(oid_result["parent_oid"] == 1800872)
        ]
        == "* ksi UMa"
    )
    assert (
        oid_result["parent_main_id"][
            np.where(oid_result["parent_oid"] == 508940)
        ]
        == "MCC 541"
    )


def test_lower_quality():
    assert lower_quality("A") == "B"
    assert lower_quality("B") == "C"
    assert lower_quality("C") == "D"
    assert lower_quality("D") == "E"
    assert lower_quality("E") == "E"


def _make_sim_objects() -> Table:
    """
    Create a minimal sim_objects table for testing joins by main_id.

    :returns: Table with main_id and ids columns.
    :rtype: astropy.table.Table
    """
    return Table(
        {
            "main_id": np.array(["A", "B"], dtype=object),
            "ids": np.array(["A|A1", "B|B1"], dtype=object),
        }
    )


def _make_sim_ident() -> Table:
    """
    Create a minimal sim_ident table for testing joins by identifier.

    :returns: Table with id and main_id columns.
    :rtype: astropy.table.Table
    """
    return Table(
        {
            "id": np.array(["A","A1", "B1","B", "C1","C"], dtype=object),
            "main_id": np.array(["A", "A","B", "B","C","C"], dtype=object),
        }
    )


def test_distance_cut_by_main_id(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Keep only rows whose main_id exists in sim_objects.

    Verifies:
    - unmatched rows are filtered out,
    - the returned table preserves the original schema (no temp columns).
    """
    # Input catalog with one non-matching row ("X")
    cat = Table(
        {
            "main_id": np.array(["A", "X", "B"], dtype=object),
            "value": np.array([1, 2, 3]),
        }
    )

    def fake_load(names):
        # utils.distance_cut loads ["sim_objects"] in main_id=True mode
        assert names == ["sim_objects"]
        return [_make_sim_objects()]

    monkeypatch.setattr(utils_module, "load", fake_load)

    out = distance_cut(cat, colname="main_id", main_id=True)

    # Only A and B remain
    assert set(out["main_id"]) == set(["A", "B"])
    assert set(out["value"]) == set([1, 3])

    # No temp columns leaked
    assert "temp1" not in out.colnames
    assert "temp2" not in out.colnames


def test_distance_cut_by_identifier(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Keep only rows whose identifier exists in sim_ident.

    Verifies:
    - unmatched rows are filtered out,
    - the returned table contains main_id coming from sim_ident.
    """
    # Input catalog with identifiers; "ID-X" is not present in sim_ident
    cat = Table(
        {
            "ext_id": np.array(["A", "X", "C1"], dtype=object),
            "note": np.array(["ok", "drop", "ok"], dtype=object),
        }
    )

    def fake_load(names):
        # utils.distance_cut loads ["sim_ident"] in main_id=False mode
        assert names == ["sim_ident"]
        return [_make_sim_ident()]

    monkeypatch.setattr(utils_module, "load", fake_load)

    out = distance_cut(cat, colname="ext_id", main_id=False)

    # Only rows with matching IDs remain; main_id is added by the join
    assert list(out["ext_id"]) == ["A", "C1"]
    assert list(out["note"]) == ["ok", "ok"]
    assert list(out["main_id"]) == ["A", "C"]

    # Temp column removed
    assert "temp1" not in out.colnames

