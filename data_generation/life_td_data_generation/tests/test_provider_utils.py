from provider.utils import (
    fill_sources_table,
    create_sources_table,
    create_provider_table,
    fetch_main_id,
    lower_quality,
    IdentifierCreator,
    OidCreator
)

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
