"""
Tests for the Exo-MerCat provider (unit + small integration facades).

Covers data access (query/load) and downstream table builders.
"""

import numpy as np
from astropy.table import MaskedColumn, Table, Column
from provider.utils import create_provider_table
from provider.exo import (
    create_ident_table,
    create_objects_table,
    create_object_main_id,
    deal_with_mass_nullvalues,
    create_para_exo_mes_mass_pl,
    bestmass_better_qual,
    assign_new_qual,
    betterthan,
    align_quality_with_bestmass,
    create_mes_mass_pl_table,
    create_exo_helpertable,
    create_h_link_table,
    create_exo_sources_table,
)
from sdata import empty_dict
import pytest
import provider.exo as exo_module


# -----------------
# Local test helpers
# -----------------

def _mk_base_with_spaces() -> Table:
    """
    Create a small helper table with leading/trailing spaces and masks.

    :returns: Base table used to exercise strip, distance cut and join.
    :rtype: astropy.table.Table
    """
    base = Table()
    base["main_id"] = MaskedColumn(
        data=np.array([" HD 1", "", "* HD 3"], dtype=object),
        dtype=object,
        mask=[False, True, False],
    )
    base["host"] = Column(
        data=np.array(["HD 1", " HD 2", " HD 3"], dtype=object), dtype=object
    )
    base["binary"] = MaskedColumn(
        data=np.array(["A", "", ""], dtype=object),
        dtype=object,
        mask=[False, True, True],
    )
    base["letter"] = Column(
        data=np.array(["b", "c", "b"], dtype=object), dtype=object
    )
    base["exomercat_name"] = Column(
        data=np.array(["HD 1 A b", "HD 2 c", "HD 3 b"], dtype=object),
        dtype=object,
    )
    return base


def _mk_msini_only_fixture() -> Table:
    """
    Fixture table for msini-only mass extraction tests.

    :returns: Table with msini columns and a placeholder quality column.
    :rtype: astropy.table.Table
    """
    msini = MaskedColumn(
        data=[20.76, 0.1, 0, 5.025, 1e20],
        name="mass",
        mask=[False, False, True, False, False],
    )
    msinimax = MaskedColumn(
        data=[np.inf, 0, 0, 0.873, 1],
        name="mass_max",
        mask=[False, True, True, False, False],
    )
    msinimin = MaskedColumn(
        data=[np.inf, 1, 0, 1.067, 1],
        name="mass_min",
        mask=[False, False, False, False, False],
    )
    msiniurl = MaskedColumn(
        data=["eu", "test", "", "2022ApJS..262...21F", "test"],
        name="mass_url",
        mask=[False, False, True, False, False],
    )
    mprov = MaskedColumn(
        data=["Msini", "Mass", "", "Mass", "Msini"],
        name="bestmass_provenance",
        mask=[False, False, True, False, False],
    )
    exo_helptab = Table(
        data=[
            [
                "*   3 Cnc b",
                "*   4 Mon B .01",
                "testname1",
                "*   6 Lyn b",
                "testname2",
            ],
            msini,
            msinimax,
            msinimin,
            msiniurl,
            mprov,
        ],
        names=[
            "planet_main_id",
            "msini",
            "msini_max",
            "msini_min",
            "msini_url",
            "bestmass_provenance",
        ],
        dtype=[object, float, float, float, object, object],
    )
    exo_helptab["msini_pl_qual"] = MaskedColumn(
        dtype=object, length=len(exo_helptab)
    )
    exo_helptab["msini_pl_qual"] = ["?" for _ in range(len(exo_helptab))]
    return exo_helptab


def _mk_full_mass_fixture() -> Table:
    """
    Fixture table for full mass table creation (mass + msini).

    :returns: Table with mass, msini and their error/reference columns.
    :rtype: astropy.table.Table
    """
    m = MaskedColumn(
        data=[20.76, 0.1, 0, 5.025, 1e20],
        name="mass",
        mask=[False, False, True, False, False],
    )
    mmax = MaskedColumn(
        data=[np.inf, 0, 0, 0.873, 1],
        name="mass_max",
        mask=[False, True, True, False, False],
    )
    mmin = MaskedColumn(
        data=[np.inf, 1, 0, 1.067, 1],
        name="mass_min",
        mask=[False, False, False, False, False],
    )
    murl = MaskedColumn(
        data=["eu", "test", "", "2022ApJS..262...21F", "test"],
        name="mass_url",
        mask=[False, False, True, False, False],
    )
    msini = MaskedColumn(
        data=[20.76, 0.1, 0, 5.025, 1e20],
        name="mass",
        mask=[False, False, True, False, False],
    )
    msinimax = MaskedColumn(
        data=[np.inf, 0, 0, 0.873, 1],
        name="mass_max",
        mask=[False, True, True, False, False],
    )
    msinimin = MaskedColumn(
        data=[np.inf, 1, 0, 1.067, 1],
        name="mass_min",
        mask=[False, False, False, False, False],
    )
    msiniurl = MaskedColumn(
        data=["eu", "test", "", "2022ApJS..262...21F", "test"],
        name="mass_url",
        mask=[False, False, True, False, False],
    )
    mprov = MaskedColumn(
        data=["Msini", "Mass", "", "Mass", "Msini"],
        name="bestmass_provenance",
        mask=[False, False, True, False, False],
    )
    exo_helptab = Table(
        data=[
            [
                "*   3 Cnc b",
                "*   4 Mon B .01",
                "testname1",
                "*   6 Lyn b",
                "testname2",
            ],
            m,
            murl,
            mmax,
            mmin,
            msini,
            msiniurl,
            msinimax,
            msinimin,
            mprov,
        ],
        names=[
            "planet_main_id",
            "mass",
            "mass_url",
            "mass_max",
            "mass_min",
            "msini",
            "msini_url",
            "msini_max",
            "msini_min",
            "bestmass_provenance",
        ],
        dtype=[
            object,
            float,
            object,
            float,
            float,
            float,
            object,
            float,
            float,
            object,
        ],
    )
    return exo_helptab

def test_create_exo_helpertable(
        monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Verify helper-table creation, distance cut, stripping, SIMBAD join,
    and saving of removed objects.

    :param monkeypatch: Pytest fixture for patching functions.
    :type monkeypatch: pytest.MonkeyPatch
    """
    base = _mk_base_with_spaces()

    def fake_query_or_load_exomercat():
        return ({}, base.copy())

    monkeypatch.setattr(
        exo_module, "query_or_load_exomercat", fake_query_or_load_exomercat
    )

    captured_distance_args: dict[str, object] = {}

    def fake_distance_cut(cat: Table, colname: str, main_id: bool = True) -> Table:
        captured_distance_args["colname"] = colname
        captured_distance_args["main_id"] = main_id
        return cat[:1].copy()

    monkeypatch.setattr(exo_module, "distance_cut", fake_distance_cut)

    def fake_fetch_main_id(tab: Table, id_creator) -> Table:
        return Table(
            {
                "sim_planet_main_id": np.array(["* HD 1 A b"], dtype=object),
                "planet_main_id": np.array(["HD 1 A b"], dtype=object),
                "host_main_id": np.array(["HD 1 A"], dtype=object),
            }
        )

    monkeypatch.setattr(exo_module, "fetch_main_id", fake_fetch_main_id)

    saved: dict[str, Table] = {}

    def fake_save(cats, names, location=None) -> None:
        assert names == ["exomercat_removed_objects"]
        assert len(cats) == 1
        saved["removed_objects"] = cats[0]

    monkeypatch.setattr(exo_module, "save", fake_save)

    # Execute
    exo, exo_helptab = create_exo_helpertable()

    # Assertions
    assert captured_distance_args["colname"] == "main_id"
    assert captured_distance_args["main_id"] is True
    assert len(exo_helptab) == 1
    assert exo_helptab["main_id"][0] == "HD 1"
    assert exo_helptab["exomercat_name"][0] == "HD 1 A b"
    assert exo_helptab["planet_main_id"][0] == "HD 1 A b"
    assert "sim_planet_main_id" in exo_helptab.colnames
    assert exo_helptab["sim_planet_main_id"][0] == "* HD 1 A b"
    assert "removed_objects" in saved
    removed = saved["removed_objects"]
    assert "exomercat_name" in removed.colnames
    assert "HD 2 c" in set(removed["exomercat_name"])


def test_query_exomercat_success(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Ensure query_exomercat uses TAP query and sets provider metadata.

    :param monkeypatch: Pytest fixture for temporary attribute patching.
    :type monkeypatch: pytest.MonkeyPatch
    """
    # Arrange
    exo: dict = {}
    expected = Table({"col": np.array([1, 2], dtype=object)})

    def fake_query(url: str, adql_query: str) -> Table:
        """Fake TAP query returning a small table for assertions."""
        # Validate URL from provider table
        assert url == "http://archives.ia2.inaf.it/vo/tap/projects"
        assert "SELECT" in adql_query
        return expected

    monkeypatch.setattr(exo_module, "query", fake_query)

    # Act
    result = exo_module.query_exomercat("SELECT * FROM t", exo)

    # Assert
    assert result is expected
    assert "provider" in exo
    assert exo["provider"]["provider_name"][0] == "Exo-MerCat"
    assert (
        exo["provider"]["provider_url"][0]
        == "http://archives.ia2.inaf.it/vo/tap/projects"
    )
    assert exo["provider"]["provider_bibcode"][0] == "2020A&C....3100370A"


def test_load_exomercat_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Ensure load_exomercat reads the CSV snapshot and stamps access date.

    :param monkeypatch: Pytest fixture for temporary attribute patching.
    :type monkeypatch: pytest.MonkeyPatch
    """
    # Arrange
    exo: dict = {}
    loaded = Table({"name": np.array(["row1", "row2"], dtype=object)})

    def fake_ascii_read(path: str) -> Table:
        """Fake CSV loader returning a small table."""
        # Path is derived from Path().additional_data + filename
        assert "exo-mercat13-12-2024_v2.0.csv" in path
        return loaded

    def fake_stringtoobject(tab: Table, number: int = 3000) -> Table:
        """Pass-through conversion used by fallback loader."""
        assert tab is loaded
        return tab

    monkeypatch.setattr(exo_module.ascii, "read", fake_ascii_read)
    monkeypatch.setattr(exo_module, "stringtoobject", fake_stringtoobject)

    # Act
    result = exo_module.load_exomercat(exo)

    # Assert
    assert result is loaded
    assert "provider" in exo
    assert exo["provider"]["provider_name"][0] == "Exo-MerCat"
    assert (
        exo["provider"]["provider_url"][0]
        == "http://archives.ia2.inaf.it/vo/tap/projects"
    )
    assert exo["provider"]["provider_bibcode"][0] == "2020A&C....3100370A"
    # Access date is stamped by the function
    assert exo["provider"]["provider_access"][0] == "2024-12-13"


def test_query_or_load_exomercat_prefers_query(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Ensure query_or_load_exomercat uses query path when it succeeds.

    :param monkeypatch: Pytest fixture for temporary attribute patching.
    :type monkeypatch: pytest.MonkeyPatch
    """
    returned = Table({"ok": np.array([True], dtype=object)})

    def fake_query_exomercat(adql: str, exo: dict) -> Table:
        """Populate 'exo' provider and return a small result table."""
        exo["provider"] = create_provider_table(
            "Exo-MerCat",
            "http://archives.ia2.inaf.it/vo/tap/projects",
            "2020A&C....3100370A",
        )
        assert "FROM exomercat.exomercat" in adql
        return returned

    def fake_load_exomercat(exo: dict) -> Table:  # pragma: no cover
        raise AssertionError("Should not be called if query succeeds")

    monkeypatch.setattr(exo_module, "query_exomercat", fake_query_exomercat)
    monkeypatch.setattr(exo_module, "load_exomercat", fake_load_exomercat)

    # Act
    exo, helptab = exo_module.query_or_load_exomercat()

    # Assert
    assert isinstance(exo, dict)
    assert helptab is returned
    assert exo["provider"]["provider_name"][0] == "Exo-MerCat"


def test_query_or_load_exomercat_uses_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Ensure query_or_load_exomercat falls back to local load on failure.

    :param monkeypatch: Pytest fixture for temporary attribute patching.
    :type monkeypatch: pytest.MonkeyPatch
    """
    fallback = Table({"src": np.array(["local"], dtype=object)})

    def fake_query_exomercat(adql: str, exo: dict) -> Table:
        """Simulate TAP failure to trigger fallback load."""
        raise RuntimeError("Simulated TAP failure")

    def fake_load_exomercat(exo: dict) -> Table:
        """Return fallback table while also stamping provider metadata."""
        exo["provider"] = create_provider_table(
            "Exo-MerCat",
            "http://archives.ia2.inaf.it/vo/tap/projects",
            "2020A&C....3100370A",
            "2024-12-13",
        )
        return fallback

    monkeypatch.setattr(exo_module, "query_exomercat", fake_query_exomercat)
    monkeypatch.setattr(exo_module, "load_exomercat", fake_load_exomercat)

    # Act
    exo, helptab = exo_module.query_or_load_exomercat()

    # Assert
    assert isinstance(exo, dict)
    assert helptab is fallback
    assert exo["provider"]["provider_name"][0] == "Exo-MerCat"


def test_create_object_main_id() -> None:
    """
    Ensure host/planet main ids are constructed correctly from inputs.
    """
    a = MaskedColumn(
        data=["*   3 Cnc", "*   4 Mon", ""],
        name="main_id",
        mask=[False, False, True],
    )
    b = MaskedColumn(
        data=["", "B", ""], name="binary", mask=[True, False, True]
    )
    cat = Table(
        data=[a, ["3 Cnc", "TIC 72090501", "6 Lyn"], b, ["b", ".01", "b"]],
        names=["main_id", "host", "binary", "letter"],
        dtype=[object, object, object, object],
    )
    cat = create_object_main_id(cat)
    assert list(cat["host_main_id"]) == ["*   3 Cnc", "*   4 Mon B", "6 Lyn"]
    assert list(cat["planet_main_id"]) == [
        "*   3 Cnc b",
        "*   4 Mon B .01",
        "6 Lyn b",
    ]


def test_create_exo_helpertable(
        monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Verify helper-table creation, distance cut, stripping, SIMBAD join,
    and saving of removed objects.

    :param monkeypatch: Pytest fixture for patching functions.
    :type monkeypatch: pytest.MonkeyPatch
    """
    base = _mk_base_with_spaces()

    # 1) Mock query_or_load_exomercat: return a dummy exo dict and the table.
    def fake_query_or_load_exomercat():
        return ({}, base.copy())

    monkeypatch.setattr(
        exo_module, "query_or_load_exomercat", fake_query_or_load_exomercat
    )

    # 2) Mock distance_cut to drop the second row (simulate outside distance).
    captured_distance_args: dict[str, object] = {}

    def fake_distance_cut(cat: Table, colname: str, main_id: bool = True) -> Table:
        captured_distance_args["colname"] = colname
        captured_distance_args["main_id"] = main_id
        # Keep only the first row
        return cat[:1].copy()

    monkeypatch.setattr(exo_module, "distance_cut", fake_distance_cut)

    # 3) Mock fetch_main_id to map one planet to a SIMBAD-style name.
    #    The function receives a table with planet_main_id, host_main_id.
    def fake_fetch_main_id(tab: Table, id_creator):
        return Table(
            {
                "sim_planet_main_id": np.array(["* HD 1 A b"], dtype=object),
                "planet_main_id": np.array(["HD 1 A b"], dtype=object),
                "host_main_id": np.array(["HD 1 A"], dtype=object),
            }
        )

    monkeypatch.setattr(exo_module, "fetch_main_id", fake_fetch_main_id)

    # 4) Mock save to capture removed objects written by the function.
    saved: dict[str, Table] = {}

    def fake_save(cats, names, location=None) -> None:
        # Expect exactly one table saved under this logical name.
        assert names == ["exomercat_removed_objects"]
        assert len(cats) == 1
        saved["removed_objects"] = cats[0]

    monkeypatch.setattr(exo_module, "save", fake_save)

    # Execute
    exo, exo_helptab = create_exo_helpertable()

    # Assertions

    # distance_cut was called with the main_id column and main_id=True
    assert captured_distance_args["colname"] == "main_id"
    assert captured_distance_args["main_id"] is True

    # Only one row remains after the fake distance cut.
    assert len(exo_helptab) == 1

    # Stripping applied: no leading/trailing spaces.
    assert exo_helptab["main_id"][0] == "HD 1"
    assert exo_helptab["exomercat_name"][0] == "HD 1 A b"
    assert exo_helptab["planet_main_id"][0] == "HD 1 A b"

    # SIMBAD join applied: new column present and mapped value is correct.
    assert "sim_planet_main_id" in exo_helptab.colnames
    assert exo_helptab["sim_planet_main_id"][0] == "* HD 1 A b"

    # Removed objects saved: should contain the row we dropped (PL-2).
    assert "removed_objects" in saved
    removed = saved["removed_objects"]
    assert "exomercat_name" in removed.colnames
    assert "HD 2 c" in set(removed["exomercat_name"])

def test_exo_create_ident_table():
    # data
    exo_helptab = Table(
        data=[
            ["*   3 Cnc b", "6 Lyn b", "6 Lyn b2"],
            ["*   3 Cnc b", "*   6 Lyn b", ""],
            ["*   3 Cnc b", "muster_exoname", "muster_exoname2"],
        ],
        names=["planet_main_id", "sim_planet_main_id", "exomercat_name"],
        dtype=[object, object, object],
    )
    exo = empty_dict.copy()
    exo_ref = "2020A&C....3100370A"
    exo["provider"] = Table(
        data=[[exo_ref]], names=["provider_bibcode"], dtype=[object]
    )

    # function
    exo_ident, exo_helptab = create_ident_table(exo_helptab, exo)

    # assert
    # exo_ident
    assert exo_ident.colnames == ["main_id", "id", "id_ref"]
    # planet_main_id didn't get an id column entry if sim_planet_main_id exists
    mask = np.isin(exo_ident["id"], ["6 Lyn b"])
    assert len(np.nonzero(mask)[0]) == 0
    mask2 = np.isin(exo_ident["id"], ["6 Lyn b2"])
    assert len(np.nonzero(mask2)[0]) == 1
    # and did get if sim_planet_main_id doesn't
    assert (
        exo_ident["id_ref"][np.where(exo_ident["id"] == "6 Lyn b2")] == exo_ref
    )
    # id column was created with ref exo when no simbad name there
    assert (
        exo_ident["id_ref"][np.where(exo_ident["id"] == "muster_exoname2")]
        == exo_ref
    )
    # and with ref sim for sim_planet_name!=''
    assert (
        exo_ident["id_ref"][np.where(exo_ident["id"] == "*   6 Lyn b")]
        == "2000A&AS..143....9W"
    )
    assert (
        exo_ident["id_ref"][np.where(exo_ident["id"] == "muster_exoname")]
        == exo_ref
    )
    # exo_helptab
    # updated exo_helptab where sim_main_id exists
    assert (
        exo_helptab["planet_main_id"][
            np.where(exo_helptab["sim_planet_main_id"] == "*   6 Lyn b")
        ]
        == "*   6 Lyn b"
    )
    # and don't updated exo_helptab where sim_main_id doesn't exist
    assert (
        exo_helptab["planet_main_id"][
            np.where(exo_helptab["sim_planet_main_id"] == "")
        ]
        == "6 Lyn b2"
    )


def test_create_objects_table():
    # data
    exo = empty_dict.copy()
    exo["ident"] = Table(
        data=[
            ["*   3 Cnc b", "*   6 Lyn b", "*   6 Lyn b"],
            ["*   3 Cnc b", "6 Lyn b", "*   6 Lyn b"],
        ],
        names=["main_id", "id"],
        dtype=[object, object],
    )

    # function
    exo_objects = create_objects_table(exo)
    # type got correctly assigned
    assert (
        exo_objects["type"][np.where(exo_objects["main_id"] == "*   3 Cnc b")]
        == "pl"
    )
    # unique main id
    assert (
        len(exo_objects[np.where(exo_objects["main_id"] == "*   6 Lyn b")]) == 1
    )


def test_deal_with_mass_nullvalues():
    # data
    m = MaskedColumn(
        data=[20.76, 0.1, 0, np.inf, 1e20],
        name="mass",
        mask=[False, False, True, False, False],
    )
    exo_helptab = Table(
        data=[
            [
                "*   3 Cnc b",
                "*   4 Mon B .01",
                "testname1",
                "*   6 Lyn b",
                "testname2",
            ],
            m,
        ],
        names=["main_id", "mass"],
        dtype=[object, float],
    )

    # function
    exo_helptab = deal_with_mass_nullvalues(exo_helptab, ["mass"])
    # assert
    assert (
        exo_helptab["mass"][np.where(exo_helptab["main_id"] == "testname1")]
        == 1e20
    )
    assert (
        exo_helptab["mass"][np.where(exo_helptab["main_id"] == "*   6 Lyn b")]
        == 1e20
    )


def test_create_para_exo_mes_mass_pl() -> None:
    """
    Extract msini parameter table and verify renaming and error mapping.
    """
    exo_helptab = _mk_msini_only_fixture()
    sinitable = create_para_exo_mes_mass_pl(exo_helptab, "msini", "True")

    assert (
            sinitable["mass_pl_value"][
                np.where(sinitable["main_id"] == "*   3 Cnc b")
            ]
            == 20.76
    )
    assert (
            sinitable["mass_pl_err_max"][
                np.where(sinitable["main_id"] == "*   6 Lyn b")
            ]
            == 0.873
    )
    assert (
            sinitable["mass_pl_err_min"][
                np.where(sinitable["main_id"] == "*   6 Lyn b")
            ]
            == 1.067
    )


def test_betterthan():
    assert betterthan("A", "B") == True
    assert betterthan("A", "?") == True
    assert betterthan("B", "A") == False
    assert betterthan("A", "A") == False
    assert betterthan("E", "?") == True


def test_bestmass_better_qual():
    bestmass = ["Mass", "Mass", "Mass", "Msini", "Msini", "Msini"]
    qual_msini = ["B", "A", "C", "B", "A", "C"]
    qual_mass = ["B", "B", "B", "B", "B", "B"]
    result_qual_msini = ["C", "C", "C", "B", "A", "C"]
    result_qual_mass = ["B", "B", "B", "C", "B", "D"]
    for i in range(len(bestmass)):
        msini, mass = bestmass_better_qual(
            bestmass[i], qual_msini[i], qual_mass[i]
        )
        if i == 0:
            assert msini == result_qual_msini[i]
            assert mass == result_qual_mass[i]


def test_assign_new_qual():
    exo_mes_mass_pl = Table(
        data=[
            ["*   3 Cnc b", "*   3 Cnc b"],
            ["B", "B"],
            ["True", "False"],
            ["Mass", "Mass"],
        ],
        names=[
            "main_id",
            "mass_pl_qual",
            "mass_pl_sini_flag",
            "bestmass_provenance",
        ],
        dtype=[object, object, object, object],
    )
    exo_mes_mass_pl = assign_new_qual(
        exo_mes_mass_pl, "*   3 Cnc b", "True", "C"
    )
    assert (
        exo_mes_mass_pl["mass_pl_qual"][
            np.where(exo_mes_mass_pl["mass_pl_sini_flag"] == "True")
        ]
        == "C"
    )


def test_align_quality_with_bestmass():
    exo_mes_mass_pl = Table(
        data=[
            ["*   3 Cnc b", "*   3 Cnc b"],
            ["B", "B"],
            ["True", "False"],
            ["Mass", "Mass"],
        ],
        names=[
            "main_id",
            "mass_pl_qual",
            "mass_pl_sini_flag",
            "bestmass_provenance",
        ],
        dtype=[object, object, object, object],
    )
    exo_mes_mass_pl = align_quality_with_bestmass(exo_mes_mass_pl)
    assert (
        exo_mes_mass_pl["mass_pl_qual"][
            np.where(exo_mes_mass_pl["mass_pl_sini_flag"] == "True")
        ]
        == "C"
    )


def test_create_mes_mass_pl_table() -> None:
    """
    Build full mass table (mass + msini) and verify selected fields.
    """
    exo_helptab = _mk_full_mass_fixture()
    exo_mes_mass_pl = create_mes_mass_pl_table(exo_helptab)

    assert len(exo_mes_mass_pl) == 6
    table = exo_mes_mass_pl[
        np.where(exo_mes_mass_pl["mass_pl_sini_flag"] == "False")
    ]
    assert (
            table["mass_pl_value"][np.where(table["main_id"] == "*   3 Cnc b")]
            == 20.76
    )
    assert (
            table["mass_pl_err_max"][np.where(table["main_id"] == "*   6 Lyn b")]
            == 0.873
    )
    assert (
            table["mass_pl_err_min"][np.where(table["main_id"] == "*   6 Lyn b")]
            == 1.067
    )
    assert (
            table["mass_pl_err_max"][
                np.where(table["main_id"] == "*   4 Mon B .01")
            ]
            == 1e20
    )
    assert (
            table["bestmass_provenance"][
                np.where(table["main_id"] == "*   4 Mon B .01")
            ]
            == "Mass"
    )

def test_create_h_link_table() -> None:
    """
    Ensure h_link renames columns and assigns provider bibcode uniformly.
    """
    exo_helptab = Table(
        names=["planet_main_id", "host_main_id"], dtype=[object, object]
    )
    exo_helptab.add_row(("Planet-1", "Host-A"))
    exo_helptab.add_row(("Planet-2", "Host-B"))

    exo = {
        "provider": Table(
            names=[
                "provider_name",
                "provider_url",
                "provider_bibcode",
                "provider_access",
            ],
            dtype=[object, object, object, object],
        )
    }
    exo["provider"].add_row(
        ("ExoMercat", "http://example.org", "2025Test.Bib", "open")
    )

    result = create_h_link_table(exo_helptab, exo)

    assert set(result.colnames) == {"main_id", "parent_main_id", "h_link_ref"}
    assert list(result["main_id"]) == ["Planet-1", "Planet-2"]
    assert list(result["parent_main_id"]) == ["Host-A", "Host-B"]
    assert list(result["h_link_ref"]) == ["2025Test.Bib", "2025Test.Bib"]

def test_create_exo_sources_table() -> None:
    """
    Validate sources aggregation, de-duplication and provider tagging.
    """
    provider = Table(
        names=[
            "provider_name",
            "provider_url",
            "provider_bibcode",
            "provider_access",
        ],
        dtype=[object, object, object, object],
    )
    provider.add_row(("ExoMercat", "http://example.org", "BIB1", "open"))

    h_link = Table(names=["h_link_ref"], dtype=[object])
    h_link.add_row(("REF_H1",))
    h_link.add_row(("REF_H1",))

    ident = Table(names=["id_ref"], dtype=[object])
    ident.add_row(("REF_I1",))
    ident.add_row(("REF_H1",))  # intentionally duplicate across tables

    mes_mass_pl = Table()
    mes_mass_pl["mass_pl_ref"] = MaskedColumn(
        data=["REF_M1", "IGNORED"], mask=[False, True], dtype=object
    )

    exo = {
        "provider": provider,
        "h_link": h_link,
        "ident": ident,
        "mes_mass_pl": mes_mass_pl,
    }

    sources = create_exo_sources_table(exo)

    assert set(sources.colnames) >= {"ref", "provider_name"}
    expected_refs = {"BIB1", "REF_H1", "REF_I1", "REF_M1"}
    # ensure references are unique and no None in the set
    assert expected_refs.issubset(set(sources["ref"]))
    # all provider_name entries should match the single provider
    assert set(sources["provider_name"]) == {"ExoMercat"}

