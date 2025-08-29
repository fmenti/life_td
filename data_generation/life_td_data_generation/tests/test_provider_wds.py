from provider.wds import assign_names
from astropy.table import Table
from provider.wds import query


def test_wds_columns_queryable():
    wds_provider = Table()
    wds_provider["provider_url"] = [
        "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"
    ]
    adql_query = [
        """SELECT TOP 10
                    WDS as wds_name, Comp as wds_comp,
                    sep1 as wds_sep1, sep2 as wds_sep2,
                    Obs1 as wds_obs1, Obs2 as wds_obs2
                    FROM "B/wds/wds" """
    ]
    wds = query(wds_provider["provider_url"][0], adql_query[0])
    assert wds.colnames == [
        "wds_name",
        "wds_comp",
        "wds_sep1",
        "wds_sep2",
        "wds_obs1",
        "wds_obs2",
    ]


def test_assign_names_trivial_binaries():
    """
    Test assign_names on trivial binaries (empty 'wds_comp').
    """
    # Input table with trivial binaries
    wds_helptab = Table(
        data=[
            ["12345+6789"],  # wds_name
            [""],  # wds_comp
            [None],  # system_name
            [None],  # primary
            [None],  # secondary
        ],
        names=["wds_name", "wds_comp", "system_name", "primary", "secondary"],
        dtype=[object, object, object, object, object],
    )

    # Apply function
    updated_table = assign_names(wds_helptab)

    # Check output
    assert updated_table["system_name"][0] == "WDS J12345+6789AB"
    assert updated_table["primary"][0] == "WDS J12345+6789A"
    assert updated_table["secondary"][0] == "WDS J12345+6789B"


def test_assign_names_higher_order_multiples_two_components():
    """
    Test assign_names on higher-order multiples with two components ('wds_comp' length = 2).
    """
    # Input table with two components in 'wds_comp'
    wds_helptab = Table(
        data=[
            ["12345+6789"],  # wds_name
            ["CD"],  # wds_comp
            [None],  # system_name
            [None],  # primary
            [None],  # secondary
        ],
        names=["wds_name", "wds_comp", "system_name", "primary", "secondary"],
        dtype=[object, object, object, object, object],
    )

    # Apply function
    updated_table = assign_names(wds_helptab)

    # Check output
    assert updated_table["system_name"][0] == "WDS J12345+6789CD"
    assert updated_table["primary"][0] == "WDS J12345+6789C"
    assert updated_table["secondary"][0] == "WDS J12345+6789D"


def test_assign_names_higher_order_multiples_multiple_components():
    """
    Test assign_names on higher-order multiples with multiple components ('wds_comp' contains a comma).
    """
    # Input table with multiple components in 'wds_comp'
    wds_helptab = Table(
        data=[
            ["12345+6789"],  # wds_name
            ["C,D"],  # wds_comp
            [None],  # system_name
            [None],  # primary
            [None],  # secondary
        ],
        names=["wds_name", "wds_comp", "system_name", "primary", "secondary"],
        dtype=[object, object, object, object, object],
    )

    # Apply function
    updated_table = assign_names(wds_helptab)

    # Check output
    assert updated_table["system_name"][0] == "WDS J12345+6789C,D"
    assert updated_table["primary"][0] == "WDS J12345+6789C"
    assert updated_table["secondary"][0] == "WDS J12345+6789D"


def test_assign_names_mixed_rows():
    """
    Test assign_names with a mix of trivial binaries and higher-order multiples.
    """
    # Input table with mixed rows
    wds_helptab = Table(
        data=[
            ["12345+6789", "98765-4321", "44444+1111"],  # wds_name
            ["", "AB", "X,Y"],  # wds_comp
            [None, None, None],  # system_name
            [None, None, None],  # primary
            [None, None, None],  # secondary
        ],
        names=["wds_name", "wds_comp", "system_name", "primary", "secondary"],
        dtype=[object, object, object, object, object],
    )

    # Apply function
    updated_table = assign_names(wds_helptab)

    # Check outputs
    # Row 1
    assert updated_table["system_name"][0] == "WDS J12345+6789AB"
    assert updated_table["primary"][0] == "WDS J12345+6789A"
    assert updated_table["secondary"][0] == "WDS J12345+6789B"

    # Row 2
    assert updated_table["system_name"][1] == "WDS J98765-4321AB"
    assert updated_table["primary"][1] == "WDS J98765-4321A"
    assert updated_table["secondary"][1] == "WDS J98765-4321B"

    # Row 3
    assert updated_table["system_name"][2] == "WDS J44444+1111X,Y"
    assert updated_table["primary"][2] == "WDS J44444+1111X"
    assert updated_table["secondary"][2] == "WDS J44444+1111Y"
