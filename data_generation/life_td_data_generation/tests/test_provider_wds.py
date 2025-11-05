from astropy.table import Table
from provider.wds import query, create_objects_table, assign_names
import numpy as np


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


class TestCreateObjectsTable:
    """Tests for the create_objects_table function."""

    def test_basic_functionality(self):
        """Test that the function creates a proper objects table with correct columns."""
        # Setup: Create mock ident and h_link tables
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B", "System AB"],
                ["Star B", "Star B", "System AB"],
            ]
        )

        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B"],
                ["System AB", "System AB"],
            ]
        )

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        wds_helptab = Table()  # Not used in this function
        test_objects = []

        # Execute
        result = create_objects_table(wds_helptab, wds, test_objects)

        # Assert
        assert "main_id" in result.colnames
        assert "ids" in result.colnames
        assert "type" in result.colnames
        assert len(result) == 3

    def test_system_classification(self):
        """Test that objects with children are classified as systems ('sy')."""
        # Setup
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B", "System AB"],
                ["Star B", "Star B", "System AB"],
            ]
        )

        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B"],
                ["System AB", "System AB"],
            ]
        )

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        # Execute
        result = create_objects_table(Table(), wds, [])

        # Assert: System AB has children, so it should be type 'sy'
        system_row = result[result["main_id"] == "System AB"]
        assert len(system_row) == 1
        assert system_row["type"][0] == "sy"

    def test_star_classification(self):
        """Test that objects without children are classified as stars ('st')."""
        # Setup
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B", "System AB"],
                ["Star A", "Star B", "System AB"],
            ]
        )

        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B"],
                ["System AB", "System AB"],
            ]
        )

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        # Execute
        result = create_objects_table(Table(), wds, [])

        # Assert: Star A and Star B have no children, so they should be type 'st'
        star_a_row = result[result["main_id"] == "Star A"]
        star_b_row = result[result["main_id"] == "Star B"]

        assert len(star_a_row) == 1
        assert star_a_row["type"][0] == "st"

        assert len(star_b_row) == 1
        assert star_b_row["type"][0] == "st"

    def test_all_stars_no_systems(self):
        """Test case where all objects are stars with no hierarchical structure."""
        # Setup
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["Star 1", "Star 2", "Star 3"],
                ["Star 2", "Star 2", "Star 3"],
            ]
        )

        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
        )  # Empty h_link table - no parent-child relationships

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        # Execute
        result = create_objects_table(Table(), wds, [])

        # Assert: All should be classified as stars
        assert len(result) == 3
        assert all(result["type"] == "st")

    def test_all_systems(self):
        """Test case where all objects in ident table are parent systems."""
        # Setup
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["System 1", "System 2"],
                ["System 1", "System 2"],
            ]
        )

        # Mock h_link where both systems have children (not in ident table)
        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
            data=[
                ["Child A", "Child B" ],
                ["System 1", "System 2"],
            ]
        )

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        # Execute
        result = create_objects_table(Table(), wds, [])

        # Assert: All should be classified as systems
        assert len(result) == 2
        assert all(result["type"] == "sy")

    def test_complex_hierarchy(self):
        """Test a more complex hierarchical structure with multiple levels."""
        # Setup: Triple system with binary subsystem
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B", "Star C", "Binary AB", "Triple ABC"],
                ["Star A", "Star B", "Star C", "Binary AB", "Triple ABC"],
            ]
        )

        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B","Binary AB","Star C"],
                ["Binary AB", "Binary AB","Triple ABC","Triple ABC"],
            ]
        )

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        # Execute
        result = create_objects_table(Table(), wds, [])

        # Assert
        assert len(result) == 5

        # Individual stars should be type 'st'
        for star in ["Star A", "Star B", "Star C"]:
            star_row = result[result["main_id"] == star]
            assert star_row["type"][
                       0] == "st", f"{star} should be classified as star"

        # Systems should be type 'sy'
        for system in ["Binary AB", "Triple ABC"]:
            system_row = result[result["main_id"] == system]
            assert system_row["type"][
                       0] == "sy", f"{system} should be classified as system"

    def test_with_test_objects(self, capsys):
        """Test that test_objects are properly reported when present."""
        # Setup
        wds_ident = Table(
            names=["main_id", "id"],
            dtype=[object, object],
            data=[
                ["Star A", "Star B"],
                ["Star A", "Star B"],
            ]
        )

        wds_h_link = Table(
            names=["main_id", "parent_main_id"],
            dtype=[object, object],
        )

        wds = {
            "ident": wds_ident,
            "h_link": wds_h_link,
        }

        test_objects = np.array(["Star A", "Star X"])

        # Execute
        result = create_objects_table(Table(), wds, test_objects)

        # Assert
        captured = capsys.readouterr()
        assert "Star A" in captured.out
        assert "Star X" not in captured.out  # Star X is not in the result
