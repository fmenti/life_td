from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table, setdiff, MaskedColumn
from provider.assign_quality_funcs import teff_st_spec_assign_quality, assign_quality_elementwise, assign_quality, exo_assign_quality


def test_teff_st_spec_assign_quality():
    gaia_mes_teff_st_spec = Table(
        data=[
            ["3 Cnc", "TIC 72090501", "6 Lyn", "HD 166348"],
            [3060, 4250, 4000, 4108.0],
            [
                "00000000000000000000000000000990099999995",
                "00000000000000000000000009900001099999995",
                "09009090000090090900900009900991099129995",
                "00100131000019999999999999999999999999999",
            ],
        ],
        names=["main_id", "teff_gspspec", "flags_gspspec"],
        dtype=[object, float, object],
    )
    gaia_mes_teff_st_spec = teff_st_spec_assign_quality(gaia_mes_teff_st_spec)

    assert gaia_mes_teff_st_spec["teff_st_qual"][0] == "B"
    assert gaia_mes_teff_st_spec["teff_st_qual"][2] == "C"
    assert gaia_mes_teff_st_spec["teff_st_qual"][3] == "D"


def test_assign_quality_elementwise():
    # data
    exo_helptab = Table(
        data=[[1e20, 1e20, 1], [1e20, 1, 1]],
        names=["mass_max", "mass_min"],
        dtype=[float, float],
    )
    exo_helptab["mass_pl_qual"] = MaskedColumn(
        dtype=object, length=len(exo_helptab)
    )

    # function
    qual_b = assign_quality_elementwise(exo_helptab, "mass", 2)
    qual_c = assign_quality_elementwise(exo_helptab, "mass", 1)
    qual_d = assign_quality_elementwise(exo_helptab, "mass", 0)

    # assert
    assert qual_b == "B"
    assert qual_c == "C"
    assert qual_d == "D"


def test_exo_assign_quality():
    # data
    exo_helptab = Table(
        data=[[1e20, 1e20, 1], [1e20, 1, 1], [1e20, 1e20, 1], [1e20, 1, 1]],
        names=["mass_max", "mass_min", "msini_max", "msini_min"],
        dtype=[float, float, float, float],
    )
    # function
    exo_helptab = exo_assign_quality(exo_helptab)
    # assert
    assert exo_helptab["mass_pl_qual"][0] == "D"
    assert exo_helptab["mass_pl_qual"][1] == "C"
    assert exo_helptab["mass_pl_qual"][2] == "B"
    assert exo_helptab["msini_pl_qual"][0] == "D"
    assert exo_helptab["msini_pl_qual"][1] == "C"
    assert exo_helptab["msini_pl_qual"][2] == "B"


@pytest.fixture
def table_with_binary_flag():
    """Fixture for test table with binary_flag column."""
    return Table({"binary_flag": ["True", "False", "True"]})


@pytest.fixture
def table_with_column_data():
    """Fixture for a generic table with additional columns."""
    return Table({"some_column": [1, 2, 3]})


def test_teff_st_spec_special_mode(table_with_column_data):
    """Test 'teff_st_spec' special mode."""
    with patch(
        "provider.assign_quality_funcs.teff_st_spec_assign_quality"
    ) as mock_func:
        mock_func.return_value = (
            table_with_column_data  # Mock the returned table
        )
        result = assign_quality(
            table_with_column_data, special_mode="teff_st_spec"
        )
        mock_func.assert_called_once_with(table_with_column_data)
        assert len(setdiff(result, table_with_column_data)) == 0


def test_exo_special_mode(table_with_column_data):
    """Test 'exo' special mode."""
    with patch("provider.assign_quality_funcs.exo_assign_quality") as mock_func:
        mock_func.return_value = (
            table_with_column_data  # Mock the returned table
        )
        result = assign_quality(table_with_column_data, special_mode="exo")
        mock_func.assert_called_once_with(table_with_column_data)
        assert len(setdiff(result, table_with_column_data)) == 0


def test_gaia_binary_special_mode(table_with_binary_flag):
    """Test 'gaia_binary' special mode."""
    result = assign_quality(
        table_with_binary_flag, column="quality", special_mode="gaia_binary"
    )
    assert list(result["quality"]) == [
        "B",
        "E",
        "B",
    ]  # Verify correct 'B' or 'E' assignment


def test_assign_quality_wds_sep1():
    # Mock table
    table = Table(
        {
            "sep_ang_obs_date": [1, 2, np.ma.core.MaskedConstant()],
        }
    )

    result = assign_quality(table, column="quality", special_mode="wds_sep1")
    assert list(result["quality"]) == ["C", "C", "E"]  # Check the results


def test_assign_quality_wds_sep2():
    # Mock table
    table = Table(
        {
            "sep_ang_obs_date": [1, np.ma.core.MaskedConstant(), 5],
        }
    )

    result = assign_quality(table, column="quality", special_mode="wds_sep2")
    assert list(result["quality"]) == ["B", "E", "B"]  # Check the results


def test_default_fallback_logic(table_with_column_data):
    """Test default fallback logic for special_mode and column behavior."""
    # Test 'coo_gal_qual'
    result_coo_gal = assign_quality(
        table_with_column_data, column="coo_gal_qual"
    )
    assert result_coo_gal["coo_gal_qual"].tolist() == [
        "?" for _ in range(len(table_with_column_data))
    ]

    # Test 'teff_st_phot'
    result_teff = assign_quality(
        table_with_column_data, column="teff_st_qual", special_mode="teff_st_phot"
    )
    assert result_teff["teff_st_qual"].tolist() == [
        "?",
        "?",
        "?",
    ]  # Replace with expected values if different

    # Test 'model'
    result_model = assign_quality(
        table_with_column_data, column="quality", special_mode="model"
    )
    assert result_model["quality"].tolist() == [
        "C",
        "C",
        "C",
    ]  # Replace with expected values if different

    # Test 'sim_binary'
    result_sim = assign_quality(
        table_with_column_data, column="quality", special_mode="sim_binary"
    )
    assert result_sim["quality"].tolist() == [
        "D",
        "D",
        "D",
    ]  # Replace with expected values if different

    # Test unknown special mode or column
    result_unknown = assign_quality(
        table_with_column_data,
        column="unknown_column",
        special_mode="unknown_mode",
    )
    assert result_unknown["unknown_column"].tolist() == [
        "?" for _ in range(len(table_with_column_data))
    ]
