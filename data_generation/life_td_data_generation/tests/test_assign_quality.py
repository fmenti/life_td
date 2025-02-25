from provider.assign_quality import *
from astropy.table import Table
import pytest
import numpy as np
from unittest.mock import patch

#import pytest


def test_teff_st_spec_assign_quality():
    gaia_mes_teff_st_spec = Table(data=[
        ['3 Cnc', 'TIC 72090501', '6 Lyn'],
        [3060, 4250, 4000],
        ['00000000000000000000000000000990099999995',
         '00000000000000000000000009900001099999995',
         '09009090000090090900900009900991099129995']],
        names=['main_id', 'teff_gspspec', 'flags_gspspec'],
        dtype=[object, float, object])
    gaia_mes_teff_st_spec = assign_quality(gaia_mes_teff_st_spec)

    assert gaia_mes_teff_st_spec['teff_st_qual'][0] == 'B'
    assert gaia_mes_teff_st_spec['teff_st_qual'][2] == 'C'


@pytest.fixture
def table_with_binary_flag():
    """Fixture for test table with binary_flag column."""
    return {'binary_flag': ['True', 'False', 'True']}


@pytest.fixture
def table_with_column_data():
    """Fixture for a generic table with additional columns."""
    return {'some_column': [1, 2, 3]}


def test_teff_st_spec_special_mode(table_with_column_data):
    """Test 'teff_st_spec' special mode."""
    with patch('your_module.teff_st_spec_assign_quality') as mock_func:
        mock_func.return_value = table_with_column_data  # Mock the returned table
        result = assign_quality(table_with_column_data, special_mode='teff_st_spec')
        mock_func.assert_called_once_with(table_with_column_data)
        assert result == table_with_column_data


def test_exo_special_mode(table_with_column_data):
    """Test 'exo' special mode."""
    with patch('your_module.exo_assign_quality') as mock_func:
        mock_func.return_value = table_with_column_data  # Mock the returned table
        result = assign_quality(table_with_column_data, special_mode='exo')
        mock_func.assert_called_once_with(table_with_column_data)
        assert result == table_with_column_data


def test_gaia_binary_special_mode(table_with_binary_flag):
    """Test 'gaia_binary' special mode."""
    result = assign_quality(table_with_binary_flag, column='quality', special_mode='gaia_binary')
    assert result['quality'] == ['B', 'E', 'B']  # Verify correct 'B' or 'E' assignment


def test_assign_quality_wds_sep1():
    # Mock table
    table = Table({
        'sep_ang_obs_date': [1, 2, np.ma.core.MaskedConstant()],
    })

    result = assign_quality(table, column='quality', special_mode='wds_sep1')
    assert list(result['quality']) == ['C', 'C', 'E']  # Check the results


def test_assign_quality_wds_sep2():
    # Mock table
    table = Table({
        'sep_ang_obs_date': [1, np.ma.core.MaskedConstant(), 5],
    })

    result = assign_quality(table, column='quality', special_mode='wds_sep2')
    assert list(result['quality']) == ['B', 'E', 'B']  # Check the results



def test_default_fallback_logic(table_with_column_data):
    """Test default fallback logic for special_mode and column behavior."""
    # Test 'coo_gal_qual'
    result_coo_gal = assign_quality(table_with_column_data, column='coo_gal_qual')
    assert result_coo_gal['coo_gal_qual'] == ['?' for _ in range(len(table_with_column_data['some_column']))]

    # Test 'teff_st_phot'
    result_teff = assign_quality(table_with_column_data, column='quality', special_mode='teff_st_phot')
    assert result_teff['quality'] == ['B', 'B', 'B']

    # Test 'model'
    result_model = assign_quality(table_with_column_data, column='quality', special_mode='model')
    assert result_model['quality'] == ['C', 'C', 'C']

    # Test 'sim_binary'
    result_sim = assign_quality(table_with_column_data, column='quality', special_mode='sim_binary')
    assert result_sim['quality'] == ['D', 'D', 'D']

    # Test unknown special mode or column
    result_unknown = assign_quality(table_with_column_data, column='unknown_column', special_mode='unknown_mode')
    assert result_unknown['unknown_column'] == ['?' for _ in range(len(table_with_column_data['some_column']))]
