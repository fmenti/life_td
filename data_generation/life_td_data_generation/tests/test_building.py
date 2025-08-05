from building import *
from astropy.table import Table, MaskedColumn, setdiff
from sdata import empty_dict_wit_columns, table_names
import pytest

def test_idsjoin_no_mask():
    """Test idsjoin with no masked values."""
    cat = Table(data=[['* 61 Cyg b|61 Cyg b', '', 'GCRV 13273|LTT 16180'],
                      ['61 Cyg b|markmuster', 'HD 201091|LSPM J2106+3844N', 'USNO-B1.0 1287-00443364']],
                names=['ids1', 'ids2'],
                dtype=[object, object])

    result = idsjoin(cat, 'ids1', 'ids2')

    # Assert no empty values
    for ID in result['ids']:
        assert ID.split('|').count('') == 0

    # Assert unique merged identifiers
    assert result['ids'][0].split('|').count('61 Cyg b') == 1
    assert result['ids'][0].split('|').count('* 61 Cyg b') == 1
    assert result['ids'][0].split('|').count('markmuster') == 1


def test_idsjoin_with_mask():
    """Test idsjoin with masked values."""
    a = MaskedColumn(data=['* 61 Cyg b|61 Cyg b', '', 'GCRV 13273|LTT 16180'], name='ids1', mask=[True, False, False])
    cat = Table(data=[a,
                      ['61 Cyg b|markmuster', 'HD 201091|LSPM J2106+3844N', 'USNO-B1.0 1287-00443364']],
                names=['ids1', 'ids2'],
                dtype=[object, object])

    result = idsjoin(cat, 'ids1', 'ids2')

    # Assert no empty values
    for ID in result['ids']:
        assert ID.split('|').count('') == 0

    # Assert unique merged identifiers
    assert result['ids'][0].split('|').count('61 Cyg b') == 1
    assert result['ids'][0].split('|').count('markmuster') == 1


def test_idsjoin_duplicates():
    """Test idsjoin with duplicate values across columns."""
    cat = Table(data=[['ID1|ID2', 'ID4|ID5'],
                      ['ID2|ID3', 'ID4|ID6']],
                names=['ids1', 'ids2'],
                dtype=[object, object])

    result = idsjoin(cat, 'ids1', 'ids2')

    # Assert duplicates are removed
    assert result['ids'][0] == 'ID1|ID2|ID3'
    assert result['ids'][1] == 'ID4|ID5|ID6'


def test_idsjoin_empty_columns():
    """Test idsjoin with empty and masked columns."""
    cat = Table(data=[['', ''],
                      ['', '']],
                names=['ids1', 'ids2'],
                dtype=[object, object])

    result = idsjoin(cat, 'ids1', 'ids2')

    # Assert resulting 'ids' column is empty
    assert all(val == '' for val in result['ids'])


def test_idsjoin_large_table():
    """Test idsjoin on a large table."""
    n_rows = 10_000
    ids1 = ['ID{}'.format(i) for i in range(n_rows)]
    ids2 = ['ID{}'.format(i + 1) for i in range(n_rows)]

    cat = Table({'ids1': ids1, 'ids2': ids2})
    result = idsjoin(cat, 'ids1', 'ids2')

    # Assert all rows processed correctly
    assert len(result) == n_rows
    for i, id_val in enumerate(result['ids']):
        assert f'ID{i}' in id_val
        assert f'ID{i + 1}' in id_val


def test_best_para_id():
    mes_table = Table(data=[[1, 1, 2, 2, 3, 3, 4, 4],
                            ['sim_id_obj1', 'irgendein_id_obj1',
                             'exomercat_id_obj2', 'gaia_id_obj2',
                             'gaia_id_obj3', 'sim_id_obj3',
                             'sim_id_obj4', 'sim_id_obj4'],
                            ['2000A&AS..143....9W', '1925AnHar.100...17C',
                             '2020A&C....3100370A', '2016A&A...595A...1G',
                             '2016A&A...595A...1G', '2000A&AS..143....9W',
                             '2016A&A...595A...1G', '2000A&AS..143....9W']],
                      names=['object_idref', 'id', 'id_ref'],
                      dtype=[int, object, object])
    best_para_table = best_para_id(mes_table)

    # non provider reverences don't get included
    assert 'irgendein_id_obj1' not in best_para_table['id']
    # ids of the same object are all unique
    assert len(best_para_table[np.where(best_para_table['object_idref'] == 4)]) == 1


def test_best_para_of_id():
    mes_table = Table(data=[[1, 1, 2, 2, 3, 3, 4, 4],
                            ['sim_id_obj1', 'irgendein_id_obj1',
                             'exomercat_id_obj2', 'gaia_id_obj2',
                             'gaia_id_obj3', 'sim_id_obj3',
                             'sim_id_obj4', 'sim_id_obj4'],
                            ['2000A&AS..143....9W', '1925AnHar.100...17C',
                             '2020A&C....3100370A', '2016A&A...595A...1G',
                             '2016A&A...595A...1G', '2000A&AS..143....9W',
                             '2016A&A...595A...1G', '2000A&AS..143....9W']],
                      names=['object_idref', 'id', 'id_ref'],
                      dtype=[int, object, object])
    best_para_table = best_para('id', mes_table)
    assert type(best_para_table) == type(Table())


def test_best_para_membership():
    mes_table = Table(data=[[1, 1, 3, 3, 1, 1],
                            [2, 2, 4, 4, 5, 5],
                            [100, 50, 999999, 20, 999999, 999999]],
                      names=['child_object_idref', 'parent_object_idref', 'membership'],
                      dtype=[int, int, int])
    best_para_table = best_para_membership(mes_table)
    wanted_table = Table(data=[[1, 3, 1],
                               [2, 4, 5],
                               [100, 20, 999999]],
                         names=['child_object_idref', 'parent_object_idref', 'membership'],
                         dtype=[int, int, int])
    assert len(setdiff(wanted_table, mes_table)) == 0


def test_best_para():
    """
    Test the best_para function with different parameter types and edge cases.
    """
    # Test Case 1: Temperature measurements (multiple values, same quality)
    mes_table_teff = Table(data=[
        ['*  61 Cyg A', '*  61 Cyg A', '*  61 Cyg A'],
        [4192.0, 4353.743, 4440.0],
        ['B', 'B', 'D'],
        [1e+20, 1e+20, 1e+20],
        [3, 2, 3]
    ], names=[
        'main_id', 'teff_st_value', 'teff_st_qual',
        'teff_st_err', 'teff_st_source_idref'
    ], dtype=[object, float, object, float, int])

    result_teff = best_para('teff_st', mes_table_teff)
    assert len(result_teff) == 1
    assert result_teff['teff_st_value'][0] == 4192.0  # First occurrence wins
    assert result_teff['teff_st_qual'][0] == 'B'
    assert result_teff['teff_st_source_idref'][0] == 3

    # Test Case 2: Planet mass measurements (complex parameter)
    mes_table_mass = Table(data=[
        ['* 61 Cyg b', '* 61 Cyg b', '* 61 Cyg b'],
        [10, 1.5, 10],
        ['B', 'B', 'D'],
        [1e+20, 0.5, 1e+20],
        [1e+20, 0.5, 1e+20],
        ['?', '?', '?'],
        ['True', 'False', 'True'],
        [1, 2, 3]
    ], names=[
        'main_id', 'mass_pl_value', 'mass_pl_qual', 'mass_pl_err_max',
        'mass_pl_err_min', 'mass_pl_rel', 'mass_pl_sini_flag',
        'mass_pl_source_idref'
    ], dtype=[object, float, object, float, float, object, object, int])

    result_mass = best_para('mass_pl', mes_table_mass)
    assert len(result_mass) == 1
    assert result_mass['mass_pl_value'][0] == 10
    assert result_mass['mass_pl_qual'][0] == 'B'
    assert result_mass['mass_pl_sini_flag'][0] == 'True'

    # Test Case 3: Multiple objects with different qualities
    test_data = Table({
        'main_id': ['star1', 'star1', 'star2', 'star2'],
        'teff_st_value': [5000, 5100, 6000, 6100],
        'teff_st_err': [100, 150, 100, 150],
        'teff_st_qual': ['C', 'A', 'B', 'C'],
        'teff_st_source_idref': [1, 2, 4, 5]
    })

    result_multi = best_para('teff_st', test_data)
    assert len(result_multi) == 2
    star1_row = result_multi[result_multi['main_id'] == 'star1']
    assert star1_row['teff_st_qual'][0] == 'A'
    star2_row = result_multi[result_multi['main_id'] == 'star2']
    assert star2_row['teff_st_qual'][0] == 'B'

    # Test Case 4: Binary parameter handling
    binary_data = Table({
        'main_id': ['star1', 'star1', 'star2'],
        'binary_flag': ['True', 'True', 'False'],
        'binary_qual': ['B', 'A', 'C'],
        'binary_source_idref': [1, 2, 3]
    })

    binary_result = best_para('binary', binary_data)
    assert len(binary_result) == 2
    mask = binary_result['main_id'] == 'star1'
    star1_rows = binary_result[mask]
    star1_qual = star1_rows['binary_qual'][0]
    assert star1_qual == 'A'


def test_assign_type_with_masked_constant():
    # Create a mock Table
    cat = Table(
        {
            'type_1': ['A', 'B', 'C'],
            'type_2': [np.ma.core.MaskedConstant(), 'None', 'D'],
            'type': [''] * 3,  # Initialize 'type' column as empty
        }
    )

    # Test the function for each row in the table
    assert assign_type(cat, 0) == 'A'  # type_2 is masked, fallback to type_1
    assert assign_type(cat, 1) == 'B'  # type_2 is 'None', fallback to type_1
    assert assign_type(cat, 2) == 'D'  # type_2 is valid, use it


def test_assign_type_with_no_masked_constant():
    # Create a mock Table
    cat = Table(
        {
            'type_1': ['X'],
            'type_2': ['Y'],
            'type': [''],  # Initialize 'type' column as empty
        }
    )

    # Test the function for the first (and only) row in the table
    assert assign_type(cat, 0) == 'Y'  # type_2 is valid, use it



def test_assign_source_idref():
    # Data setup
    # Create a catalog table with reference columns
    cat = Table({
        'main_id': ['star1', 'star2', 'star3'],
        'temp_ref': ['ref1', 'ref2', 'None'],
        'temp_value': MaskedColumn([100, 200, 300], mask=[False, False, True]),
        'mass_ref': ['ref1', 'None', 'ref2']
    })

    # Create a sources table
    sources = Table({
        'ref': ['ref1', 'ref2'],
        'source_id': [1, 2],
        'provider_name': ['provider1', 'provider1']
    })

    # Parameters to process
    paras = ['temp', 'mass']
    provider = 'provider1'

    # Execute function
    result = assign_source_idref(cat, sources, paras, provider)

    # Assertions
    # Check if source_idref columns were created
    assert 'temp_source_idref' in result.colnames
    assert 'mass_source_idref' in result.colnames

    # Check if the correct source_ids are assigned for each reference
    temp_ref1_rows = result[result['temp_ref'] == 'ref1']
    temp_ref2_rows = result[result['temp_ref'] == 'ref2']
    assert all(temp_ref1_rows['temp_source_idref'] == 1)
    assert all(temp_ref2_rows['temp_source_idref'] == 2)

    mass_ref1_rows = result[result['mass_ref'] == 'ref1']
    mass_ref2_rows = result[result['mass_ref'] == 'ref2']
    assert all(mass_ref1_rows['mass_source_idref'] == 1)
    assert all(mass_ref2_rows['mass_source_idref'] == 2)

    # Check if masked values get source_id 999999
    masked_temp_rows = result[result['temp_value'].mask]
    assert all(masked_temp_rows['temp_source_idref'] == 999999)

    # Check number of rows remained the same
    assert len(result) == len(cat)



