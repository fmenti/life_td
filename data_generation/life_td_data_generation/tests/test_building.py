from building import *
from astropy.table import Table, MaskedColumn, setdiff
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
    mes_table_teff = Table(data=[['*  61 Cyg A', '*  61 Cyg A', '*  61 Cyg A'],
                            [4192.0, 4353.743, 4440.0],
                            ['B', 'B', 'D'],
                            [1e+20, 1e+20, 1e+20],
                            [3, 2, 3]],
                      names=['main_id', 'teff_st_value', 'teff_st_qual',
                             'teff_st_err', 'teff_st_source_idref'],
                      dtype=[object, float, object, float, int])
    best_para_table_teff = best_para('teff_st', mes_table_teff)
    wanted_table_teff = Table(data=[['*  61 Cyg A'],
                               [4192.0],
                               ['B'],
                               [1e+20],
                               [3]],
                         names=['main_id', 'teff_st_value', 'teff_st_qual',
                                'teff_st_err', 'teff_st_source_idref'],
                         dtype=[object, float, object, float, int])
    assert len(setdiff(wanted_table_teff, best_para_table_teff)) == 0

    mes_table_mass = Table(data=[['* 61 Cyg b', '* 61 Cyg b', '* 61 Cyg b'],
                            [10, 1.5, 10],
                            ['B', 'B', 'D'],
                            [1e+20, 0.5, 1e+20],
                            [1e+20, 0.5, 1e+20],
                            ['?', '?', '?'],
                            ['True', 'False', 'True'],
                            [1, 2, 3]],
                      names=['main_id', 'mass_pl_value', 'mass_pl_qual',
                             'mass_pl_err_max','mass_pl_err_min',
                             'mass_pl_rel','mass_pl_sini_flag','mass_pl_source_idref'],
                      dtype=[object, float, object, float,float,object, object,int])
    best_para_table_mass = best_para('mass_pl', mes_table_mass)
    wanted_table_mass = Table(data=[['* 61 Cyg b'],
                               [10],
                               ['B'],
                               [1e+20],
                                [1e+20],
                                ['?'],
                                ['True'],
                               [1]],
                        names=['main_id', 'mass_pl_value', 'mass_pl_qual',
                                'mass_pl_err_max', 'mass_pl_err_min',
                                'mass_pl_rel','mass_pl_sini_flag','mass_pl_source_idref'],
                         dtype=[object, float, object, float, float,object,object,int])
    assert len(setdiff(wanted_table_mass, best_para_table_mass)) == 0


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



