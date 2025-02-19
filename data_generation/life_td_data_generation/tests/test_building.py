from building import *
from astropy.table import Table, MaskedColumn, setdiff


def test_idsjoin_no_mask():
    cat = Table(data=[['* 61 Cyg b|61 Cyg b', '', 'GCRV 13273|LTT 16180'],
                      ['61 Cyg b|markmuster', 'HD 201091|LSPM J2106+3844N', 'USNO-B1.0 1287-00443364']],
                names=['ids1', 'ids2'],
                dtype=[object, object])
    cat = idsjoin(cat, 'ids1', 'ids2')

    for ID in cat['ids']:
        assert ID.split('|').count('') == 0

    assert cat['ids'][0].split('|').count('61 Cyg b') == 1
    assert cat['ids'][0].split('|').count('* 61 Cyg b') == 1
    assert cat['ids'][0].split('|').count('markmuster') == 1


def test_idsjoin_with_mask():
    a = MaskedColumn(data=['* 61 Cyg b|61 Cyg b', '', 'GCRV 13273|LTT 16180'], name='ids1', mask=[True, False, False])
    cat = Table(data=[a,
                      ['61 Cyg b|markmuster', 'HD 201091|LSPM J2106+3844N', 'USNO-B1.0 1287-00443364']],
                names=['ids1', 'ids2'],
                dtype=[object, object])
    cat = idsjoin(cat, 'ids1', 'ids2')

    for ID in cat['ids']:
        assert ID.split('|').count('') == 0

    assert cat['ids'][0].split('|').count('61 Cyg b') == 1
    assert cat['ids'][0].split('|').count('markmuster') == 1


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

    #non provider reverences don't get included
    assert 'irgendein_id_obj1' not in best_para_table['id']
    #ids of the same object are all unique
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
    mes_table = Table(data=[['*  61 Cyg A', '*  61 Cyg A', '*  61 Cyg A'],
                            [4192.0, 4353.743, 4440.0],
                            ['B', 'B', 'D'],
                            [1e+20, 1e+20, 1e+20],
                            [3, 2, 3]],
                      names=['main_id', 'teff_st_value', 'teff_st_qual',
                             'teff_st_err', 'teff_st_source_idref'],
                      dtype=[object, float, object, float, int])
    best_para_table = best_para('teff_st', mes_table)
    wanted_table = Table(data=[['*  61 Cyg A'],
                               [4192.0],
                               ['B'],
                               [1e+20],
                               [3]],
                         names=['main_id', 'teff_st_value', 'teff_st_qual',
                                'teff_st_err', 'teff_st_source_idref'],
                         dtype=[object, float, object, float, int])
    assert len(setdiff(wanted_table, mes_table)) == 0

#def test_building():
#input_data=
#building(prov_tables_list,table_names,list_of_tables)
#    assert False
