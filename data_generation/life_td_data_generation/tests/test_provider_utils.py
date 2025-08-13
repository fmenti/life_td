from provider.utils import *


#import pytest


def test_create_provider_table_date_given():
    gk_provider = create_provider_table('Grant Kennedy Disks',
                                        'http://drgmk.com/sdb/',
                                        'priv. comm.', '2024-02-09')
    assert gk_provider['provider_name'] == 'Grant Kennedy Disks'
    assert gk_provider['provider_url'] == 'http://drgmk.com/sdb/'
    assert gk_provider['provider_bibcode'] == 'priv. comm.'
    assert gk_provider['provider_access'] == '2024-02-09'


def test_create_provider_table_no_date_given():
    gk_provider = create_provider_table('Grant Kennedy Disks',
                                        'http://drgmk.com/sdb/',
                                        'priv. comm.')
    assert gk_provider['provider_access'] == datetime.now().strftime('%Y-%m-%d')


#@pytest.mark.dependency()
def test_fetch_main_id():
    #data
    id_cat = Table(data=[['HD 166', '6 Cet']],
                names=['sdb_host_main_id'],
                dtype=[object])
    oid_cat = Table(data=[[1800872, 508940]],
                names=['parent_oid'],
                dtype=[int])

    #function
    id_result = fetch_main_id(id_cat, id_creator=IdentifierCreator(name='main_id', colname='sdb_host_main_id'))
    oid_result = fetch_main_id(oid_cat, id_creator=OidCreator(name='parent_main_id', colname='parent_oid'))

    #assert
    assert id_result['main_id'][np.where(id_result['sdb_host_main_id'] == 'HD 166')] == 'HD    166'
    assert id_result['main_id'][np.where(id_result['sdb_host_main_id'] == '6 Cet')] == '*   6 Cet'
    assert oid_result['parent_main_id'][np.where(oid_result['parent_oid'] == 1800872)] == '* ksi UMa'
    assert oid_result['parent_main_id'][np.where(oid_result['parent_oid'] == 508940)] == 'MCC 541'


def test_lower_quality():
    assert lower_quality('A') == 'B'
    assert lower_quality('B') == 'C'
    assert lower_quality('C') == 'D'
    assert lower_quality('D') == 'E'
    assert lower_quality('E') == 'E'

