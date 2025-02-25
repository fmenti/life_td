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
    cat = Table(data=[['HD 166', '6 Cet']],
                names=['sdb_host_main_id'],
                dtype=[object])

    #function
    result = fetch_main_id(cat, id_creator=IdentifierCreator(name='main_id', colname='sdb_host_main_id'))

    #assert
    assert result['main_id'][np.where(result['sdb_host_main_id'] == 'HD 166')] == 'HD    166'
    assert result['main_id'][np.where(result['sdb_host_main_id'] == '6 Cet')] == '*   6 Cet'


def test_lower_quality():
    assert lower_quality('A') == 'B'
    assert lower_quality('B') == 'C'
    assert lower_quality('C') == 'D'
    assert lower_quality('D') == 'E'
    assert lower_quality('E') == 'E'

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
