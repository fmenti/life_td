from provider.utils import *


def test_create_provider_table_date_given():
    gk_provider = create_provider_table('Grant Kennedy Disks',
                                        'http://drgmk.com/sdb/',
                                        'priv. comm.','2024-02-09')
    assert gk_provider['provider_name'] == 'Grant Kennedy Disks'
    assert gk_provider['provider_url'] == 'http://drgmk.com/sdb/'
    assert gk_provider['provider_bibcode'] == 'priv. comm.'
    assert gk_provider['provider_access'] == '2024-02-09'
    
def test_create_provider_table_no_date_given():
    gk_provider = create_provider_table('Grant Kennedy Disks',
                                        'http://drgmk.com/sdb/',
                                        'priv. comm.')
    assert gk_provider['provider_access'] == datetime.now().strftime('%Y-%m-%d')    

