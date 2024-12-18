from provider.utils import *
from astropy.table import Table
import numpy as np

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

def test_query():
    link='http://dc.zah.uni-heidelberg.de/tap'
    adql_query="""
    SELECT TOP 10 main_id as Child_main_id, object_id as child_object_id
    FROM life_td.h_link
    JOIN life_td.ident as p on p.object_idref=parent_object_idref
    JOIN life_td.object on object_id=child_object_idref
    WHERE p.id = '* alf Cen'
    """
    table=query(link,adql_query)
    
    assert '* alf Cen A' in table['child_main_id']
    
    adql_query2="""
        SELECT TOP 10 *
        FROM TAP_UPLOAD.t1 
        WHERE child_main_id= '* alf Cen A'
        """
    upload_table=[table]
    table2=query(link,adql_query2,upload_table)
    
    assert table2['child_main_id'][0]=='* alf Cen A'
