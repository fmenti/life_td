from provider.utils import *
from astropy.table import Table
import numpy as np


def test_query():
    link='http://dc.zah.uni-heidelberg.de/tap'
    adql_query="""
    SELECT main_id as Child_main_id, object_id as child_object_id
    FROM life_td.h_link
    JOIN life_td.ident as p on p.object_idref=parent_object_idref
    JOIN life_td.object on object_id=child_object_idref
    WHERE p.id = '* alf Cen'
    """
    table=query(link,adql_query)
    
    assert '* alf Cen A' in table['child_main_id']
    
    adql_query2="""
        SELECT *
        FROM TAP_UPLOAD.t1 
        WHERE child_main_id= '* alf Cen A'
        """
    upload_table=[table]
    table2=query(link,adql_query2,upload_table)
    
    assert table2['child_main_id'][0]=='* alf Cen A'