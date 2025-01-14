from provider.utils import *
from astropy.table import Table, setdiff
import numpy as np

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

def test_fetch_main_id():
    #data
    exo_helptab = Table(data=[['*   3 Cnc b','6 Lyn b',    '*  14 Her b'],
                              ['*   3 Cnc', '*   6 Lyn','*  14 Her']],
                 names=['planet_main_id','host_main_id'], 
                 dtype=[object, object])

    result = Table(data=[['*   6 Lyn b','*  14 Her b'],
                         ['6 Lyn b',    '*  14 Her b'],
                              [ '*   6 Lyn','*  14 Her']],
                 names=['sim_planet_main_id','planet_main_id','host_main_id'], 
                 dtype=[object,object, object])

    #function
    exo_helptab2=fetch_main_id(exo_helptab['planet_main_id','host_main_id'],
                             #host_main_id just present to create table in contrast to column
                             IdentifierCreator(name='sim_planet_main_id',colname='planet_main_id'))

    #assert
    assert len(setdiff(exo_helptab2,result)) == 0


