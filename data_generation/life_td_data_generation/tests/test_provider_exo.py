from provider.exo2 import *
from astropy.table import Table,join,setdiff
from provider.utils import query
from sdata import empty_dict
import numpy as np


def test_exo_queryable():
    exo_provider=Table()
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    adql_query=["""SELECT top 10 *
                FROM exomercat.exomercat """]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
    print(exo.colnames)
    assert type(exo)==type(Table())

def test_exo_columns_queryable():  
    exo_provider=Table()           
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    adql_query=["""SELECT TOP 10
                    main_id, host, exomercat_name,
                    binary,letter
                    FROM exomercat.exomercat"""]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
    assert exo.colnames==['main_id', 'host', 'exomercat_name', 'binary', 'letter']

def test_exo_main_object_ids():
    a=MaskedColumn(data=['*   3 Cnc', '*   4 Mon', ''],name='main_id',mask=[False,False,True])
    b=MaskedColumn(data=['','B',''],name='binary',mask=[True,False,True])
    cat = Table(data=[a,
                       ['3 Cnc', 'TIC 72090501', '6 Lyn'],
                     b,
                     ['b','.01','b']],
                 names=['main_id','host','binary','letter'], 
                 dtype=[object, object,object,object])
    cat=create_object_main_id(cat)
    assert list(cat['host_main_id'])==['*   3 Cnc','*   4 Mon B','6 Lyn']
    assert list(cat['planet_main_id'])==['*   3 Cnc b','*   4 Mon B .01','6 Lyn b']

def test_exo_create_ident_table():
    # data
    exo_helptab = Table(data=[['*   3 Cnc b','*   4 Mon B .01','6 Lyn b'],
                               ['*   3 Cnc b', '', '*   6 Lyn b'],
                     ['*   3 Cnc b','*   4 Mon B .01','muster_exoname']],
                 names=['planet_main_id','sim_planet_main_id','exomercat_name'], 
                 dtype=[object, object,object])
    exo = empty_dict.copy()
    exo['provider'] = Table(data = [['2020A&C....3100370A']],
                            names = ['provider_bibcode'],
                            dtype=[object])
                        
    # function
    exo_ident, exo_helptab = create_ident_table(exo_helptab,exo)

    # assert

    assert exo_ident.colnames == ['main_id','id','id_ref']
    #this is a complicated way for assessing that planet_main_id column was created when sim_planet_name==''
    main_id_table=Table(data=[['*   4 Mon B .01'],['*   4 Mon B .01'],['2020A&C....3100370A']],
                        names=['main_id','id','id_ref'],
                        dtype=[object,object,object])
    assert  len(join(exo_ident,main_id_table,keys='id'))==1
    # and not for sim_planet_name!=''
    #### here not working yet### because exomercat_name same as sim_planet_main_id...
    main_id_table=Table(data=[['*   6 Lyn b'],['*   6 Lyn b'],['2020A&C....3100370A']],
                        names=['main_id','id','id_ref'],
                        dtype=[object,object,object])
    assert  len(join(exo_ident,main_id_table,keys='id'))==0
    # updated exo_helptab where sim_main_id exists
    assert exo_helptab['planet_main_id'][np.where(exo_helptab['sim_planet_main_id']=='*   6 Lyn b')] == '*   6 Lyn b'
    assert exo_helptab['planet_main_id'][np.where(exo_helptab['sim_planet_main_id']=='')] == '*   4 Mon B .01'

