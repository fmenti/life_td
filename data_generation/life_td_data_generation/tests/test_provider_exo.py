from provider.exo2 import *
from astropy.table import Table,join,setdiff
from sdata import empty_dict
import numpy as np

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
    exo_helptab = Table(data=[['*   3 Cnc b','6 Lyn b',    '6 Lyn b2'],
                              ['*   3 Cnc b', '*   6 Lyn b',''],
                              ['*   3 Cnc b','muster_exoname','muster_exoname2']],
                 names=['planet_main_id','sim_planet_main_id','exomercat_name'], 
                 dtype=[object, object,object])
    exo = empty_dict.copy()
    exo_ref='2020A&C....3100370A'
    exo['provider'] = Table(data = [[exo_ref]],
                            names = ['provider_bibcode'],
                            dtype=[object])
                        
    # function
    exo_ident, exo_helptab = create_ident_table(exo_helptab,exo)

    # assert
    # exo_ident
    assert exo_ident.colnames == ['main_id','id','id_ref']
    # planet_main_id didn't get an id column entry if sim_planet_main_id exists
    mask=np.isin(exo_ident['id'],['6 Lyn b']) 
    assert len(np.nonzero(mask)[0])==0
    mask2=np.isin(exo_ident['id'],['6 Lyn b2']) 
    assert len(np.nonzero(mask2)[0])==1
    # and did get if sim_planet_main_id doesn't
    assert  exo_ident['id_ref'][np.where(exo_ident['id'] == '6 Lyn b2')] == exo_ref
    # id column was created with ref exo when no simbad name there
    assert exo_ident['id_ref'][np.where(exo_ident['id'] == 'muster_exoname2')] == exo_ref
    # and with ref sim for sim_planet_name!=''
    assert  exo_ident['id_ref'][np.where(exo_ident['id'] == '*   6 Lyn b')] == '2000A&AS..143....9W'
    assert  exo_ident['id_ref'][np.where(exo_ident['id'] == 'muster_exoname')] == exo_ref
    # exo_helptab
    # updated exo_helptab where sim_main_id exists
    assert exo_helptab['planet_main_id'][np.where(exo_helptab['sim_planet_main_id']=='*   6 Lyn b')] == '*   6 Lyn b'
    # and don't updated exo_helptab where sim_main_id doesn't exist
    assert exo_helptab['planet_main_id'][np.where(exo_helptab['sim_planet_main_id']=='')] == '6 Lyn b2'

def test_exo_create_objects_table():
    # data
    exo = empty_dict.copy()
    exo['ident'] = Table(data = [['*   3 Cnc b','*   6 Lyn b','*   6 Lyn b'],
                                ['*   3 Cnc b','6 Lyn b','*   6 Lyn b']],
                            names = ['main_id','id'],
                            dtype=[object,object])
    
    # function
    exo_objects = create_objects_table(exo)
    # type got correctly assigned
    assert exo_objects['type'][np.where(exo_objects['main_id']=='*   3 Cnc b')] == 'pl'
    # unique main id
    assert len(exo_objects[np.where(exo_objects['main_id']=='*   6 Lyn b')]) == 1

def test_create_mes_mass_pl_table():
    #data
    m=MaskedColumn(data=[20.76,0,5.025,1e+20],name='mass',mask=[False,True,False,False])
    mmax=MaskedColumn(data=[0.73,0,0.873,1],name='mass_max',mask=[False,True,False,False])
    mmin=MaskedColumn(data=[0.73,0,1.067,1],name='mass_min',mask=[False,True,False,False])
    murl=MaskedColumn(data=['eu','','2022ApJS..262...21F','test'],name='mass_url',mask=[False,True,False,False])
    exo_helptab= Table(data = [['*   3 Cnc b','*   4 Mon B .01','*   6 Lyn b','testname'],
                               m,murl,mmax,mmin],
                            names = ['planet_main_id','mass','mass_url','mass_max','mass_min'],
                            dtype=[object,float,object,float,float])
    #function
    exo_mes_mass_pl=create_mes_mass_pl_table(exo_helptab)

    #assert
    # keep only non masked entries
    assert len(exo_mes_mass_pl) == 2
    assert exo_mes_mass_pl['mass_pl_value'][np.where(exo_mes_mass_pl['main_id']=='*   3 Cnc b')] == 20.76
    assert exo_mes_mass_pl['mass_pl_err'][np.where(exo_mes_mass_pl['main_id']=='*   6 Lyn b')] == 1.067

def test_create_mes_mass_pl_table_new():
    #data
    m=MaskedColumn(data=[20.76,0,5.025,1e+20],name='mass',mask=[False,True,False,False])
    mmax=MaskedColumn(data=[0.73,0,0.873,1],name='mass_max',mask=[False,True,False,False])
    mmin=MaskedColumn(data=[0.73,0,1.067,1],name='mass_min',mask=[False,True,False,False])
    murl=MaskedColumn(data=['eu','','2022ApJS..262...21F','test'],name='mass_url',mask=[False,True,False,False])
    exo_helptab= Table(data = [['*   3 Cnc b','*   4 Mon B .01','*   6 Lyn b','testname'],
                               m,murl,mmax,mmin],
                            names = ['planet_main_id','mass','mass_url','mass_max','mass_min'],
                            dtype=[object,float,object,float,float])
    #function
    exo_mes_mass_pl=create_mes_mass_pl_table(exo_helptab)

    #assert
    # keep only non masked entries
    assert len(exo_mes_mass_pl) == 2
    assert exo_mes_mass_pl['mass_pl_value'][np.where(exo_mes_mass_pl['main_id']=='*   3 Cnc b')] == 20.76
    assert exo_mes_mass_pl['mass_pl_err_max'][np.where(exo_mes_mass_pl['main_id']=='*   6 Lyn b')] == 0.873
    assert exo_mes_mass_pl['mass_pl_err_main'][np.where(exo_mes_mass_pl['main_id']=='*   6 Lyn b')] == 1.067
    
