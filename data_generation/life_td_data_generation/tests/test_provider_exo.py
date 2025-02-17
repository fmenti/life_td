from provider.exo import *
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

def test_assign_quality_elementwise():
    #data
    exo_helptab= Table(data = [[ 1e+20,1e+20,1],
                              [1e+20,1,1]],
                            names = ['mass_max','mass_min'],
                            dtype=[float,float])
    exo_helptab['mass_pl_qual']=MaskedColumn(dtype=object,length=len(exo_helptab))

    #function
    qual_b=assign_quality_elementwise(exo_helptab,'mass',2)
    qual_c=assign_quality_elementwise(exo_helptab,'mass',1)
    qual_d=assign_quality_elementwise(exo_helptab,'mass',0)

    #assert
    assert qual_b == 'B'
    assert qual_c == 'C'
    assert qual_d == 'D'

def test_assign_quality():
    # data
    exo_helptab= Table(data = [[ 1e+20,1e+20,1],
                              [1e+20,1,1],
                              [ 1e+20,1e+20,1],
                              [1e+20,1,1]],
                            names = ['mass_max','mass_min','msini_max','msini_min'],
                            dtype=[float,float,float,float])
    # function
    exo_helptab=assign_quality(exo_helptab,['mass','msini'])
    # assert
    assert exo_helptab['mass_pl_qual'][0] == 'D'
    assert exo_helptab['mass_pl_qual'][1] == 'C'
    assert exo_helptab['mass_pl_qual'][2] == 'B'
    assert exo_helptab['msini_pl_qual'][0] == 'D'
    assert exo_helptab['msini_pl_qual'][1] == 'C'
    assert exo_helptab['msini_pl_qual'][2] == 'B'

def test_deal_with_mass_nullvalues():
    # data
    m=MaskedColumn(data=[20.76,0.1,0,np.inf,1e+20],
                   name='mass',mask=[False,False,True,False,False])
    exo_helptab= Table(data = [['*   3 Cnc b','*   4 Mon B .01','testname1','*   6 Lyn b','testname2'],
                               m],
                            names = ['main_id','mass'],
                            dtype=[object,float])
    
    # function
    exo_helptab = deal_with_mass_nullvalues(exo_helptab,['mass'])
    # assert
    assert exo_helptab['mass'][np.where(exo_helptab['main_id']=='testname1')]==1e+20
    assert exo_helptab['mass'][np.where(exo_helptab['main_id']=='*   6 Lyn b')]==1e+20

def test_create_para_exo_mes_mass_pl():
    # data
    msini=MaskedColumn(data=[20.76,0.1,0,5.025,1e+20],
                   name='mass',mask=[False,False,True,False,False])
    msinimax=MaskedColumn(data=[np.inf,0,0,0.873,1],
                      name='mass_max',mask=[False,True,True,False,False])
    msinimin=MaskedColumn(data=[np.inf,1,0,1.067,1],
                      name='mass_min',mask=[False,False,False,False,False])
    msiniurl=MaskedColumn(data=['eu','test','','2022ApJS..262...21F','test'],
                      name='mass_url',mask=[False,False,True,False,False])
    mprov=MaskedColumn(data=['Msini','Mass','','Mass','Msini'],
                      name='bestmass_provenance',mask=[False,False,True,False,False])
    exo_helptab= Table(data = [['*   3 Cnc b','*   4 Mon B .01','testname1','*   6 Lyn b','testname2'],
                               msini,msinimax,msinimin,msiniurl,mprov],
                            names = ['planet_main_id','msini','msini_max','msini_min',
                                    'msini_url','bestmass_provenance'],
                            dtype=[object,float,float,float,object,object])
    exo_helptab['msini_pl_qual']=MaskedColumn(dtype=object,length=len(exo_helptab))
    exo_helptab['msini_pl_qual']=['?' for j in range(len(exo_helptab))]
    
    # function
    sinitable = create_para_exo_mes_mass_pl(exo_helptab,'msini','True')
    # assert
    assert sinitable['mass_pl_value'][np.where(sinitable['main_id']=='*   3 Cnc b')] == 20.76
    assert sinitable['mass_pl_err_max'][np.where(sinitable['main_id']=='*   6 Lyn b')] == 0.873
    assert sinitable['mass_pl_err_min'][np.where(sinitable['main_id']=='*   6 Lyn b')] == 1.067

def test_betterthan():
    assert betterthan('A','B') == True
    assert betterthan('A','?') == True
    assert betterthan('B','A') == False
    assert betterthan('A','A') == False
    assert betterthan('E','?') == True

def test_bestmass_better_qual():
    bestmass=['Mass','Mass','Mass','Msini','Msini','Msini']
    qual_msini=['B','A','C','B','A','C']
    qual_mass=['B','B','B','B','B','B']
    result_qual_msini=['C','C','C','B','A','C']
    result_qual_mass=['B','B','B','C','B','D']
    for i in range(len(bestmass)):
        msini,mass=bestmass_better_qual(bestmass[i],qual_msini[i],qual_mass[i])
        if i == 0:
            assert msini==result_qual_msini[i]
            assert mass==result_qual_mass[i]

def test_assign_new_qual():
    exo_mes_mass_pl = Table(data = [['*   3 Cnc b','*   3 Cnc b'],
                                   ['B','B'],
                                   ['True','False'],
                                   ['Mass','Mass']],
                            names = ['main_id','mass_pl_qual','mass_pl_sini_flag','bestmass_provenance'],
                            dtype=[object,object,object,object])
    exo_mes_mass_pl = assign_new_qual(exo_mes_mass_pl,'*   3 Cnc b','True','C')
    assert exo_mes_mass_pl['mass_pl_qual'][np.where(exo_mes_mass_pl['mass_pl_sini_flag']=='True')]=='C'

def test_align_quality_with_bestmass():
    exo_mes_mass_pl = Table(data = [['*   3 Cnc b','*   3 Cnc b'],
                                   ['B','B'],
                                   ['True','False'],
                                   ['Mass','Mass']],
                            names = ['main_id','mass_pl_qual','mass_pl_sini_flag','bestmass_provenance'],
                            dtype=[object,object,object,object])
    exo_mes_mass_pl=align_quality_with_bestmass(exo_mes_mass_pl)
    assert exo_mes_mass_pl['mass_pl_qual'][np.where(exo_mes_mass_pl['mass_pl_sini_flag']=='True')]=='C'

def test_create_mes_mass_pl_table():
    #data
    m=MaskedColumn(data=[20.76,0.1,0,5.025,1e+20],
                   name='mass',mask=[False,False,True,False,False])
    mmax=MaskedColumn(data=[np.inf,0,0,0.873,1],
                      name='mass_max',mask=[False,True,True,False,False])
    mmin=MaskedColumn(data=[np.inf,1,0,1.067,1],
                      name='mass_min',mask=[False,False,False,False,False])
    murl=MaskedColumn(data=['eu','test','','2022ApJS..262...21F','test'],
                      name='mass_url',mask=[False,False,True,False,False])
    msini=MaskedColumn(data=[20.76,0.1,0,5.025,1e+20],
                   name='mass',mask=[False,False,True,False,False])
    msinimax=MaskedColumn(data=[np.inf,0,0,0.873,1],
                      name='mass_max',mask=[False,True,True,False,False])
    msinimin=MaskedColumn(data=[np.inf,1,0,1.067,1],
                      name='mass_min',mask=[False,False,False,False,False])
    msiniurl=MaskedColumn(data=['eu','test','','2022ApJS..262...21F','test'],
                      name='mass_url',mask=[False,False,True,False,False])
    mprov=MaskedColumn(data=['Msini','Mass','','Mass','Msini'],
                      name='bestmass_provenance',mask=[False,False,True,False,False])
    exo_helptab= Table(data = [['*   3 Cnc b','*   4 Mon B .01','testname1','*   6 Lyn b','testname2'],
                               m,murl,mmax,mmin,msini,msiniurl,msinimax,msinimin,mprov],
                            names = ['planet_main_id','mass','mass_url','mass_max','mass_min',
                                    'msini','msini_url','msini_max','msini_min','bestmass_provenance'],
                            dtype=[object,float,object,float,float,float,object,float,float,object])
    #function
    exo_mes_mass_pl=create_mes_mass_pl_table(exo_helptab)

    #assert
    # keep only non masked entries
    assert len(exo_mes_mass_pl) == 6
    table=exo_mes_mass_pl[np.where(exo_mes_mass_pl['mass_pl_sini_flag']=='False')]
    assert table['mass_pl_value'][np.where(table['main_id']=='*   3 Cnc b')] == 20.76
    assert table['mass_pl_err_max'][np.where(table['main_id']=='*   6 Lyn b')] == 0.873
    assert table['mass_pl_err_min'][np.where(table['main_id']=='*   6 Lyn b')] == 1.067
    assert table['mass_pl_err_max'][np.where(table['main_id']=='*   4 Mon B .01')] == 1e+20
    assert table['bestmass_provenance'][np.where(table['main_id']=='*   4 Mon B .01')] == 'Mass'
    

    
    
