""" 
Generates the data for the database for the provider Exo-Mercat. 
Distance cut is performed by joining on life_db simbad objects.
"""

import numpy as np #arrays
from astropy import io
from astropy.table import Table, Column, MaskedColumn, join, setdiff, unique, vstack
from datetime import datetime

#self created modules
from utils.io import save, load, Path, stringtoobject
from provider.utils import fetch_main_id, IdentifierCreator, fill_sources_table, create_sources_table, query, distance_cut, ids_from_ident, create_provider_table, nullvalues, replace_value, lower_quality
from sdata import empty_dict

def query_or_load_exomercat():
    """
    Queries or loads exomercat table.
    
    If the exomercat server is not online a temporary method to 
    ingest old exomercat data is used.
    
    :returns: Dictionary of database table names and tables.
    :rtype: dict(str,astropy.table.table.Table)
    """
    exo = empty_dict.copy()
    adql_query="""SELECT *
                  FROM exomercat.exomercat"""
    
    try:
        exo['provider'] = create_provider_table('Exo-MerCat',
                                        "http://archives.ia2.inaf.it/vo/tap/projects",
                                        '2020A&C....3100370A')#change bibcode once paper on ads
        exo_helptab=query(exo['provider']['provider_url'][0],adql_query)
    except:   
        exo['provider'] = create_provider_table('Exo-MerCat',
                                        "http://archives.ia2.inaf.it/vo/tap/projects",
                                        '2020A&C....3100370A','2024-12-13')        
        exo_helptab=io.ascii.read(
                Path().additional_data+"exo-mercat13-12-2024_v2.0.csv")
        exo_helptab=stringtoobject(exo_helptab,3000)
        exo['provider']['provider_access']=['2024-12-13']
        print('loading exomercat version from',exo['provider']['provider_access'])
    return exo, exo_helptab

def create_object_main_id(exo_helptab):
    # initializing column
    exo_helptab['planet_main_id']=Column(dtype=object,
                                            length=len(exo_helptab))
    exo_helptab['host_main_id']=exo_helptab['main_id']

    for i in range(len(exo_helptab)):
        # if there is a main id use that
        if type(exo_helptab['main_id'][i])!=np.ma.core.MaskedConstant:
            hostname=exo_helptab['main_id'][i]
        # else use host column entry
        else:
            hostname=exo_helptab['host'][i]
        # if there is a binary entry
        if type(exo_helptab['binary'][i])!=np.ma.core.MaskedConstant:
            exo_helptab['host_main_id'][i]=hostname+' '+exo_helptab['binary'][i]
        else:
            exo_helptab['host_main_id'][i]=hostname
        exo_helptab['planet_main_id'][i]=exo_helptab[
                    'host_main_id'][i]+' '+exo_helptab['letter'][i]
    return exo_helptab

def create_exo_helpertable():
    """
    Creates helper table.
    
    :returns: Helper table and dictionary of database table names and tables.
    :rtype: astropy.table.table.Table, dict(str,astropy.table.table.Table)
    """

    exo, exo_helptab = query_or_load_exomercat() 
        
    exo_helptab=create_object_main_id(exo_helptab)
    #join exo_helptab on host_main_id and sim_objects main_id
    # Unfortunately exomercat does not provide distance measurements so we
    #   relie on matching it to simbad for enforcing the database cutoff of 20 pc.
    #   The downside is, that the match does not work so well due to different
    #   identifier notation which means we loose some objects. That is, however,
    #   preferrable to having to do the work of checking the literature.
    #   A compromise is to keep the list of objects I lost for later improvement.
    exo_helptab_before_distance_cut=exo_helptab
    exo_helptab=distance_cut(exo_helptab,'main_id')

    # removing whitespace in front of main_id and name.
    # done after distance_cut function to prevent missing values error
    for i in range(len(exo_helptab)):
        exo_helptab['planet_main_id'][i]=exo_helptab['planet_main_id'][i].strip()
        exo_helptab['main_id'][i]=exo_helptab['main_id'][i].strip()
        exo_helptab['exomercat_name'][i]=exo_helptab['exomercat_name'][i].strip()
        
    #fetching simbad main_id for planet since sometimes exomercat planet main id is not the same
    exo_helptab2=fetch_main_id(exo_helptab['planet_main_id','host_main_id'],
                             #host_main_id just present to create table in contrast to column
                             IdentifierCreator(name='sim_planet_main_id',colname='planet_main_id'))

    
    notinsimbad=exo_helptab['planet_main_id'][np.where(np.isin(
            exo_helptab['planet_main_id'],exo_helptab2['planet_main_id'],
            invert=True))]
    #I use a left join as otherwise I would loose some objects that are not in simbad
    exo_helptab=join(exo_helptab,
                            exo_helptab2['sim_planet_main_id','planet_main_id'],
                            keys='planet_main_id',join_type='left')

    #show which elements from exo_helptab were not found in sim_objects
    exo_helptab_before_distance_cut['exomercat_name']=exo_helptab_before_distance_cut['exomercat_name'].astype(object)
    removed_objects=setdiff(exo_helptab_before_distance_cut,exo_helptab,keys=['exomercat_name'])
    save([removed_objects],['exomercat_removed_objects'])
    return exo, exo_helptab

def create_ident_table(exo_helptab,exo):
    """
    Creates identifier table.
    
    :param exo_helptab: Exo-Mercat helper table.
    :type exo_helptab: astropy.table.table.Table
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Identifier table.
    :rtype: astropy.table.table.Table
    """
    exo_ident=Table(names=['main_id','id','id_ref'],dtype=[object,object,object])
    
    for i in range(len(exo_helptab)):
        if exo_helptab['sim_planet_main_id'][i]!='':
            main_id=exo_helptab['sim_planet_main_id'][i]
            exo_helptab['planet_main_id'][i]=exo_helptab['sim_planet_main_id'][i] 
            # issue: changing exo_helptab but not returning it
            [simbad_ref]=load(['sim_provider'])
            idref=simbad_ref['provider_bibcode'][0]
        else: 
            main_id = exo_helptab['planet_main_id'][i]
            # in else because otherwise will result in double sim main id 
            # from sim provider in case of sim being main id
            idref=exo['provider']['provider_bibcode'][0]
        exo_ident.add_row([main_id,main_id,idref])
        exo_ident.add_row([main_id,exo_helptab['exomercat_name'][i],
                           exo['provider']['provider_bibcode'][0]])
        
    exo_ident=unique(exo_ident)
    
    return exo_ident, exo_helptab

def create_objects_table(exo):
    """
    Creates objects table.
    
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Objects table.
    :rtype: astropy.table.table.Table
    """
    # tbd at one point: add hosts to object because currently parent not in sim is not in objects table
    
    exo_objects=Table(names=['main_id','ids'],dtype=[object,object])
    exo_objects=ids_from_ident(exo['ident']['main_id','id'],exo_objects)
    exo_objects['type']=['pl' for j in range(len(exo_objects))]
    return exo_objects

def assign_quality_elementwise(exo_helptab,para,i):
    qual='B'
    if exo_helptab[para+'_max'][i]==1e+20:
        qual=lower_quality(qual)
    if exo_helptab[para+'_min'][i]==1e+20:
        qual=lower_quality(qual)
    return qual

def assign_quality(exo_helptab,parameter):
    for para in parameter:
        exo_helptab[para+'_pl_qual']=MaskedColumn(dtype=object,length=len(exo_helptab))
        exo_helptab[para+'_pl_qual']=['?' for j in range(len(exo_helptab))]
        for i in range(len(exo_helptab)):
            exo_helptab[para+'_pl_qual'][i] = assign_quality_elementwise(exo_helptab,para,i)
    return exo_helptab

def deal_with_mass_nullvalues(exo_helptab,cols):
    for colname in cols:
        exo_helptab=nullvalues(exo_helptab,colname,1e+20,verbose=False)
        exo_helptab=replace_value(exo_helptab,colname,np.inf,1e+20)
    return exo_helptab

def create_para_exo_mes_mass_pl(exo_helptab,para,sini_flag):
    table=exo_helptab['planet_main_id',para,para+'_max',para+'_min',para+'_url',
                            para+'_pl_qual','bestmass_provenance']
    table.rename_columns(['planet_main_id',para,para+'_url',
                              para+'_max',para+'_min'],
                        ['main_id','mass_pl_value','mass_pl_ref',
                                'mass_pl_err_max','mass_pl_err_min'])
    if para=='msini':
        table.rename_column(para+'_pl_qual','mass_pl_qual')
    table['mass_pl_sini_flag']=[sini_flag for j in range(len(table))]
    #remove null values
    table=table[np.where(table['mass_pl_value']!=1e+20)]
    return table

def betterthan(qual1,qual2):
    #is it usefull to do all this work with the quality when I could just assign b for bestmass and c for otherone?
    #nearly done, finish this then do other one as well. call one exomercatqual other dbqual
    quals=['A','B','C','D','E','?']
    for i in [0,1,2,3,4]:
        if qual1 == quals[i] and qual2 in quals[i+1:]:
            result=True
            return result
        else:
            result=False 
    return result

def bestmass_better_qual(bestmass,qual_msini,qual_mass):
    if bestmass=='Mass':
        if betterthan(qual_msini,qual_mass) or qual_msini==qual_mass:
            qual_msini=lower_quality(qual_msini)
            qual_msini, qual_mass=bestmass_better_qual(bestmass,qual_msini,qual_mass)
    elif bestmass == 'Msini':
        if betterthan(qual_mass,qual_msini) or qual_msini==qual_mass:
            qual_mass=lower_quality(qual_mass)
            qual_msini, qual_mass=bestmass_better_qual(bestmass,qual_msini,qual_mass)
    return qual_msini, qual_mass

def assign_new_qual(exo_mes_mass_pl,main_id,flag,new_qual):
    temp=[np.where(exo_mes_mass_pl['main_id']==main_id)]
    for i in temp[0][0]:
        if exo_mes_mass_pl['mass_pl_sini_flag'][i]==flag:
            exo_mes_mass_pl['mass_pl_qual'][i]=new_qual
    return exo_mes_mass_pl

def align_quality_with_bestmass(exo_mes_mass_pl):
    grouped_mass_pl = exo_mes_mass_pl.group_by('main_id')
    ind=grouped_mass_pl.groups.indices
    for i in range(len(ind)-1):#range(len(ind)-1):
        if ind[i+1]-ind[i]==2:
            bestmass=grouped_mass_pl['bestmass_provenance'][ind[i]]
            if grouped_mass_pl['mass_pl_sini_flag'][ind[i]]=='True':
                qual_msini=grouped_mass_pl['mass_pl_qual'][ind[i]]
                qual_mass=grouped_mass_pl['mass_pl_qual'][ind[i]+1]
            else:
                qual_msini=grouped_mass_pl['mass_pl_qual'][ind[i]+1]
                qual_mass=grouped_mass_pl['mass_pl_qual'][ind[i]]
            new_qual_msini, new_qual_mass=bestmass_better_qual(bestmass,qual_msini,qual_mass)
            #assign them the right place
            #I think I need a where here
            if new_qual_msini!=qual_msini:
                main_id=grouped_mass_pl['main_id'][ind[i]]
                exo_mes_mass_pl=assign_new_qual(exo_mes_mass_pl,main_id,'True',new_qual_msini)
            if new_qual_mass!=qual_mass:
                main_id=grouped_mass_pl['main_id'][ind[i]]
                exo_mes_mass_pl=assign_new_qual(exo_mes_mass_pl,main_id,'False',new_qual_mass)
                
    return exo_mes_mass_pl
    

def create_mes_mass_pl_table(exo_helptab):
    """
    Creates planetary mass measurement table.
    
    :param exo_helptab: Exo-Mercat helper table.
    :type exo_helptab: astropy.table.table.Table
    :returns: Planetary mass measurement table.
    :rtype: astropy.table.table.Table
    """    

    # include other parameters r, a, e, i, p, status
    cols = ['mass_max','mass_min', 'mass','msini', 'msini_max', 'msini_min']
    exo_helptab = deal_with_mass_nullvalues(exo_helptab,cols)
        
    exo_helptab=assign_quality(exo_helptab,['mass','msini'])
        

    exo_mes_mass_pl1=create_para_exo_mes_mass_pl(exo_helptab,'mass','False')
    exo_mes_mass_pl2=create_para_exo_mes_mass_pl(exo_helptab,'msini','True')
        
    exo_mes_mass_pl=vstack([exo_mes_mass_pl1,exo_mes_mass_pl2])
    exo_mes_mass_pl=unique(exo_mes_mass_pl)
    exo_mes_mass_pl=align_quality_with_bestmass(exo_mes_mass_pl)  
    return exo_mes_mass_pl

def create_h_link_table(exo_helptab,exo):
    """
    Creates hierarchical link table.
    
    :param exo_helptab: Exo-Mercat helper table.
    :type exo_helptab: astropy.table.table.Table
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Hierarchical link table.
    :rtype: astropy.table.table.Table
    """
    exo_h_link=exo_helptab['planet_main_id', 'host_main_id']
    exo_h_link.rename_columns(['planet_main_id','host_main_id'],
                              ['main_id','parent_main_id'])
    exo_h_link['h_link_ref']=[exo['provider']['provider_bibcode'][0] for j in range(len(exo_h_link))]
    return exo_h_link

def create_exo_sources_table(exo):
    """
    Creates sources table.
    
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Sources table.
    :rtype: astropy.table.table.Table
    """
    tables=[exo['provider'],exo['h_link'],exo['ident'],exo['mes_mass_pl']]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['h_link_ref'],['id_ref'],['mass_pl_ref']]
    exo_sources=create_sources_table(tables,ref_columns,
                              exo['provider']['provider_name'][0])
    return exo_sources  

def provider_exo():
    """
    This function obtains and arranges exomercat data.
    
    :returns: Dictionary with names and astropy tables containing
        reference data, object data, identifier data, object to object
        relation data, basic planetary data and planetary mass measurement 
        data.
    :rtype: dict(str,astropy.table.table.Table)
    """
    
    exo,exo_helptab=create_exo_helpertable()

    exo['ident'],exo_helptab=create_ident_table(exo_helptab,exo)
    exo['objects']=create_objects_table(exo)
    exo['mes_mass_pl']=create_mes_mass_pl_table(exo_helptab)
    exo['h_link']=create_h_link_table(exo_helptab,exo)
    exo['sources']=create_exo_sources_table(exo)

    save(list(exo.values()),['exo_'+ element for element in list(exo.keys())])
    return exo
