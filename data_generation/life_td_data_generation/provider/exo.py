""" 
Generates the data for the database for the provider Exo-Mercat. 
"""

import numpy as np #arrays
from astropy import io
from astropy.table import Table, Column, MaskedColumn, join, setdiff
from datetime import datetime

#self created modules
from utils.io import save, Path, stringtoobject
from provider.utils import fetch_main_id, IdentifierCreator, fill_sources_table, create_sources_table, query, distance_cut, ids_from_ident, create_provider_table
from sdata import empty_dict

def create_exo_helpertable(exo):
    """
    Creates helper table.
    
    If the exomercat server is not online a temporary method to 
    ingest old exomercat data is used.
    
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Helper table and dictionary of database table names and tables.
    :rtype: astropy.table.table.Table, dict(str,astropy.table.table.Table)
    """
    #---------------define query----------------------------------------
    adql_query="""SELECT *
                  FROM exomercat.exomercat_legacy"""
    #---------------obtain data-----------------------------------------
    try:
        exo['provider'] = create_provider_table('Exo-MerCat',
                                        "http://archives.ia2.inaf.it/vo/tap/projects",
                                        '2020A&C....3100370A')
        exo_helptab=query(exo['provider']['provider_url'][0],adql_query)
    except:   
        exo['provider'] = create_provider_table('Exo-MerCat',
                                        "http://archives.ia2.inaf.it/vo/tap/projects",
                                        '2020A&C....3100370A','2023-02-05')        
        exo_helptab=io.ascii.read(
                Path().additional_data+"exo-mercat05-02-2023_v2.0.csv")
        exo_helptab=stringtoobject(exo_helptab,3000)
        exo['provider']['provider_access']=['2023-02-05']
        print('loading exomercat version from',exo['provider']['provider_access'])

    
        
    #----------------putting object main identifiers together-----------
    
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
    #join exo_helptab on host_main_id and sim_objects main_id
    # Unfortunately exomercat does not provide distance measurements so we
    # relie on matching it to simbad for enforcing the database cutoff of 20 pc.
    # The downside is, that the match does not work so well due to different
    # identifier notation which means we loose some objects. That is, however,
    # preferrable to having to do the work of checking the literature.
    # A compromise is to keep the list of objects I lost for later improvement.
    
    exo_helptab_before_distance_cut=exo_helptab
    exo_helptab=distance_cut(exo_helptab,'main_id')

    # removing whitespace in front of main_id and name.
    # done after distance_cut function to prevent missing values error
    for i in range(len(exo_helptab)):
        exo_helptab['planet_main_id'][i]=exo_helptab['planet_main_id'][i].strip()
        exo_helptab['main_id'][i]=exo_helptab['main_id'][i].strip()
        exo_helptab['name'][i]=exo_helptab['name'][i].strip()
        
    #fetching simbad main_id for planet since sometimes exomercat planet main id is not the same
    exo_helptab2=fetch_main_id(exo_helptab['planet_main_id','host_main_id'],
                             #host_main_id just present to create table in contrast to column
                             IdentifierCreator(name='sim_planet_main_id',colname='planet_main_id'))

    
    notinsimbad=exo_helptab['planet_main_id'][np.where(np.in1d(
            exo_helptab['planet_main_id'],exo_helptab2['planet_main_id'],
            invert=True))]
    #I use a left join as otherwise I would loose some objects that are not in simbad
    exo_helptab=join(exo_helptab,
                            exo_helptab2['sim_planet_main_id','planet_main_id'],
                            keys='planet_main_id',join_type='left')

    #show which elements from exo_helptab were not found in sim_objects
    exo_helptab_before_distance_cut['name']=exo_helptab_before_distance_cut['name'].astype(object)
    removed_objects=setdiff(exo_helptab_before_distance_cut,exo_helptab,keys=['name'])
    save([removed_objects],['exomercat_removed_objects'])
    return exo_helptab,exo

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
    exo_ident=Table(names=['main_id','id'],dtype=[object,object])
    exo_helptab['old_planet_main_id']=exo_helptab['planet_main_id']
    #why not found??
    for i in range(len(exo_helptab)):
        if exo_helptab['sim_planet_main_id'][i]!='':
            #not included as already in simbad provider
            #exo_ident.add_row([exo_helptab['sim_planet_main_id'][i],
                               #exo_helptab['sim_planet_main_id'][i]]) 
            if exo_helptab['planet_main_id'][i]!=exo_helptab['sim_planet_main_id'][i]:
                exo_ident.add_row([exo_helptab['sim_planet_main_id'][i],
                               exo_helptab['planet_main_id'][i]])
            if exo_helptab['planet_main_id'][i]!=exo_helptab['name'][i] and \
                    exo_helptab['sim_planet_main_id'][i]!=exo_helptab['name'][i]:
                exo_ident.add_row([exo_helptab['sim_planet_main_id'][i],
                               exo_helptab['name'][i]])
            exo_helptab['planet_main_id'][i]=exo_helptab['sim_planet_main_id'][i]
        else:
            exo_ident.add_row([exo_helptab['planet_main_id'][i],
                               exo_helptab['planet_main_id'][i]])
            if exo_helptab['planet_main_id'][i]!=exo_helptab['name'][i]:
                exo_ident.add_row([exo_helptab['planet_main_id'][i],
                               exo_helptab['name'][i]])
    exo_ident['id_ref']=[exo['provider']['provider_bibcode'][0] for j in range(len(exo_ident))]

    # TBD: I had a wrong double object though currently not any longer
    #print("""TBD: I have a wrong double object because of different amount of white
    #      spaces between catalog and number""")
    #print(exo_ident[np.where(exo_ident['main_id']=='Wolf  940 b')])
    return exo_ident

def create_objects_table(exo):
    """
    Creates objects table.
    
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Objects table.
    :rtype: astropy.table.table.Table
    """
    # tbd at one point: I think I want to add hosts to object
    exo_objects=Table(names=['main_id','ids'],dtype=[object,object])
    exo_objects=ids_from_ident(exo['ident']['main_id','id'],exo_objects)
    #grouped_exo_ident=exo_ident.group_by('main_id')
    #ind=grouped_exo_ident.groups.indices
    #for i in range(len(ind)-1):
    # -1 is needed because else ind[i+1] is out of bonds
     #   ids=[]
      #  for j in range(ind[i],ind[i+1]):
       #     ids.append(grouped_exo_ident['id'][j])
        #ids="|".join(ids)
        #exo_objects.add_row([grouped_exo_ident['main_id'][ind[i]],ids])
    exo_objects['type']=['pl' for j in range(len(exo_objects))]
    return exo_objects

def create_mes_mass_pl_table(exo_helptab,exo):
    """
    Creates planetary mass measurement table.
    
    :param exo_helptab: Exo-Mercat helper table.
    :type exo_helptab: astropy.table.table.Table
    :param exo: Dictionary of database table names and tables.
    :type exo: dict(str,astropy.table.table.Table)
    :returns: Planetary mass measurement table.
    :rtype: astropy.table.table.Table
    """
    #initialize columns exo_helptab['mass_pl_rel'] and exo_helptab['mass_pl_err']
    exo_helptab['mass_pl_err']=Column(dtype=float,length=len(exo_helptab))
    exo_helptab['mass_pl_rel']=Column(dtype=object,length=len(exo_helptab))
    exo_helptab['mass_pl_qual']=MaskedColumn(dtype=object,length=len(exo_helptab))
    exo_helptab['mass_pl_qual']=['?' for j in range(len(exo_helptab))]
    #transforming mass errors from upper (mass_max) and lower (mass_min) error
    # into instead error (mass_error) as well as relation (mass_pl_rel)
    for i in range(len(exo_helptab)):
        if type(exo_helptab['mass_max'][i])==np.ma.core.MaskedConstant or \
                  exo_helptab['mass_max'][i]==np.inf:
            if type(exo_helptab['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exo_helptab['mass_min'][i]==np.inf:
                exo_helptab['mass_pl_rel'][i]=None
                exo_helptab['mass_pl_err'][i]=1e+20
            else:
                exo_helptab['mass_pl_rel'][i]='<'
                exo_helptab['mass_pl_err'][i]=exo_helptab['mass_min'][i]
                exo_helptab['mass_pl_qual'][i]='C'
                # tbd: check if relation is correct in case of maximum error 
                # on a lower limit value')
        else:
            if type(exo_helptab['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exo_helptab['mass_min'][i]==np.inf:
                exo_helptab['mass_pl_rel'][i]='>'
                exo_helptab['mass_pl_err'][i]=exo_helptab['mass_max'][i]
                exo_helptab['mass_pl_qual'][i]='C'
            else:
                exo_helptab['mass_pl_rel'][i]='='
                exo_helptab['mass_pl_err'][i]=max(exo_helptab['mass_max'][i],
                                        exo_helptab['mass_min'][i])
                exo_helptab['mass_pl_qual'][i]='B'
    exo_mes_mass_pl=exo_helptab['planet_main_id','mass','mass_pl_err','mass_url',
                            'mass_pl_rel','mass_pl_qual']
    exo_mes_mass_pl.rename_columns(['planet_main_id','mass','mass_url'],
                                    ['main_id','mass_pl_value','mass_pl_ref'])
    #remove masked rows
    exo_mes_mass_pl.remove_rows(exo_mes_mass_pl['mass_pl_value'].mask.nonzero()[0])
    #remove null values
    exo_mes_mass_pl=exo_mes_mass_pl[np.where(exo_mes_mass_pl['mass_pl_value']!=1e+20)]
    #tbd: include masssini measurements from exomercat
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
    exo = empty_dict.copy()
    
    #wait, how do I do the distance cut if I don't have that in here?
    #there is a distance cut function using simbad objects
    exo_helptab,exo=create_exo_helpertable(exo)

    exo['ident']=create_ident_table(exo_helptab,exo)
    exo['objects']=create_objects_table(exo)
    exo['mes_mass_pl']=create_mes_mass_pl_table(exo_helptab,exo)
    exo['h_link']=create_h_link_table(exo_helptab,exo)
    exo['sources']=create_exo_sources_table(exo)

    save(list(exo.values()),['exo_'+ element for element in list(exo.keys())])
    return exo
