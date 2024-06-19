""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from astropy import io
from astropy.table import Table, Column, MaskedColumn, join, setdiff
from datetime import datetime

#self created modules
from utils.io import save, Path
from provider.utils import fetch_main_id, IdentifierCreator, sources_table, query, distance_cut, ids_from_ident
import sdata as sdc


def provider_exo(table_names,exo_list_of_tables,temp=False):
    """
    This function obtains and arranges exomercat data.
    
    If the exomercat server is not online a temporary method to 
    ingest old exomercat data was implemented and can be accessed by 
    setting temp=True as argument.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param exo_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :type exo_list_of_tables: list(astropy.table.table.Table)
    :param bool temp: Defaults to True.  Determines if the exomercat
        data gets queried (False) or loaded from an old version (True).
    :returns: List of astropy table containing
        reference data, object data, identifier data, object to object
        relation data, basic planetary data and planetary mass measurement 
        data.
    :rtype: list(astropy.table.table.Table)
    """
    
    #---------------define provider-------------------------------------
    exo_provider=Table()
    exo_provider['provider_name']=['Exo-MerCat']
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    exo_provider['provider_bibcode']=['2020A&C....3100370A']
    
    
    print('Creating ',exo_provider['provider_name'][0],' tables ...')
    #---------------define query----------------------------------------
    adql_query="""SELECT *
                  FROM exomercat.exomercat"""
    #---------------obtain data-----------------------------------------
    if temp:
        exomercat=io.ascii.read(
                Path().additional_data+"exo-mercat05-02-2023_v2.0.csv")
        exomercat=stringtoobject(exomercat,3000)
        exo_provider['provider_access']=['2023-02-05']

    else:
        exomercat=query(exo_provider['provider_url'][0],adql_query)
        exo_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    #----------------putting object main identifiers together-----------
    
    # initializing column
    exomercat['planet_main_id']=Column(dtype=object,
                                                length=len(exomercat))
    exomercat['host_main_id']=exomercat['main_id']

    for i in range(len(exomercat)):
        # if there is a main id use that
        if type(exomercat['main_id'][i])!=np.ma.core.MaskedConstant:
            hostname=exomercat['main_id'][i]
        # else use host column entry
        else:
            hostname=exomercat['host'][i]
        # if there is a binary entry
        if type(exomercat['binary'][i])!=np.ma.core.MaskedConstant:
            exomercat['host_main_id'][i]=hostname+' '+exomercat['binary'][i]
        else:
            exomercat['host_main_id'][i]=hostname
        exomercat['planet_main_id'][i]=exomercat[
                        'host_main_id'][i]+' '+exomercat['letter'][i]
    #join exomercat on host_main_id and sim_objects main_id
    # Unfortunately exomercat does not provide distance measurements so we
    # relie on matching it to simbad for enforcing the database cutoff of 20 pc.
    # The downside is, that the match does not work so well due to different
    # identifier notation which means we loose some objects. That is, however,
    # preferrable to having to do the work of checking the literature.
    # A compromise is to keep the list of objects I lost for later improvement.
    
    exo=exomercat
    exomercat=distance_cut(exomercat,'main_id')

    # removing whitespace in front of main_id and name.
    # done after distance_cut function to prevent missing values error
    for i in range(len(exomercat)):
        exomercat['planet_main_id'][i]=exomercat['planet_main_id'][i].strip()
        exomercat['main_id'][i]=exomercat['main_id'][i].strip()
        exomercat['name'][i]=exomercat['name'][i].strip()
        
    #fetching simbad main_id for planet since sometimes exomercat planet main id is not the same
    exomercat2=fetch_main_id(exomercat['planet_main_id','host_main_id'],
                             #host_main_id just present to create table in contrast to column
                             IdentifierCreator(name='sim_planet_main_id',colname='planet_main_id'))

    
    notinsimbad=exomercat['planet_main_id'][np.where(np.in1d(
            exomercat['planet_main_id'],exomercat2['planet_main_id'],
            invert=True))]
    #I use a left join as otherwise I would loose some objects that are not in simbad
    exomercat=join(exomercat,
                            exomercat2['sim_planet_main_id','planet_main_id'],
                            keys='planet_main_id',join_type='left')

    #show which elements from exomercat were not found in sim_objects
    exo['name']=exo['name'].astype(object)
    removed_objects=setdiff(exo,exomercat,keys=['name'])
    save([removed_objects],['exomercat_removed_objects'])

    #-------------exo_ident---------------
    exo_ident=Table(names=['main_id','id'],dtype=[object,object])
    exomercat['old_planet_main_id']=exomercat['planet_main_id']
    #why not found??
    for i in range(len(exomercat)):
        if exomercat['sim_planet_main_id'][i]!='':
            #not included as already in simbad provider
            #exo_ident.add_row([exomercat['sim_planet_main_id'][i],
                               #exomercat['sim_planet_main_id'][i]]) 
            if exomercat['planet_main_id'][i]!=exomercat['sim_planet_main_id'][i]:
                exo_ident.add_row([exomercat['sim_planet_main_id'][i],
                               exomercat['planet_main_id'][i]])
            if exomercat['planet_main_id'][i]!=exomercat['name'][i] and \
                    exomercat['sim_planet_main_id'][i]!=exomercat['name'][i]:
                exo_ident.add_row([exomercat['sim_planet_main_id'][i],
                               exomercat['name'][i]])
            exomercat['planet_main_id'][i]=exomercat['sim_planet_main_id'][i]
        else:
            exo_ident.add_row([exomercat['planet_main_id'][i],
                               exomercat['planet_main_id'][i]])
            if exomercat['planet_main_id'][i]!=exomercat['name'][i]:
                exo_ident.add_row([exomercat['planet_main_id'][i],
                               exomercat['name'][i]])
    exo_ident['id_ref']=[exo_provider['provider_bibcode'][0] for j in range(len(exo_ident))]

    # TBD: I had a wrong double object though currently not any longer
    #print("""TBD: I have a wrong double object because of different amount of white
    #      spaces between catalog and number""")
    #print(exo_ident[np.where(exo_ident['main_id']=='Wolf  940 b')])
    
    #-------------exo_objects-------------------------------------------
    # tbd at one point: I think I want to add hosts to object
    exo_objects=Table(names=['main_id','ids'],dtype=[object,object])
    exo_objects=ids_from_ident(exo_ident['main_id','id'],exo_objects)
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
    
    #-------------------exo_mes_mass_pl---------------------
    #initialize columns exomercat['mass_pl_rel'] and exomercat['mass_pl_err']
    exomercat['mass_pl_err']=Column(dtype=float,length=len(exomercat))
    exomercat['mass_pl_rel']=Column(dtype=object,length=len(exomercat))
    exomercat['mass_pl_qual']=MaskedColumn(dtype=object,length=len(exomercat))
    exomercat['mass_pl_qual']=['?' for j in range(len(exomercat))]
    #transforming mass errors from upper (mass_max) and lower (mass_min) error
    # into instead error (mass_error) as well as relation (mass_pl_rel)
    for i in range(len(exomercat)):
        if type(exomercat['mass_max'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_max'][i]==np.inf:
            if type(exomercat['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_min'][i]==np.inf:
                exomercat['mass_pl_rel'][i]=None
                exomercat['mass_pl_err'][i]=1e+20
            else:
                exomercat['mass_pl_rel'][i]='<'
                exomercat['mass_pl_err'][i]=exomercat['mass_min'][i]
                exomercat['mass_pl_qual'][i]='C'
                # tbd: check if relation is correct in case of maximum error 
                # on a lower limit value')
        else:
            if type(exomercat['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_min'][i]==np.inf:
                exomercat['mass_pl_rel'][i]='>'
                exomercat['mass_pl_err'][i]=exomercat['mass_max'][i]
                exomercat['mass_pl_qual'][i]='C'
            else:
                exomercat['mass_pl_rel'][i]='='
                exomercat['mass_pl_err'][i]=max(exomercat['mass_max'][i],
                                        exomercat['mass_min'][i])
                exomercat['mass_pl_qual'][i]='B'
    exo_mes_mass_pl=exomercat['planet_main_id','mass','mass_pl_err','mass_url',
                            'mass_pl_rel','mass_pl_qual']
    exo_mes_mass_pl.rename_columns(['planet_main_id','mass','mass_url'],
                                    ['main_id','mass_pl_value','mass_pl_ref'])
    #remove masked rows
    exo_mes_mass_pl.remove_rows(exo_mes_mass_pl['mass_pl_value'].mask.nonzero()[0])
    #remove null values
    exo_mes_mass_pl=exo_mes_mass_pl[np.where(exo_mes_mass_pl['mass_pl_value']!=1e+20)]
    #tbd: include masssini measurements from exomercat

    #-------------exo_h_link---------------
    exo_h_link=exomercat['planet_main_id', 'host_main_id']
    exo_h_link.rename_columns(['planet_main_id','host_main_id'],
                              ['main_id','parent_main_id'])
    exo_h_link['h_link_ref']=[exo_provider['provider_bibcode'][0] for j in range(len(exo_h_link))]
    #-------------exo_sources---------------
    ref_columns=[['provider_bibcode'],['h_link_ref'],['id_ref'],['mass_pl_ref']]
    exo_sources=Table()
    tables=[exo_provider,exo_h_link,exo_ident,exo_mes_mass_pl]
    for cat,ref in zip(tables,ref_columns):
        exo_sources=sources_table(cat,ref,exo_provider['provider_name'][0],exo_sources)
  
    for i in range(len(table_names)):
        if table_names[i]=='sources': exo_list_of_tables[i]=exo_sources
        if table_names[i]=='provider': exo_list_of_tables[i]=exo_provider 
        if table_names[i]=='objects': exo_list_of_tables[i]=exo_objects
        if table_names[i]=='ident': exo_list_of_tables[i]=exo_ident
        if table_names[i]=='h_link': exo_list_of_tables[i]=exo_h_link
        if table_names[i]=='mes_mass_pl': exo_list_of_tables[i]=exo_mes_mass_pl
        save([exo_list_of_tables[i][:]],['exo_'+table_names[i]])
    return exo_list_of_tables