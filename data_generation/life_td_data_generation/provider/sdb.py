""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from astropy import io
from astropy.table import Table
from datetime import datetime

#self created modules
from utils.io import stringtoobject, save, Path
from sdata import empty_cat
from provider.utils import fetch_main_id, IdentifierCreator, fill_sources_table, create_provider_table

def create_sdb_helpertable(distance_cut_in_pc,sdb_ref):
    plx_in_mas_cut=1000./distance_cut_in_pc

    #loading table obtained via direct communication from Grant Kennedy
    sdb_helptab=io.votable.parse_single_table(
        Path().additional_data+"sdb_30pc_09_02_2024.xml").to_table()
    #transforming from string type into object to have variable length
    sdb_helptab=stringtoobject(sdb_helptab,212)
    #removing objects with plx_value=='None' or masked entries
    if len(sdb_helptab['plx_value'].mask.nonzero()[0])>0:
        print('careful, masked entries in plx_value')

    #sorting out everything with plx_value too big
    sdb_helptab=sdb_helptab[np.where(sdb_helptab['plx_value']>plx_in_mas_cut)]
    #adds the column for object type
    sdb_helptab['type']=['di' for j in range(len(sdb_helptab))]
    sdb_helptab['disks_ref']=[sdb_ref
                        for j in range(len(sdb_helptab))]
    #making sure identifiers are unique
    ind=sdb_helptab.group_by('id').groups.indices
    for i in range(len(ind)-1):
        l=ind[i+1]-ind[i]
        if l==2:
            sdb_helptab['id'][ind[i]]=sdb_helptab['id'][ind[i]]+'a'
            sdb_helptab['id'][ind[i]+1]=sdb_helptab['id'][ind[i]+1]+'b'
        if l>2:
            print('more than two disks with same name')
    #fetching updated main identifier of host star from simbad
    sdb_helptab.rename_column('main_id','sdb_host_main_id')
    sdb_helptab=fetch_main_id(sdb_helptab,
                           IdentifierCreator(name='main_id',colname='sdb_host_main_id'))
    return sdb_helptab

def create_h_link_table(sdb_helptab):
    h_link=sdb_helptab['id','main_id','disks_ref']
    h_link.rename_columns(['main_id','disks_ref'],
                             ['parent_main_id','h_link_ref'])
    h_link.rename_column('id','main_id')
    return h_link

def create_objects_table(sdb_helptab):        
    objects_table=sdb_helptab['id','ids','type']
    objects_table.rename_column('id','main_id')
    return objects_table

def create_ident_table(sdb_helptab):
    ident_table=sdb_helptab['ids','id','disks_ref']
    # would prefer to use id instad of ids paremeter but this raises an error
    # so I use ids which has the same content as id
    ident_table.rename_columns(['ids','disks_ref'],['main_id','id_ref'])
    return ident_table

def create_sources_table(sdb_helptab,provider_table):
    tables=[provider_table,sdb_helptab]
    sources_table=Table()
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['disks_ref']]
    for cat,ref in zip(tables,ref_columns):
        sources_table=fill_sources_table(cat,ref,provider_table['provider_name'][0],
                                sources_table)
    return sources_table

def provider_sdb(distance_cut_in_pc):
    """
    Optains and arranges disk data.

    :returns: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data and basic disk data.
    :rtype: list(astropy.table.table.Table)
    """
    
    sdb = empty_cat.copy()
    sdb['provider'] = create_provider_table('Grant Kennedy Disks',
                                        'http://drgmk.com/sdb/',
                                        'priv. comm.','2024-02-09')
    
    sdb_helptab = create_sdb_helpertable(distance_cut_in_pc,
                          sdb['provider']['provider_bibcode'][0])

    sdb['h_link']=create_h_link_table(sdb_helptab)
    sdb_helptab['ids']=sdb_helptab['id']#because only single id per source given
    sdb['objects']=create_objects_table(sdb_helptab)
    sdb['ident']=create_ident_table(sdb_helptab)
    sdb['sources']=create_sources_table(sdb_helptab,sdb['provider'])
    #--------------creating disk_basic table ------------------
    sdb['disk_basic']=sdb_helptab['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    
    for column in ['rdisk_bb','e_rdisk_bb']:
        if len(sdb_helptab[column].mask.nonzero()[0])>0:
            print('careful, masked entries in ', column)
    sdb['disk_basic'].rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    sdb['disk_basic']=sdb['disk_basic'][np.where(np.isfinite(
                                                sdb['disk_basic']['rad_value']))]

    save(list(sdb.values()),['sdb_'+ element for element in list(sdb.keys())])
    return sdb