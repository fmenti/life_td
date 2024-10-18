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
from provider.utils import fetch_main_id, IdentifierCreator, sources_table, create_provider_table

def create_gk_helpertable(distance_cut_in_pc,gk_ref):
    plx_in_mas_cut=1000./distance_cut_in_pc

    #loading table obtained via direct communication from Grant Kennedy
    gk_disks=io.votable.parse_single_table(
        Path().additional_data+"sdb_30pc_09_02_2024.xml").to_table()
    #transforming from string type into object to have variable length
    gk_disks=stringtoobject(gk_disks,212)
    #removing objects with plx_value=='None' or masked entries
    if len(gk_disks['plx_value'].mask.nonzero()[0])>0:
        print('careful, masked entries in plx_value')

    #sorting out everything with plx_value too big
    gk_disks=gk_disks[np.where(gk_disks['plx_value']>plx_in_mas_cut)]
    #adds the column for object type
    gk_disks['type']=['di' for j in range(len(gk_disks))]
    gk_disks['disks_ref']=[gk_ref
                        for j in range(len(gk_disks))]
    #making sure identifiers are unique
    ind=gk_disks.group_by('id').groups.indices
    for i in range(len(ind)-1):
        l=ind[i+1]-ind[i]
        if l==2:
            gk_disks['id'][ind[i]]=gk_disks['id'][ind[i]]+'a'
            gk_disks['id'][ind[i]+1]=gk_disks['id'][ind[i]+1]+'b'
        if l>2:
            print('more than two disks with same name')
    #fetching updated main identifier of host star from simbad
    gk_disks.rename_column('main_id','gk_host_main_id')
    gk_disks=fetch_main_id(gk_disks,
                           IdentifierCreator(name='main_id',colname='gk_host_main_id'))
    return gk_disks

def provider_gk(distance_cut_in_pc):
    """
    Optains and arranges disk data.

    :returns: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data and basic disk data.
    :rtype: list(astropy.table.table.Table)
    """
    
    gk = empty_cat.copy()
    gk['provider'] = create_provider_table('Grant Kennedy Disks',
                                        'http://drgmk.com/sdb/',
                                        'priv. comm.','2024-02-09')
    
    gk_disks = create_gk_helpertable(distance_cut_in_pc,
                          gk['provider']['provider_bibcode'][0])

    #--------------creating h_link table ---------------------
    gk['h_link']=gk_disks['id','main_id','disks_ref']
    gk['h_link'].rename_columns(['main_id','disks_ref'],
                             ['parent_main_id','h_link_ref'])
    gk['h_link'].rename_column('id','main_id')
    #--------------creating objects table  --------------------
    gk_disks['ids']=gk_disks['id']#because only single id per source given
    gk['objects']=gk_disks['id','ids','type']
    gk['objects'].rename_column('id','main_id')
    #--------------creating ident table ----------------------
    gk['ident']=gk_disks['ids','id','disks_ref']
    # would prefer to use id instad of ids paremeter but this raises an error
    # so I use ids which has the same content as id
    gk['ident'].rename_columns(['ids','disks_ref'],['main_id','id_ref'])
    #--------------creating sources table --------------------
    
    tables=[gk['provider'],gk_disks]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['disks_ref']]
    for cat,ref in zip(tables,ref_columns):
        gk['sources']=sources_table(cat,ref,gk['provider']['provider_name'][0],
                                gk['sources'])
    #--------------creating disk_basic table ------------------
    gk['disk_basic']=gk_disks['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    
    for column in ['rdisk_bb','e_rdisk_bb']:
        if len(gk_disks[column].mask.nonzero()[0])>0:
            print('careful, masked entries in ', column)
    gk['disk_basic'].rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    gk['disk_basic']=gk['disk_basic'][np.where(np.isfinite(
                                                gk['disk_basic']['rad_value']))]

    save(list(gk.values()),['gk_'+ element for element in list(gk.keys())])
    return gk