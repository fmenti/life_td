""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from astropy import io
from astropy.table import Table
from datetime import datetime

#self created modules
from utils.utils import stringtoobject, save
import sdata as sdc
from provider.utils import fetch_main_id, sources_table

additional_data_path='../../additional_data/'


def provider_gk(table_names,gk_list_of_tables,distance_cut_in_pc):
    """
    Optains and arranges disk data.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param gk_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :type gk_list_of_tables: list(astropy.table.table.Table)
    :returns: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data and basic disk data.
    :rtype: list(astropy.table.table.Table)
    """
    
    plx_in_mas_cut=1000./distance_cut_in_pc
    #---------------define provider-------------------------------------
    gk_provider=Table()
    gk_provider['provider_name']=['Grant Kennedy Disks']
    gk_provider['provider_url']=['http://drgmk.com/sdb/']
    gk_provider['provider_bibcode']=['priv. comm.']
    gk_provider['provider_access']=['2024-02-09']
    
    print('Creating ',gk_provider['provider_name'][0],' tables ...')
    #loading table obtained via direct communication from Grant Kennedy
    gk_disks=io.votable.parse_single_table(
        additional_data_path+"sdb_30pc_09_02_2024.xml").to_table()
    #transforming from string type into object to have variable length
    gk_disks=stringtoobject(gk_disks,212)
    #removing objects with plx_value=='None' or masked entries
    if len(gk_disks['plx_value'].mask.nonzero()[0])>0:
        print('careful, masked entries in plx_value')
    #gk_disks['plx_value']=gk_disks['plx_value'].filled('None')
    #gk_disks=gk_disks[np.where(gk_disks['plx_value']!='None')]
    #gk_disks['plx_value']=gk_disks['plx_value'].astype(float)
    
    #sorting out everything with plx_value too big
    gk_disks=gk_disks[np.where(gk_disks['plx_value']>plx_in_mas_cut)]
    #adds the column for object type
    gk_disks['type']=['di' for j in range(len(gk_disks))]
    gk_disks['disks_ref']=['Grant Kennedy'
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
    gk_disks=fetch_main_id(gk_disks,colname='gk_host_main_id',name='main_id',
                           oid=False)
    #--------------creating output table gk_h_link ---------------------
    gk_h_link=gk_disks['id','main_id','disks_ref']
    gk_h_link.rename_columns(['main_id','disks_ref'],
                             ['parent_main_id','h_link_ref'])
    gk_h_link.rename_column('id','main_id')
    #--------------creating output table gk_objects --------------------
    gk_disks['ids']=gk_disks['id']#because only single id per source given
    gk_objects=gk_disks['id','ids','type']
    gk_objects.rename_column('id','main_id')
    #--------------creating output table gk_ident ----------------------
    gk_ident=gk_disks['ids','id','disks_ref']
    # would prefer to use id instad of ids paremeter but this raises an error
    # so I use ids which has the same content as id
    gk_ident.rename_columns(['ids','disks_ref'],['main_id','id_ref'])
    #--------------creating output table gk_sources --------------------
    gk_sources=Table()
    tables=[gk_provider,gk_disks]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['disks_ref']]
    for cat,ref in zip(tables,ref_columns):
        gk_sources=sources_table(cat,ref,gk_provider['provider_name'][0],
                                gk_sources)
    #--------------creating output table gk_disk_basic------------------
    gk_disk_basic=gk_disks['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    
    for column in ['rdisk_bb','e_rdisk_bb']:
        if len(gk_disks[column].mask.nonzero()[0])>0:
            print('careful, masked entries in ', column)
        #replacing 'None' with 'nan' as the first one is not float convertible
        #gk_disk_basic=replace_value(gk_disk_basic,column,'None','nan')
        #gk_disk_basic[column].fill_value='nan' #because defeault is None and not float convertible
        #though this poses the issue that the float default float fill_value is 1e20
        #gk_disk_basic[column].filled()
        #gk_disk_basic[column]=gk_disk_basic[column].astype(float)
    gk_disk_basic.rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    gk_disk_basic=gk_disk_basic[np.where(np.isfinite(
                                                gk_disk_basic['rad_value']))]
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': gk_list_of_tables[i]=gk_sources
        if table_names[i]=='provider': gk_list_of_tables[i]=gk_provider
        if table_names[i]=='objects': gk_list_of_tables[i]=gk_objects
        if table_names[i]=='ident': gk_list_of_tables[i]=gk_ident
        if table_names[i]=='h_link': gk_list_of_tables[i]=gk_h_link
        if table_names[i]=='disk_basic': gk_list_of_tables[i]=gk_disk_basic
        save([gk_list_of_tables[i][:]],['gk_'+table_names[i]])
    return gk_list_of_tables