""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from astropy.table import Table, unique, vstack, join
from datetime import datetime

#self created modules
from utils.io import save, load
from provider.utils import fill_sources_table, create_sources_table, query, ids_from_ident, distance_cut, create_provider_table
from sdata import empty_cat

def create_wds_helpertable(wds,temp,test_objects):    
    # define queries
    adql_query=["""SELECT
                    WDS as wds_name, Comp as wds_comp,
                    sep1 as wds_sep1, sep2 as wds_sep2, 
                    Obs1 as wds_obs1, Obs2 as wds_obs2
                    FROM "B/wds/wds" """]
    
    print('Creating ',wds['provider']['provider_name'][0],' tables ...')
    #perform query for objects with parallax >50mas
    test_objects=np.array(test_objects)
    if temp:
        print(' loading...')
        [wds_helptab]=load(['wds_helptab'])
        # currently temp=True not giving same result because 
        # wds['system_main_id'][j] are '' and not masked
        for col in ['system_main_id','primary_main_id','secondary_main_id']:
            wds_helptab[col][np.where(wds_helptab[col]=='')]=np.ma.masked
        # tbd: add provider_access of last query
    else:
        print(' querying VizieR for WDS...')
        wds_helptab=query(wds['provider']['provider_url'][0],adql_query[0])
        
        # I need to match the wds objects with the simbad ones to inforce the
        # distance cut since wds does not have distance information.
        
        # initializing and setting type for object comparison in later join 
        for col in ['sim_wds_id','system_name','primary','secondary']:
            wds_helptab[col]=wds_helptab['wds_name'].astype(object)
        
        # assigning correct name of system, primary and secondary for each wds object
        for j in range(len(wds_helptab)):
            if wds_helptab['wds_comp'][j]=='':#trivial binaries
                wds_helptab['system_name'][j]='WDS J'+wds_helptab['wds_name'][j]+'AB'
                #AB added since apparently simbad calls trivial binary system AB too
                wds_helptab['primary'][j]='WDS J'+wds_helptab['wds_name'][j]+'A'
                wds_helptab['secondary'][j]='WDS J'+wds_helptab['wds_name'][j]+'B'
            else:#higer order multiples
                wds_helptab['system_name'][j]='WDS J'+wds_helptab['wds_name'][j]+wds_helptab['wds_comp'][j]
                if len(wds_helptab['wds_comp'][j])==2:
                    wds_helptab['primary'][j]='WDS J'+wds_helptab['wds_name'][j]+wds_helptab['wds_comp'][j][0]
                    wds_helptab['secondary'][j]='WDS J'+wds_helptab['wds_name'][j]+wds_helptab['wds_comp'][j][1]
                else:
                    components=wds_helptab['wds_comp'][j].split(',')
                    wds_helptab['primary'][j]='WDS J'+wds_helptab['wds_name'][j]+components[0]
                    wds_helptab['secondary'][j]='WDS J'+wds_helptab['wds_name'][j]+components[1]
        # print('number of trivial binary systems:',
        #   len(wds[np.where(wds['wds_comp']=='')]))
                
        if len(test_objects)>0:
            print('in wds as system_name', 
                    test_objects[np.where(np.in1d(test_objects,
                                                    wds_helptab['system_name']))])
            print('in wds as primary',
                    test_objects[np.where(np.in1d(test_objects,
                                                    wds_helptab['primary']))])
            print('in wds as secondary',
                    test_objects[np.where(np.in1d(test_objects,
                                                    wds_helptab['secondary']))])


    # an alternative would be to query simbad for the main id and then cut by distance
    # this however takes way longer as it joins 150'000 elements
    #    wds=fetch_main_id(wds,colname='wds_full_name',name='main_id',oid=False)
    #    wds=distance_cut(wds,colname='wds_full_name',main_id=True)
        print(' performing distance cut...')
        
        #assigning main_id for system using sim_hlink and cutting on the system
        #  or the components
        wds_system_cut=distance_cut(wds_helptab,colname='system_name',main_id=False)
        wds_system_cut.rename_column('main_id','system_main_id')
        
        wds_primary_cut=distance_cut(wds_helptab,colname='primary',main_id=False)
        
        wds_secondary_cut=distance_cut(wds_helptab,colname='secondary',main_id=False)
        [sim_h_link]=load(['sim_h_link'])
        #joining parent object
        wds_primary_cut=join(wds_primary_cut,sim_h_link['main_id','parent_main_id'],
                                  keys='main_id',join_type='left')
        wds_primary_cut.rename_columns(['main_id','parent_main_id'],['primary_main_id','system_main_id'])
        
        wds_secondary_cut=join(wds_secondary_cut,sim_h_link['main_id','parent_main_id'],
                                  keys='main_id',join_type='left')
        wds_secondary_cut.rename_columns(['main_id','parent_main_id'],['secondary_main_id','system_main_id'])
        #here some empty ones when child is known in simbad but no parent. in 
        # this case would I want to assign system_name in system main_id? do it 
        # later
        
        wds_helptab=vstack([wds_system_cut,wds_primary_cut])
        wds_helptab=vstack([wds_helptab,wds_secondary_cut])
                        
        if len(test_objects)>0:
            print(wds_helptab['system_main_id','primary_main_id','secondary_main_id'])
            print('in wds as system_main_id',
                    test_objects[np.where(np.in1d(test_objects,
                                                    wds_helptab['system_main_id']))])
            print('in wds as primary_main_id',
                    test_objects[np.where(np.in1d(test_objects,
                                                    wds_helptab['primary_main_id']))])
            print('in wds as secondary_main_id',
                    test_objects[np.where(np.in1d(test_objects,
                                                    wds_helptab['secondary_main_id']))])
        
        save([wds_helptab],['wds_helptab'])
    
    wds_helptab['system_main_id']=wds_helptab['system_main_id'].astype(object)
    wds_helptab['system_name']=wds_helptab['system_name'].astype(object)

    return wds_helptab

def create_ident_and_h_link_table(wds_helptab,wds,test_objects):
    
    #-----------------creating output table wds_ident and wds_h_link------------
    wds_ident=Table(names=['main_id','id'],dtype=[object,object],masked=True)
    # create wds_h_link (for systems)
    wds_h_link=Table(names=['main_id','parent_main_id'],dtype=[object,object])
    #add all relevant invormation
    #about identifiers
    table_main=['system_name','system_main_id','system_main_id',
               'primary','primary_main_id','primary_main_id',
               'secondary','secondary_main_id','secondary_main_id']
    table_id=['system_name','system_main_id','system_name',
             'primary','primary_main_id','primary',
             'secondary','secondary_main_id','secondary']
    empty=Table(names=['main_id'],dtype=[object],masked=True)
    for id1,id2 in zip(table_main,table_id):
        temp=empty.copy()
        temp['main_id']=wds_helptab[id1].astype(object)
        temp['id']=wds_helptab[id2].astype(object)
        wds_ident=vstack([wds_ident,temp])
    
    #about relations of objects
    table_main_id=['primary','primary','primary_main_id','primary_main_id',
                   'secondary','secondary','secondary_main_id','secondary_main_id']
    table_parent=['system_name','system_main_id','system_name','system_main_id',
                  'system_name','system_main_id','system_name','system_main_id']
    for id1,id2 in zip(table_main_id,table_parent):
        temp=empty.copy()
        temp['main_id']=wds_helptab[id1].astype(object)
        temp['parent_main_id']=wds_helptab[id2].astype(object)
        wds_h_link=vstack([wds_h_link,temp])
    #delete all rows containing masked entries
    wds_ident.remove_rows(wds_ident['main_id'].mask.nonzero()[0])
    wds_ident.remove_rows(wds_ident['id'].mask.nonzero()[0])
    wds_h_link.remove_rows(wds_h_link['main_id'].mask.nonzero()[0])
    wds_h_link.remove_rows(wds_h_link['parent_main_id'].mask.nonzero()[0])

    #uniqueness
    wds_ident=unique(wds_ident)
    wds_h_link=unique(wds_h_link)

    #delete entries where id instead of main_id used
    not_identical_rows_id=wds_ident['id'][np.where(wds_ident['main_id']!=wds_ident['id'])]
    remove=np.in1d(wds_ident['main_id'],not_identical_rows_id)
    wds_ident.remove_rows(remove)

    #for h_link replacing instead of deleting because there can be cases where the information I need is only available this way
    #e.g. simbad query on system name results in system main_id. simbad query on primary gives primary_main_id. but 
    #simbad does not have those two objects connected through h_link. therefore I only have [primary,system_main_id] but
    #would want to have [primary_main_id,system_main_id]

    #where h_link main_id not in ident_main_id
    not_main_id=np.invert(np.in1d(wds_h_link['main_id'],wds_ident['main_id']))
    
    if len(test_objects)>0:
        print('number of test objects that are in h_link main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_h_link['main_id']))])
        print('number of test objects that are in h_link parent_main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_h_link['parent_main_id']))])
        print('number of test objects that are in main_id of ident table \n', \
              test_objects[np.where(np.in1d(test_objects,wds_ident['main_id']))])
        print('number of test objects that are in main_id of h_link but not ident table \n', \
              test_objects[np.where(np.in1d(test_objects,wds_h_link['main_id'][not_main_id]))])

    #replace it with the corresponding  ident main_id
    for j in range(len(wds_h_link['main_id'][not_main_id])):
        temp=wds_h_link['main_id'][not_main_id][j]
        #why is this line not working? am I assigning stuff to a copy instead of alias?
        wds_h_link['main_id'][np.where(wds_h_link['main_id']==temp)]=wds_ident['main_id'][np.where(
                wds_ident['id']==temp)]
    
    if len(test_objects)>0:
        print('number of test objects that are in h_link main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_h_link['main_id']))])
        print('number of test objects that are in main_id of ident table \n', \
              test_objects[np.where(np.in1d(test_objects,wds_ident['main_id']))])
    
    #where h_link parent_main_id not in ident_main_id
    not_parent_main_id=np.invert(np.in1d(wds_h_link['parent_main_id'],wds_ident['main_id']))

    #replace it with the corresponding  ident main_id
    for j in range(len(wds_h_link['parent_main_id'][not_parent_main_id])):
        if len(wds_ident['main_id'][np.where(
                    wds_ident['id']==wds_h_link['parent_main_id'][not_parent_main_id][j])])==1:
            temp=wds_h_link['parent_main_id'][not_parent_main_id][j]
            wds_h_link['parent_main_id'][np.where(wds_h_link['parent_main_id']==temp)]=wds_ident['main_id'][np.where(
                    wds_ident['id']==temp)]
        #else:
        #nestled multiples with non hierarchical measurements e.g. AC component when A and B are closest and C further away
            #print(wds_ident['main_id'][np.where(
            #        wds_ident['id']==wds_h_link['parent_main_id'][not_parent_main_id][j])])
            #print(wds_h_link['parent_main_id'][not_parent_main_id][j])
            #print(wds_h_link['main_id'][not_parent_main_id][j])
            #wds_h_link['parent_main_id'][not_parent_main_id][j]=wds_ident['main_id'][np.where(
            #        wds_ident['id']==wds_h_link['parent_main_id'][not_parent_main_id][j])][0]

    wds_h_link=unique(wds_h_link) 

    #refs
    wds_ident['id_ref']=[wds['provider']['provider_bibcode'][0] for j in range(len(wds_ident))]
    wds_h_link['h_link_ref']=[wds['provider']['provider_bibcode'][0] for j in range(len(wds_h_link))]    
    return wds_ident,wds_h_link

def create_objects_table(wds_helptab,wds,test_objects):
    #create ids
    wds_objects=Table(names=['main_id','ids'],dtype=[object,object])
    wds_objects=ids_from_ident(wds['ident']['main_id','id'],wds_objects)
    #if it has children, it is type system
    #if it has no children it can either be star or close in system
    wds_objects['type']=['sy' for j in range(len(wds_objects))]
    #change to st for those that have no children
    wds_objects['type'][np.invert(np.in1d(wds_objects['main_id'],wds['h_link']['parent_main_id']))]=['st' for j in range(len(
            [np.invert(np.in1d(wds_objects['main_id'],wds['h_link']['parent_main_id']))]))] 
    
    if len(test_objects)>0:
        print('number of test objects that are in objects main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_objects['main_id']))])    
    return wds_objects

def create_mes_binary_table(wds_helptab,wds,test_objects):
    wds_mes_binary=wds['objects']['main_id','type']#[np.where(wds_objects['type']=='sy')]
    wds_mes_binary.rename_column('type','binary_flag')
    wds_mes_binary['binary_flag']=wds_mes_binary['binary_flag'].astype(object)
    wds_mes_binary['binary_flag']=['True' for j in range(len(wds_mes_binary))]
    wds_mes_binary['binary_ref']=[wds['provider']['provider_bibcode'][0] for j in range(len(wds_mes_binary))]
    wds_mes_binary['binary_qual']=['C' for j in range(len(wds_mes_binary))]
    
    if len(test_objects)>0:
        print('number of test objects that are in mes_binary main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_mes_binary['main_id']))])    
    return wds_mes_binary

def create_mes_sep_ang_table(wds_helptab,wds,test_objects):
    #-----------------creating output table wds_mes_sep_ang------------------------
    #better join them
    wds_mes_sep_ang0=join(wds_helptab['system_name','wds_sep1','wds_obs1','wds_sep2','wds_obs2'],
                                  wds['ident']['main_id','id'],keys_left='system_name', keys_right='id')
    #replacing empty system_main_id with main_id from ident using system_name column
    #masked_system_main_id=wds['system_main_id'].mask.nonzero()[0]
    #for j in range(len(wds[masked_system_main_id])):
     #   wds['system_main_id'][masked_system_main_id][j]=wds_ident['main_id'][np.where(
      #          wds_ident['id']==wds['system_name'][masked_system_main_id][j])]
    wds_mes_sep_ang1=wds_mes_sep_ang0['main_id','wds_sep1','wds_obs1']
    wds_mes_sep_ang1.rename_columns(['wds_sep1','wds_obs1'],['sep_ang_value','sep_ang_obs_date'])
    wds_mes_sep_ang1['sep_ang_qual']= \
            ['C' if type(j)!=np.ma.core.MaskedConstant else 'E' for j in wds_mes_sep_ang1['sep_ang_obs_date']]
    #issue, what if system_main_id is empty?
    
    
    wds_mes_sep_ang2=wds_mes_sep_ang0['main_id','wds_sep2','wds_obs2']
    wds_mes_sep_ang2.rename_columns(['wds_sep2','wds_obs2'],['sep_ang_value','sep_ang_obs_date'])
    wds_mes_sep_ang2['sep_ang_qual']= \
            ['B' if type(j)!=np.ma.core.MaskedConstant else 'E' for j in wds_mes_sep_ang2['sep_ang_obs_date']]
    wds_mes_sep_ang=vstack([wds_mes_sep_ang1,wds_mes_sep_ang2])
    #add a quality to sep1 which is better than sep2. because newer measurements should be better.
    wds_mes_sep_ang['sep_ang_ref']=[wds['provider']['provider_bibcode'][0] for j in range(len(wds_mes_sep_ang))]
    #wds_mes_sep_ang.rename_column('system_main_id','main_id')
    #remove columns where sep_ang_value is masked
    wds_mes_sep_ang.remove_columns(wds_mes_sep_ang['sep_ang_value'].mask.nonzero()[0])
    #uniqueness where obs date not known 
    if len(wds_mes_sep_ang['sep_ang_obs_date'].mask.nonzero()[0])>0:
        unique_unknown_obs_date=unique(wds_mes_sep_ang[np.where(
                wds_mes_sep_ang['sep_ang_obs_date'].mask.nonzero()[0])],keys=['main_id','sep_ang_value'])
        unique_known_obs_date=unique(wds_mes_sep_ang[np.where(
                np.invert(wds_mes_sep_ang['sep_ang_obs_date'].mask.nonzero()[0]))],
                keys=['main_id','sep_ang_value','sep_ang_obs_date'])   
        wds_mes_sep_ang=vstack([unique_unknown_obs_date,unique_known_obs_date])
    else:
        wds_mes_sep_ang=unique(wds_mes_sep_ang)
        
    if len(test_objects)>0:
        print('number of test objects that are in mes_sep_ang main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_mes_sep_ang['main_id']))])

    return wds_mes_sep_ang

def create_wds_sources_table(wds):
    tables=[wds['provider'],wds['ident']]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['id_ref']]
    wds_sources=create_sources_table(tables,ref_columns,
                    wds['provider']['provider_name'][0])
    return wds_sources 

def provider_wds(temp=False,test_objects=[]):
    """
    This function obtains and arranges wds data.
    
    :param bool temp: Defaults to False. Used for debugging. saves querrying time.
    :returns: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data, basic stellar data and binarity data.
    :rtype:  list(astropy.table.table.Table)
    """
    wds = empty_cat.copy()
    wds['provider'] = create_provider_table('WDS',
                            'http://tapvizier.u-strasbg.fr/TAPVizieR/tap',
                             '2001AJ....122.3466M')

    wds_helptab = create_wds_helpertable(wds,temp,test_objects)
    wds['ident'],wds['h_link']=create_ident_and_h_link_table(wds_helptab,
                                                             wds,test_objects)    
    wds['objects']=create_objects_table(wds_helptab,wds,test_objects)
    wds['mes_binary']=create_mes_binary_table(wds_helptab,wds,test_objects)    
    wds['mes_sep_ang']=create_mes_sep_ang_table(wds_helptab,wds,test_objects)
    wds['sources']=create_wds_sources_table(wds)

    save(list(wds.values()),['wds_'+ element for element in list(wds.keys())])
    return wds
