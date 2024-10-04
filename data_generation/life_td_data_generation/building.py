""" Combines the data from the individual data providers. """

import numpy as np #arrays
from astropy.table import Column, join, column, vstack, Table, unique

#self created modules
from utils.io import save, Path
from provider.utils import nullvalues, replace_value
import sdata as sdc

def idsjoin(cat,column_ids1,column_ids2):
    """
    Merges the identifiers from two different columns into one.
    
    :param cat: Astropy table containing two identifer columns.
    :type cat: astropy.table.table.Table
    :param column_ids1: column name for the first identifier column
    :type column_ids1: str
    :param column_ids2: column name for the second identifier column
    :type column_ids2: str
    :return: Like input cat but with only one identifier
        column containing the unique identifiers of both of the
        previous identifier columns.
    :rtype: astropy.table.table.Table
    """
    
    cat['ids']=Column(dtype=object, length=len(cat))
    for column in [column_ids1,column_ids2]:
        cat=nullvalues(cat,column,'')
    for i in range(len(cat)):
        # splitting object into list of elements
        ids1=cat[column_ids1][i].split('|')
        ids2=cat[column_ids2][i].split('|')
        if ids2==['']:
            cat['ids'][i]=cat[column_ids1][i]
        elif ids1==['']:
            cat['ids'][i]=cat[column_ids2][i]
        else:
            ids=ids1+ids2#should be list
            #removing double entries
            ids=set(ids)
            #changing type back into list
            ids=list(ids)
            #joining list into object with elements separated by |
            ids="|".join(ids)
            cat['ids'][i]=ids
        cat['ids'][i]=cat['ids'][i].strip('|')
    return cat

def objectmerging(cat):
    """
    Merges the data of each object given in the different providers.
    
    The object is the same physical one but the data is provided by 
    different providers and merged into one entry.
    
    :param cat: Astropy table containing multiple entries for the same
        physical objects due to data from different providers.
    :type cat: astropy.table.table.Table
    :returns: Like cat with unique object entries.
    :rtype: astropy.table.table.Table
    """
    
    cat=idsjoin(cat,'ids_1','ids_2')
    cat.remove_columns(['ids_1','ids_2'])
    #merging types
    #initializing column
    if 'type' not in cat.colnames:#----------
        cat['type']=Column(dtype=object, length=len(cat))
        cat['type_1']=cat['type_1'].astype(object)
        cat['type_2']=cat['type_2'].astype(object)
        for i in range(len(cat)):
            if type(cat['type_2'][i])==np.ma.core.MaskedConstant or \
                    cat['type_2'][i]=='None':
                cat['type'][i]=cat['type_1'][i]
            else:
                cat['type'][i]=cat['type_2'][i]
        cat.remove_columns(['type_1','type_2'])
    return cat

def match(cat,sources,paras,provider):
    """
    This function joins the source identifiers of parameters of cat.
    
    :param cat: With empty para_source_id columns.
    :type cat: astropy.table.table.Table
    :param sources: Contains reference data.
    :type sources: astropy.table.table.Table
    :param paras: Describes a parameter in cat.
    :type paras: str
    :param provider: Name of the data provider.
    :type provider: str
    :returns: Astropy table containing para_source_id data.
    :rtype: astropy.table.table.Table
    """
    #for all parameters specified
    for para in paras:
        #if they have reference columns
        if para+'_ref' in cat.colnames:
            #if those reference columns are masked
            cat=nullvalues(cat,para+'_ref','None')
            #join to each reference parameter its source_id
            cat=join(cat,sources['ref','source_id'][np.where(
                            sources['provider_name']==provider)],
                            keys_left=para+'_ref',keys_right='ref',
                            join_type='left')
            #renaming column to specify to which parameter the source_id
            # correspond
            cat.rename_column('source_id',f'{para}_source_idref')
            #deleting double column containing reference information
            cat.remove_columns('ref')
            #in case the para_value entry is masked this if environment
            # will put the source_id entry to null
            if para+'_value' in cat.colnames:
                if type(cat[para+'_value'])==column.MaskedColumn:
                    for i in cat[para+'_value'].mask.nonzero()[0]:
                        cat[f'{para}_source_idref'][i]=999999
    return cat

def merge_table(cat1,cat2):
    """
    Merges two tables.
    
    :param cat1:
    :type cat1: astropy.table.table.Table
    :param cat2:
    :type cat2: astropy.table.table.Table
    :returns: Merged table.
    :rtype: astropy.table.table.Table
    """
    
    if len(cat1)==0 or len(cat2)==0:
        #in this case astropy join function wouldn't work
        merged_cat=vstack([cat1,cat2])
    #elif #some columns being empty, others not:
        #remove empty columns
        #make sure merged_cat has all colnames it needs
    else:
        merged_cat=join(cat1,cat2)
    
    return merged_cat

def best_para_id(mes_table):
    para='id'
    best_para_table=mes_table[:0].copy()
    grouped_mes_table=mes_table.group_by('id_ref')
    #making simbad default best para
    mask = grouped_mes_table.groups.keys['id_ref'] == '2000A&AS..143....9W' 
    best_para_table=grouped_mes_table.groups[mask]
    # TBD: use id_ref as variable from provider_bibcode 
    #        instad of constant""")
    for ref in ['2022A&A...664A..21Q','2016A&A...595A...1G','priv. comm.',
                '2020A&C....3100370A','2001AJ....122.3466M']:
        #priority of id best para: 
        mask = grouped_mes_table.groups.keys['id_ref'] == ref
        all_ref_ids=grouped_mes_table.groups[mask]
        #removing those already in best_para_table
        new_ids=all_ref_ids[np.where(np.invert(np.in1d(
                                    all_ref_ids['id'],
                                    best_para_table['id'])))]
        best_para_table=vstack([best_para_table,new_ids])
    best_para_table.remove_column('id_ref')
    return best_para_table

def best_para_membership(mes_table):
    para='membership'
    best_para_table=mes_table[:0].copy()
    grouped_mes_table=mes_table.group_by(['child_object_idref',
                                          'parent_object_idref'])
    ind=grouped_mes_table.groups.indices
    for i in range(len(ind)-1):
        l=ind[i+1]-ind[i]
        if l==1:
            best_para_table.add_row(grouped_mes_table[ind[i]])
        else:
            temp=grouped_mes_table[ind[i]:ind[i+1]]
            not_nan_temp=temp[np.where(temp[para]!=999999)]
            if len(not_nan_temp)>0:
                max_row=not_nan_temp[np.where(
                        not_nan_temp[para]==max(not_nan_temp[para]))]
                for j in range(ind[i],ind[i+1]):
                    if grouped_mes_table[para][j]==max(not_nan_temp[para]):
                        best_para_table.add_row(grouped_mes_table[j])
                        break # make sure not multiple of same max 
                               # value are added
            else:
                #if none of the objects has a membership entry 
                # then pick just first one
                best_para_table.add_row(grouped_mes_table[ind[i]])
    return best_para_table

def best_para(para,mes_table):
    """
    This function keeps only highest quality row for each object. 
    
    tried to avoid . for performance reasons
    
    :param para: Describes parameter e.g. mass
    :type para: str
    :param mes_table: Contains only columns 'main_id',
         para+'_value',para+'_err',para+'_qual' and para+'_ref'
    :type mes_table: astropy.table.table.Table
    :returns: Table like mes_table but only highest quality rows for
         each object
    :rtype: astropy.table.table.Table
    """
    
    if para=='id':
        return best_para_id(mes_table)
    elif para=='membership':
        return best_para_membership(mes_table)
    elif para=='binary':
        columns=['main_id',para+'_flag',para+'_qual',para+'_source_idref']
    elif para=='mass_pl':
        columns=['main_id',para+'_value',para+'_rel',para+'_err',para+'_qual',
                para+'_source_idref']
    elif para=='sep_ang':
        columns=['main_id',para+'_value',para+'_err',para+'_obs_date',
                para+'_qual',para+'_source_idref']
    else:
        columns=['main_id',para+'_value',para+'_err',para+'_qual',
                para+'_source_idref']
    mes_table=mes_table[columns[0:]]
    best_para_table=mes_table[columns[0:]][:0].copy()
    #group mes_table by object (=main_id)
    grouped_mes_table=mes_table.group_by('main_id')
    #take highest quality
    #quite time intensive (few minutes) could maybe be optimized using 
    # np.in1d function
    for j in range(len(grouped_mes_table.groups.keys)):
    # go through all objects
        for qual in ['A','B','C','D','E','?']:
            for i in range(len(grouped_mes_table.groups[j])):
                if grouped_mes_table.groups[j][i][para+'_qual']==qual:
                    best_para_table.add_row(grouped_mes_table.groups[j][i])
                    break # if best para found go to next main_id group
            else:
                continue # only executed if the inner loop did NOT break
            break  # only executed if the inner loop DID break 
    return best_para_table

def best_parameters_ingestion(cat_mes,cat_basic,para,columns=[]):
    """
    Takes from the multiple measurement table best parameters for basic table.

    :para cat_mes:
    :type cat_mes: astropy.table.table.Table
    :para cat_basic:
    :type cat_basic: astropy.table.table.Table
    :para str para: Parameter name
    :para columns: List of column names of the parameter.
    :type columns: list(str)
    :returns:
    :rtype: astropy.table.table.Table
    """
    best_para_cat_mes=best_para(para,cat_mes)
    if columns!=[]:
        cat_basic.remove_columns(columns)
    cat_basic=join(cat_basic,best_para_cat_mes,join_type='left')
    return cat_basic

def provider_data_merging(cat,table_names,table_name,prov_tables_list,o_merging=False,para_match=False):
    """
    Merges the data from the different providers.
    
    :para cat: 
    :type cat: list(astropy.table.table.Table)
    :para table_names:
    :type table_names: list(str)
    :para str table_name: 
    :param prov_tables_list: Containing simbad, grant kennedy, exomercat, gaia 
        and wds data.
    :type prov_tables_list: list(astropy.table.table.Table)
    :para bool o_merging:
    :para bool para_match:
    :returns:
    :rtype: astropy.table.table.Table
    """
    print(f'Building {table_name} table ...')#sources
    n=table_names.index(table_name)
    
    for j in range(len(prov_tables_list)):
        if para_match:
            # redefining of paras multiple times is not optimal but easier to 
            # read the code if I don't have to define them globally or pass
            # them to the function.
            #print('careful, I have here hardcoded parameter and table order')
            paras=[['id'],['h_link'],['coo','plx','dist_st','coo_gal','mag_i',
                                      'mag_j','mag_k','class','sptype'],
                   ['mass_pl'],['rad'],['mass_pl'],['teff_st'],
                   ['radius_st'],['mass_st'],['binary'],['sep_ang'],['h_link']]
            if len(prov_tables_list[j][n])>0:
                #prov_tables_list[j] is a table containing the two columns ref 
                # and provider name. replacing ref columns with 
                # corresponding source_idref one. issue is that order 
                # prov_tables_list and provider_name not the same
                prov_tables_list[j][n]=match(prov_tables_list[j][n],cat[0],paras[n-3],
                                        prov_tables_list[j][2]['provider_name'][0])
        if len(cat[n])>0:
            #joining data from different providers (simbad,...,wds)
            if len(prov_tables_list[j][n])>0:
                if o_merging:
                    cat[n]=join(cat[n],prov_tables_list[j][n],
                                         keys='main_id',join_type='outer')
                    cat[n]=objectmerging(cat[n])
                else:
                    cat[n]=join(cat[n],prov_tables_list[j][n],join_type='outer')
                
        else:
            cat[n]=prov_tables_list[j][n]
    return cat[n]

def unify_null_values(cat,table_names):
    print('Unifying null values...')
    # unify null values (had 'N' and '?' because of astropy default 
    # fill_value and type conversion string vs object)
    tables=[cat[table_names.index('star_basic')],
            cat[table_names.index('planet_basic')],
            cat[table_names.index('disk_basic')],
            cat[table_names.index('mes_mass_pl')],
            cat[table_names.index('mes_teff_st')],
            cat[table_names.index('mes_radius_st')],
            cat[table_names.index('mes_mass_st')],
            cat[table_names.index('mes_binary')]]
    columns=[['coo_qual','coo_gal_qual','plx_qual','dist_st_qual',
              'sep_ang_qual','teff_st_qual','radius_st_qual','binary_flag',
              'binary_qual','mass_st_qual','sptype_qual','class_temp',
              'class_temp_nr'],
             ['mass_pl_qual','mass_pl_rel'],
             ['rad_qual','rad_rel'],['mass_pl_qual','mass_pl_rel'],
             ['teff_st_qual'],['radius_st_qual'],['mass_st_qual'],
             ['binary_qual']]
    for i in range(len(tables)):
        for col in columns[i]:
            tables[i]=replace_value(tables[i],col,'N','?')
            tables[i]=replace_value(tables[i],col,'N/A','?')
    return cat

def initialize_empty_tables(list_of_tables):
    """
    Creates empty tables as needed for final database ingestion
    """
    
    empty=sdc.provider('empty')
    n_tables=len(empty.list_of_tables)
    cat=[Table() for i in range(n_tables)]    
    return empty, n_tables, cat

#------------------------provider combining----------------------------
def building(prov_tables_list,table_names,list_of_tables):
    """
    This function builds the tables for the LIFE database.
    
    :param prov_tables_list: Containing simbad, grant kennedy, exomercat, gaia 
        and wds data.
    :type prov_tables_list: list(astropy.table.table.Table)
    :param table_names: Objects correspond to the names of the 
        astropy tables contained in prov_tables_list and the return list.
    :type table_names: list of objects of type str
    :param list_of_tables: Empty astropy tables to be filled
        in and returned.
    :type list_of_tables: list(astropy.table.table.Table)
    :returns: Containing data combined from the different prov_tables_list.
    :rtype: list containing objects of type astropy.table.table.Table
    """
    
    empty, n_tables, cat = initialize_empty_tables(list_of_tables)
    
    #for the sources and objects joins tables from different prov_tables_list
    cat[0]=provider_data_merging(cat,table_names,'sources',prov_tables_list)
    
    #adding empty columns for later being able to join tables
    cat[0]=vstack([cat[0],empty.table('sources')])
    
    # keeping only unique values then create identifiers for the tables
    if len(cat[0])>0:#not needed as never no sources given
        cat[0]=unique(cat[0],silent=True)
        cat[0]['source_id']=[j+1 for j in range(len(cat[0]))]
    
    cat[1]=provider_data_merging(cat,table_names,'objects',prov_tables_list,o_merging=True)
                    
    #assigning object_id
    cat[1]['object_id']=[j+1 for j in range(len(cat[1]))]

    # At one point I would like to be able to merge objects with main_id
    # NAME Proxima Centauri b and Proxima Centauri b

    cat[2]=provider_data_merging(cat,table_names,'provider',prov_tables_list)
    
    #I do this to get those columns that are empty in the data
    cat[2]=vstack([cat[2],empty.table('provider')])
       
    
    for i in range(3,n_tables): 
    # for the tables star_basic,...,mes_mass_st
        cat[i]=provider_data_merging(cat,table_names,table_names[i],
                                     prov_tables_list,para_match=True)

        #I do this to get those columns that are empty in the data
        cat[i]=vstack([cat[i],empty.table(table_names[i])])
        cat[i]=cat[i].filled() 
        # because otherwise unique does neglect masked columns

        if 'object_idref' in cat[i].colnames and len(cat[i])>0: 
            # add object_idref
            # first remove the object_idref we got from empty 
            # initialization though I would prefer a more elegant way 
            # to do this. Is needed as empty columns don't work for join
            cat[i].remove_column('object_idref') 
            cat[i]=join(cat[i],cat[1]['object_id','main_id'],
                                 join_type='left')
            cat[i].rename_column('object_id','object_idref')
        if table_names[i]=='ident':
            cat[i]=best_para('id',cat[i])
        if table_names[i]=='h_link':
            #expanding from child_main_id to object_idref
            #first remove the child_object_idref we got from empty
            # initialization. Would prefer a more elegant way to do this
            cat[i].remove_column('child_object_idref')
            cat[i]=join(cat[i],cat[1]['object_id','main_id'],
                                 keys='main_id',join_type='left')
            cat[i].rename_columns(['object_id','main_id'],
                                  ['child_object_idref','child_main_id'])
            
            #expanding from parent_main_id to parent_object_idref
            cat[i].remove_column('parent_object_idref')
            #kick out any h_link rows where parent_main_id not in
            # objects (e.g. clusters)
            cat[i]=join(cat[i],cat[1]['object_id','main_id'],
                   keys_left='parent_main_id',keys_right='main_id')
            #removing because same as parent_main_id
            cat[i].remove_column('main_id')
            cat[i].rename_column('object_id','parent_object_idref')
            cat[table_names.index('best_h_link')]=best_para('membership',cat[table_names.index('h_link')])
        if table_names[i]=='star_basic':
            #choosing all objects with type star or system. this I use 
            # to join the object_id parameter from objects table to 
            # star_basic. what about gaia stuff where I don't know 
            # this? there I also don't have star_basic info.
            # Note: main_id was only added because I have not found out 
            # how to do join with just one column of a table
            stars=cat[1]['object_id','main_id'][np.where(
                            cat[1]['type']=='st')]
            systems=cat[1]['object_id','main_id'][np.where(
                            cat[1]['type']=='sy')]
            temp=vstack([stars,systems])
            temp.rename_column('object_id','object_idref')
            # cat[i] are all the star_cat tables from prov_tables_list where 
            # those are given the new objects are needed to join the 
            # best parameters from mes_ tables later on
            cat[i]=join(cat[i],temp,join_type='outer',
                                 keys=['object_idref','main_id'])
        if table_names[i]=='planet_basic':
            planets=cat[1]['object_id','main_id'][np.where(
                            cat[1]['type']=='pl')]
            planets.rename_column('object_id','object_idref')
            cat[i]=planets #can't use join below because cat[i] has no rows
        if table_names[i]=='mes_teff_st':
            cat[table_names.index('star_basic')]=best_parameters_ingestion(
                    cat[i], cat[table_names.index('star_basic')],
                    'teff_st',['teff_st_value','teff_st_err',
                               'teff_st_qual','teff_st_source_idref',
                               'teff_st_ref'])
        if table_names[i]=='mes_radius_st':
            cat[table_names.index('star_basic')]=best_parameters_ingestion(
                    cat[i], cat[table_names.index('star_basic')],
                    'radius_st',['radius_st_value','radius_st_err',
                                 'radius_st_qual','radius_st_source_idref',
                                 'radius_st_ref'])
        if table_names[i]=='mes_mass_st':
            cat[table_names.index('star_basic')]=best_parameters_ingestion(
                    cat[i], cat[table_names.index('star_basic')],
                    'mass_st',['mass_st_value','mass_st_err',
                               'mass_st_qual','mass_st_source_idref',
                               'mass_st_ref'])
        if table_names[i]=='mes_mass_pl':
            cat[table_names.index('planet_basic')]=best_parameters_ingestion(
                    cat[i], cat[table_names.index('planet_basic')],
                    'mass_pl')
        if table_names[i]=='mes_binary':
            cat[table_names.index('star_basic')]=best_parameters_ingestion(
                    cat[i], cat[table_names.index('star_basic')],
                    'binary',['binary_flag',
                              'binary_qual','binary_source_idref',
                               'binary_ref'])
        if table_names[i]=='mes_sep_ang':
            cat[table_names.index('star_basic')]=best_parameters_ingestion(
                    cat[i], cat[table_names.index('star_basic')],
                    'sep_ang',['sep_ang_value','sep_ang_err',
                               'sep_ang_obs_date','sep_ang_qual',
                               'sep_ang_source_idref','sep_ang_ref'])
            
        cat[i]=cat[i].filled()
        
        if len(cat[i])==0:
            print('warning: empty table',i,table_names[i])
        else:
            #only keeping unique entries
            cat[i]=unique(cat[i],silent=True)
            
    #next line is needed as multimeasurement adaptions lead to potentially masked entries
    cat[table_names.index('star_basic')]=cat[table_names.index('star_basic')].filled()
    
    cat=unify_null_values(cat,table_names)
            
    # TBD: Add exact object distance cut. So far for correct treatment
    #       of boundary objects 10% additional distance cut used""")
    
    print('Saving data...')
    save(cat,table_names,location=Path().data)
    return cat
