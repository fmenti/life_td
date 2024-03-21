import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
from datetime import datetime
import importlib #reloading external functions after modification

#self created modules
import helperfunctions as hf
importlib.reload(hf)#reload module after changing it
import provider as p
importlib.reload(p)#reload module after changing it

#-------------------initialization function------------------------------------
def initialize_database_tables(table_names,list_of_tables):
    """
    This function initializes the database tables.
    
    It does so with column name and data
    type specified but no actual data in them.
    :param table_names:
    :type table_names: list of objects of type str
    :param list_of_tables:
    :type list_of_tables:
    :returns: Initiated database tables in the order of table_names.
    :rtype: list(astropy.table.table.Table)
    """
    
    # Explanation of abbreviations: id stands for identifier, idref for
    # reference identifier and parameter_source_idref for the identifier in the
    # source table corresponding to the mentioned parameter

    sources=ap.table.Table(
        #reference,...
        names=['ref','provider_name',
               'source_id'],
        dtype=[object,object,int])

    objects=ap.table.Table(
        #object id, type of object, all identifiers, main id
        names=['object_id','type','ids','main_id'],
        dtype=[int,object,object,object])
    
    provider=ap.table.Table(
        #reference,...
        names=['provider_name','provider_url','provider_bibcode',
               'provider_access'],
        dtype=[object,object,object,object])

    #identifier table
    ident=ap.table.Table(
        #object idref, id, source idref for the id parameter
        names=['object_idref','id','id_source_idref'],
        dtype=[int,object,int])

    #hierarchical link table (which means relation between objects)
    h_link=ap.table.Table(
        #child object idref (e.g. planet X),
        #parent object idref (e.g. host star of planet X)
        #source idref of h_link parameter, h_link reference,
        #membership probability
        names=['child_object_idref','parent_object_idref',
               'h_link_source_idref','h_link_ref','membership'],
        dtype=[int,int,int,object,int])

    star_basic=ap.table.Table(
        #object idref, RA coordinate, DEC coordinate,
        #coordinate error ellypse angle, major axis and minor axis,
        #coordinate quality, source idref of coordinate parameter,
        #coordinate reference, parallax value, parallax error, parallax quality
        #source idref of parallax parameter ... same for distance parameter
        names=['object_idref','coo_ra','coo_dec','coo_err_angle',
               'coo_err_maj','coo_err_min','coo_qual',
               'coo_source_idref','coo_ref',
               'coo_gal_l','coo_gal_b','coo_gal_err_angle',
               'coo_gal_err_maj','coo_gal_err_min','coo_gal_qual',
               'coo_gal_source_idref','coo_gal_ref',
               'mag_i_value','mag_i_err','mag_i_qual','mag_i_source_idref',
               'mag_i_ref',
               'mag_j_value','mag_j_err','mag_j_qual','mag_j_source_idref',
               'mag_j_ref',
               'mag_k_value','mag_k_err','mag_k_qual','mag_k_source_idref',
               'mag_k_ref',
               'plx_value','plx_err','plx_qual','plx_source_idref',
               'plx_ref',
               'dist_st_value','dist_st_err','dist_st_qual','dist_st_source_idref',
               'dist_st_ref',
               'sptype_string','sptype_err','sptype_qual','sptype_source_idref',
               'sptype_ref',
               'class_temp','class_temp_nr','class_lum','class_source_idref',
               'class_ref',
               'teff_st_value','teff_st_err','teff_st_qual','teff_st_source_idref',
               'teff_st_ref',
               'radius_st_value','radius_st_err','radius_st_qual',
               'radius_st_source_idref','radius_st_ref',
               'mass_st_value','mass_st_err','mass_st_qual','mass_st_source_idref',
               'mass_st_ref',
               'binary_flag','binary_qual','binary_source_idref','binary_ref',
               'sep_ang_value','sep_ang_err','sep_ang_obs_date','sep_ang_qual',
               'sep_ang_source_idref','sep_ang_ref'],
        dtype=[int,float,float,float,#coo
               float,float,object,
               int,object,
               float,float,float,#coo_gal
               float,float,object,
               int,object,
               float,float,object,int,#mag_i
               object,
               float,float,object,int,#mag_j
               object,
               float,float,object,int,#mag_k
               object,
               float,float,object,int,#plx
               object,
               float,float,object,int,#dist
               object,
               object,float,object,int,#sptype
               object,
               object,object,object,int,#class
               object,
               float,float,object,int,#teff
               object,
               float,float,object,#rad
               int,object,
               float,float,object,int,#mass
               object,
               object,object,int,object,#binary
               float,float,int,object,#sep_ang
               int,object])
    
    planet_basic=ap.table.Table(
        #object idref, mass value, mass error, mass realtion (min, max, equal),
        #mass quality, source idref of mass parameter, mass reference
        names=['object_idref','mass_pl_value','mass_pl_err','mass_pl_rel',
               'mass_pl_qual','mass_pl_source_idref','mass_pl_ref'],
        dtype=[int,float,float,object,object,int,object])

    disk_basic=ap.table.Table(
        #object idref, black body radius value, bbr error,
        #bbr relation (min, max, equal), bbr quality,...
        names=['object_idref','rad_value','rad_err','rad_rel','rad_qual',
               'rad_source_idref','rad_ref'],
        dtype=[int,float,float,object,object,int,object])

    mes_mass_pl=ap.table.Table(
        names=['object_idref','mass_pl_value','mass_pl_err','mass_pl_rel',
               'mass_pl_qual','mass_pl_source_idref','mass_pl_ref'],
        dtype=[int,float,float,object,object,int,object])
    
    mes_teff_st=ap.table.Table(
        names=['object_idref','teff_st_value','teff_st_err','teff_st_qual',
               'teff_st_source_idref','teff_st_ref'],
        dtype=[int,float,float,object,int,object])
    
    mes_radius_st=ap.table.Table(
        names=['object_idref','radius_st_value','radius_st_err',
               'radius_st_qual','radius_st_source_idref','radius_st_ref'],
        dtype=[int,float,float,object,int,object])
    
    mes_mass_st=ap.table.Table(
        names=['object_idref','mass_st_value','mass_st_err','mass_st_qual',
               'mass_st_source_idref','mass_st_ref'],
        dtype=[int,float,float,object,int,object])
    
    mes_binary=ap.table.Table(
        names=['object_idref','binary_flag','binary_qual',
               'binary_source_idref','binary_ref'],
        dtype=[int,object,object,int,object])
    
    mes_sep_ang=ap.table.Table(
        names=['object_idref',
               'sep_ang_value','sep_ang_err',
               'sep_ang_obs_date','sep_ang_qual',
               'sep_ang_source_idref','sep_ang_ref'],
        dtype=[int,
               float,float,
               int,object,
               int,object])
    
    best_h_link=ap.table.Table(
        #child object idref (e.g. planet X),
        #parent object idref (e.g. host star of planet X)
        #source idref of h_link parameter, h_link reference,
        #membership probability
        names=['child_object_idref','parent_object_idref',
               'h_link_source_idref','h_link_ref','membership'],
        dtype=[int,int,int,object,int])

    for i in range(len(table_names)):
        if table_names[i]=='sources': list_of_tables[i]=sources
        if table_names[i]=='provider': list_of_tables[i]=provider
        if table_names[i]=='objects': list_of_tables[i]=objects
        if table_names[i]=='ident': list_of_tables[i]=ident
        if table_names[i]=='h_link': list_of_tables[i]=h_link
        if table_names[i]=='star_basic': list_of_tables[i]=star_basic
        if table_names[i]=='planet_basic': list_of_tables[i]=planet_basic
        if table_names[i]=='disk_basic': list_of_tables[i]=disk_basic
        if table_names[i]=='mes_mass_pl': list_of_tables[i]=mes_mass_pl
        if table_names[i]=='mes_teff_st': list_of_tables[i]=mes_teff_st
        if table_names[i]=='mes_radius_st': list_of_tables[i]=mes_radius_st
        if table_names[i]=='mes_mass_st': list_of_tables[i]=mes_mass_st
        if table_names[i]=='mes_binary': list_of_tables[i]=mes_binary
        if table_names[i]=='mes_sep_ang': list_of_tables[i]=mes_sep_ang
        if table_names[i]=='best_h_link': list_of_tables[i]=best_h_link
        hf.save([list_of_tables[i][:]],['empty_'+table_names[i]])
    return list_of_tables

def idsjoin(cat,column_ids1,column_ids2):
    """
    This function merges the identifiers from two different columns into one.
    
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
    
    cat['ids']=ap.table.Column(dtype=object, length=len(cat))
    for column in [column_ids1,column_ids2]:
        cat=p.nullvalues(cat,column,'')
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
    This function merges the data of one object from different providers.
    
    The object is the same physical one but the data is provided by different 
    providers and merged into one entry.
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
        cat['type']=ap.table.Column(dtype=object, length=len(cat))
        cat['type_1']=cat['type_1'].astype(object)
        cat['type_2']=cat['type_2'].astype(object)
        for i in range(len(cat)):
            if type(cat['type_2'][i])==np.ma.core.MaskedConstant or cat['type_2'][i]=='None':
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
            cat=p.nullvalues(cat,para+'_ref','None')
            #join to each reference parameter its source_id
            cat=ap.table.join(cat,sources['ref','source_id'][np.where(
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
                if type(cat[para+'_value'])==ap.table.column.MaskedColumn:
                    for i in cat[para+'_value'].mask.nonzero()[0]:
                        cat[f'{para}_source_idref'][i]=999999
    return cat

def best_para(para,mes_table):
    """
    This function keeps only highest quality row for each object. 
    
    tried to avoid . for performance reasons
    :param para: Describes parameter e.g. mass
    :type para: str
    :param mes_table: Contains only columns
        'main_id', para+'_value',para+'_err',para+'_qual' and para+'_ref'
    :type mes_table: astropy.table.table.Table
    :returns: Table like mes_table but only highest quality rows for each object
    :rtype: astropy.table.table.Table
    """
    
    if para=='id':
        best_para=mes_table[:0].copy()
        grouped_mes_table=mes_table.group_by('id_ref')
        mask = grouped_mes_table.groups.keys['id_ref'] == '2000A&AS..143....9W'# sim
        best_para=grouped_mes_table.groups[mask]
        print('TBD: use id_ref as variable from provider_bibcode instad of constant')
        for ref in ['2022A&A...664A..21Q','2016A&A...595A...1G','priv. comm.',
                    '2020A&C....3100370A','2001AJ....122.3466M']:
            mask = grouped_mes_table.groups.keys['id_ref'] == ref
            new_ids=grouped_mes_table.groups[mask][np.where(np.invert(np.in1d(
                                        grouped_mes_table.groups[mask]['id'],
                                        best_para['id'])))]
        best_para=ap.table.vstack([best_para,new_ids])
        best_para.remove_column('id_ref')
        return best_para
    elif para=='membership':
        best_para=mes_table[:0].copy()
        grouped_mes_table=mes_table.group_by(['child_object_idref',
                                              'parent_object_idref'])
        ind=grouped_mes_table.groups.indices
        for i in range(len(ind)-1):
            l=ind[i+1]-ind[i]
            if l==1:
                best_para.add_row(grouped_mes_table[ind[i]])
            else:
                temp=grouped_mes_table[ind[i]:ind[i+1]]
                not_nan_temp=temp[np.where(temp[para]!=999999)]
                if len(not_nan_temp)>0:
                    max_row=not_nan_temp[np.where(
                            not_nan_temp[para]==max(not_nan_temp[para]))]
                    for j in range(ind[i],ind[i+1]):
                        if grouped_mes_table[para][j]==max(not_nan_temp[para]):
                            best_para.add_row(grouped_mes_table[j])
                            break#make sure not multiple of same max value are added
                else:#if none of the objects has a membership entry then pick just first one
                    best_para.add_row(grouped_mes_table[ind[i]])
        return best_para
    elif para=='binary':
        columns=['main_id',para+'_flag',para+'_qual',para+'_source_idref']
    elif para=='mass_pl':
        columns=['main_id',para+'_value',para+'_rel',para+'_err',para+'_qual',para+'_source_idref']
    elif para=='sep_ang':
        columns=['main_id',para+'_value',para+'_err',para+'_obs_date',para+'_qual',para+'_source_idref']
    else:
        columns=['main_id',para+'_value',para+'_err',para+'_qual',para+'_source_idref']
    mes_table=mes_table[columns[0:]]
    best_para=mes_table[columns[0:]][:0].copy()
    #group mes_table by object (=main_id)
    grouped_mes_table=mes_table.group_by('main_id')
    #take highest quality
    #quite time intensive (few minutes) could maybe be optimized using np.in1d function
    for j in range(len(grouped_mes_table.groups.keys)):#go through all objects
        for qual in ['A','B','C','D','E','?']:
            for i in range(len(grouped_mes_table.groups[j])):
                if grouped_mes_table.groups[j][i][para+'_qual']==qual:
                    best_para.add_row(grouped_mes_table.groups[j][i])
                    break # if best para found go to next main_id group
            else:
                continue  # only executed if the inner loop did NOT break
            break  # only executed if the inner loop DID break 
    return best_para

#------------------------provider combining-----------------
def building(providers,table_names,list_of_tables):
    """
    This function builds the tables for the LIFE database.
    
    :param providers: Containing simbad, grant kennedy, exomercat, gaia and wds data.
    :type providers: list containing objects of type astropy.table.table.Table
    :param table_names: Objects correspond to the names of the 
        astropy tables contained in providers and the return list.
    :type table_names: list of objects of type str
    :param list_of_tables: Empty astropy tables to be filled
        in and returned.
    :type list_of_tables: list containing objects of type astropy.table.table.Table
    :returns: Containing data combined from the different providers.
    :rtype: list containing objects of type astropy.table.table.Table
    """
    
    # Creates empty tables as needed for final database ingestion
    init=initialize_database_tables(table_names,list_of_tables)
    n_tables=len(init)

    cat=[ap.table.Table() for i in range(n_tables)]
    #for the sources and objects joins tables from different providers
    
    print(f'Building {table_names[0]} table ...')#sources
    for j in range(len(providers)):
        if len(cat[0])>0:
            #joining data from different providers
            if len(providers[j][0])>0:
                cat[0]=ap.table.join(cat[0],providers[j][0],join_type='outer')
        else:
            cat[0]=providers[j][0]
        
    #I do this to get those columns that are empty in the data
    cat[0]=ap.table.vstack([cat[0],init[0]])
    
    # keeping only unique values then create identifiers for those tables
    if len(cat[0])>0:
        cat[0]=ap.table.unique(cat[0],silent=True)
        cat[0]['source_id']=[j+1 for j in range(len(cat[0]))]
    
    print(f'Building {table_names[1]} table ...')#objects
    for j in range(len(providers)):
            if len(cat[1])>0:
                #joining data from different providers
                if len(providers[j][1])>0:
                    cat[1]=ap.table.join(cat[1],providers[j][1],keys='main_id',join_type='outer')
                    cat[1]=objectmerging(cat[1])
            else:
                cat[1]=providers[j][1]
                
    #removed vstack with init to not have object_id as is empty anyways
    
    #assigning object_id
    cat[1]['object_id']=[j+1 for j in range(len(cat[1]))]

    # At one point I would like to be able to merge objects with main_id
    # NAME Proxima Centauri b and Proxima Centauri b

    
    print(f'Building {table_names[2]} table ...')#provider
    paras=[['id'],['h_link'],['coo','plx','dist_st','coo_gal','mag_i','mag_j','mag_k','class','sptype'],
           ['mass_pl'],['rad'],['mass_pl'],['teff_st'],
          ['radius_st'],['mass_st'],['binary'],['sep_ang']]
    
    #merging the different provider tables
    for j in range(len(providers)):
        if len(cat[2])>0:
            #joining data from different providers
            if len(providers[j][2])>0:
                cat[2]=ap.table.join(cat[2],providers[j][2],join_type='outer')
        else:
            cat[2]=providers[j][2]
    #I do this to get those columns that are empty in the data
    cat[2]=ap.table.vstack([cat[2],init[2]])
       
    
    for i in range(3,n_tables): # for the tables star_basic,...,mes_mass_st
        print(f'Building {table_names[i]} table ...')
        
        for j in range(len(providers)):#for the different providers (simbad,...,wds)
            if len(providers[j][i])>0:
                #providers[j] is a table containing the two columns ref and provider name
                #replacing ref columns with corresponding source_idref one
                #issue is that order providers and provider_name not the same
                providers[j][i]=match(providers[j][i],cat[0],paras[i-3],providers[j][2]['provider_name'][0])
            if len(cat[i])>0:
                #joining data from different providers
                if len(providers[j][i])>0:
                    cat[i]=ap.table.join(cat[i],providers[j][i],join_type='outer')
            else:
                cat[i]=providers[j][i]
        
        #I do this to get those columns that are empty in the data
        cat[i]=ap.table.vstack([cat[i],init[i]])
        cat[i]=cat[i].filled() #because otherwise unique does neglect masked columns
        

            #I have quite some for loops here, will be slow
            #this function is ready for testing

        if 'object_idref' in cat[i].colnames and len(cat[i])>0: #add object_idref
            #first remove the object_idref we got from empty initialization
            #though I would prefer a more elegant way to do this
            cat[i].remove_column('object_idref') 
            cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],
                                 join_type='left')
            cat[i].rename_column('object_id','object_idref')
        if table_names[i]=='ident':
            cat[i]=best_para('id',cat[i])
        if table_names[i]=='h_link':
            #expanding from child_main_id to object_idref
            #first remove the child_object_idref we got from empty
            # initialization. Would prefer a more elegant way to do this
            cat[i].remove_column('child_object_idref')
            cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],
                                 keys='main_id',join_type='left')
            cat[i].rename_columns(['object_id','main_id'],
                                  ['child_object_idref','child_main_id'])
            
            #expanding from parent_main_id to parent_object_idref
            cat[i].remove_column('parent_object_idref')
            #kick out any h_link rows where parent_main_id not in
            # objects (e.g. clusters)
            cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],
                   keys_left='parent_main_id',keys_right='main_id')
            #removing because same as parent_main_id
            cat[i].remove_column('main_id')
            cat[i].rename_column('object_id','parent_object_idref')
            #null values
            #cat[i]['membership'].fill_value=-1
            #cat[i]['membership']=cat[i]['membership'].filled()
        
        if table_names[i]=='star_basic':
            #choosing all objects with type star or system
            #this I use to join the object_id parameter from objects table to star_basic
            #what about gaia stuff where I don't know this? there I also don't have star_basic info
            #Note: main_id was only added because I have not found out how
            # to do join with just one column of a table
            stars=cat[1]['object_id','main_id'][np.where(
                            cat[1]['type']=='st')]
            systems=cat[1]['object_id','main_id'][np.where(
                            cat[1]['type']=='sy')]
            temp=ap.table.vstack([stars,systems])
            temp.rename_column('object_id','object_idref')
            #cat[i] are all the star_cat tables from providers where those are given
            #the new objects are needed to join the best parameters from mes_ tables later on
            cat[i]=ap.table.join(cat[i],temp,join_type='outer',
                                 keys=['object_idref','main_id'])
        if table_names[i]=='mes_teff_st':
            teff_st_best_para=best_para('teff_st',cat[i])
            cat[5].remove_columns(['teff_st_value','teff_st_err',
                                   'teff_st_qual','teff_st_source_idref',
                                       'teff_st_ref'])
            cat[5]=ap.table.join(cat[5],teff_st_best_para,join_type='left')
        if table_names[i]=='mes_radius_st':
            radius_st_best_para=best_para('radius_st',cat[i])
            cat[5].remove_columns(['radius_st_value','radius_st_err',
                                   'radius_st_qual','radius_st_source_idref',
                                   'radius_st_ref'])
            cat[5]=ap.table.join(cat[5],radius_st_best_para,join_type='left')
        if table_names[i]=='mes_mass_st':
            mass_st_best_para=best_para('mass_st',cat[i])
            cat[5].remove_columns(['mass_st_value','mass_st_err',
                                   'mass_st_qual','mass_st_source_idref',
                                   'mass_st_ref'])
            cat[5]=ap.table.join(cat[5],mass_st_best_para,join_type='left')
        if table_names[i]=='mes_mass_pl':
            mass_pl_best_para=best_para('mass_pl',cat[i])
            cat[6]=mass_pl_best_para
            planets=cat[1]['object_id','main_id'][np.where(
                            cat[1]['type']=='pl')]
            cat[6]=ap.table.join(cat[6],planets['main_id','object_id'])
            cat[6].rename_column('object_id','object_idref')
        if table_names[i]=='mes_binary':
            binary_best_para=best_para('binary',cat[i])
            cat[5].remove_columns(['binary_flag',
                                   'binary_qual','binary_source_idref',
                                   'binary_ref'])
            cat[5]=ap.table.join(cat[5],binary_best_para,join_type='left')
        if table_names[i]=='mes_sep_ang':
            sep_ang_best_para=best_para('sep_ang',cat[i])
            cat[5].remove_columns(['sep_ang_value','sep_ang_err','sep_ang_obs_date',
                                  'sep_ang_qual','sep_ang_source_idref',
                                  'sep_ang_ref'])
            cat[5]=ap.table.join(cat[5],sep_ang_best_para,join_type='left')
        if table_names[i]=='best_h_link':
            cat[i]=best_para('membership',cat[4])
        cat[i]=cat[i].filled()
        
        if len(cat[i])==0:
            print('warning: empty table',i,table_names[i])
        else:
            #only keeping unique entries
            cat[i]=ap.table.unique(cat[i],silent=True)
    cat[5]=cat[5].filled()
    print('Unifying null values...')
    #unify null values (had 'N' and '?' because of ap default fill_value and type conversion string vs object)
    tables=[cat[table_names.index('star_basic')],cat[table_names.index('planet_basic')],
            cat[table_names.index('disk_basic')],cat[table_names.index('mes_mass_pl')],
            cat[table_names.index('mes_teff_st')],cat[table_names.index('mes_radius_st')],
            cat[table_names.index('mes_mass_st')],cat[table_names.index('mes_binary')]]
    columns=[['coo_qual','coo_gal_qual','plx_qual','dist_st_qual',
              'sep_ang_qual','teff_st_qual','radius_st_qual','binary_flag',
              'binary_qual','mass_st_qual','sptype_qual','class_temp','class_temp_nr'],
             ['mass_pl_qual','mass_pl_rel'],
             ['rad_qual','rad_rel'],['mass_pl_qual','mass_pl_rel'],
             ['teff_st_qual'],['radius_st_qual'],['mass_st_qual'],['binary_qual']]
    for i in range(len(tables)):
        for col in columns[i]:
            tables[i]=p.replace_value(tables[i],col,'N','?')
            tables[i]=p.replace_value(tables[i],col,'N/A','?')
    print('TBD: Add exact object distance cut. So far for correct treatment \
          of boundary objects 10% additional distance cut used')
    print('Saving data...')
    hf.save(cat,table_names)
    return cat
