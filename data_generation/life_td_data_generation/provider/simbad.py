""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
import astropy as ap #votables
from datetime import datetime

#self created modules
from utils.utils import save
from provider.utils import fetch_main_id, sources_table, query, nullvalues, replace_value
import sdata as sdc

additional_data_path='../../additional_data/'


def stars_in_multiple_system(cat,sim_h_link,all_objects):
    """
    Assigns object type to special subset of stars.

    This function assignes object type 'st' to those objects that are in
    multiple systems but don't have any stellar child object.

    :param cat: Table alias containing 
    :type cat: astropy.table.table.Table
    :param sim_h_link: Table copy containing columns main_id,
        type and sptype_string.
    :type sim_h_link: astropy.table.table.Table
    :param all_objects: Table copy containing columns 
        parent_main_id and h_link_ref. Rows are all objects with child
        objects.
    :type all_objects: astropy.table.table.Table
    :returns: Table alias like param cat with desired types
        adapted.
    :rtype: astropy.table.table.Table
    """

    #all type sy objects: cat['main_id','type']
    #this should work if alias works well
    #need parent_main_id for sim_h_link here. but setdiff does 
    #not support that.
    parents=sim_h_link['parent_main_id','main_id','h_link_ref'][:]
    parents.rename_column('main_id','child_main_id')
    parents.rename_column('parent_main_id','main_id')
    sy_wo_child=ap.table.setdiff(cat['main_id','type','sptype_string'][:],
                                 parents[:],keys=['main_id'])
    #that don t have children: sy_wo_child['main_id','type']
    #list of those with children
    sy_w_child=ap.table.join(parents[:],
                            cat['main_id','type','sptype_string'][:],
                            keys=['main_id'])
    #list of those with children joined with type of child
    all_objects.rename_columns(['type','main_id'],
                                ['child_type','child_main_id'])
    sy_w_child=ap.table.join(sy_w_child[:],
                            all_objects['child_type','child_main_id'][:],
                             keys=['child_main_id'],join_type='left')
    #remove all where type child is not pl
    sy_w_child_pl=sy_w_child[np.where(sy_w_child['child_type']=='pl')]
    if len(sy_w_child_pl)==0:
        #no systems with child of type planet
        sy_wo_child_st=sy_wo_child
    else:
        #join with list of sy that dont habe children
        sy_wo_child_st=ap.table.vstack([sy_wo_child[:],sy_w_child_pl[:]])
        sy_wo_child_st.remove_column('child_type')
    #systems that don t have children except planets: sy_wo_child_st
    #no + in sptype_string because that is another indication of binarity
    temp=[len(i.split('+'))==1 for i in sy_wo_child_st['sptype_string']]
    #have it as an array of bools 
    temp=np.array(temp)
    #have it as lisit of indices 
    temp=list(np.where(temp==True)[0])
    single_sptype=sy_wo_child_st[:][temp]
    #and no + in spectral type: single_sptype['main_id','type']      
    cat['type'][np.where(np.in1d(cat['main_id'],
                                single_sptype['main_id']))]=\
              ['st' for j in range(len(cat[np.where(np.in1d(cat['main_id'],
                                        single_sptype['main_id']))]))]        
    return cat

#-----------------------------provider data ingestion-------------------
def provider_simbad(sim_list_of_tables,distance_cut_in_pc,
                    test_objects=[]):
    """
    Optains and arranges SIMBAD data.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param sim_list_of_tables: Contains empty output tables.
    :type sim_list_of_tables: list(astropy.table.table.Table)
    :returns: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data, basic stellar data and binarity data.
    :rtype: list(astropy.table.table.Table)
    """
    
    plx_in_mas_cut=1000./distance_cut_in_pc
    #making cut a bit bigger for correct treatment of objects on boundary
    plx_cut=plx_in_mas_cut-plx_in_mas_cut/10.
    #---------------define provider-------------------------------------
    sdc_simbad=sdc.provider('simbad')
    table_names=sdc_simbad.table_names
    sdc_simbad.table('provider').add_row()
    sdc_simbad.table('provider')['provider_name']='SIMBAD',
    sdc_simbad.table('provider')['provider_url']=\
            "http://simbad.u-strasbg.fr:80/simbad/sim-tap",
    sdc_simbad.table('provider')['provider_bibcode']='2000A&AS..143....9W'
    sdc_simbad.table('provider')['provider_access']= \
            datetime.now().strftime('%Y-%m-%d')
    
    sim_provider=sdc_simbad.table('provider')
    #---------------define queries--------------------------------------
    select="""SELECT b.main_id,b.ra AS coo_ra,b.dec AS coo_dec,
        b.coo_err_angle, b.coo_err_maj, b.coo_err_min,b.oid,
        b.coo_bibcode AS coo_ref, b.coo_qual,b.sp_type AS sptype_string,
        b.sp_qual AS sptype_qual, b.sp_bibcode AS sptype_ref,
        b.plx_err, b.plx_value, b.plx_bibcode AS plx_ref,b.plx_qual,
        h_link.membership, h_link.parent AS parent_oid,
        h_link.link_bibcode AS h_link_ref, a.otypes,ids.ids,
        f.I as mag_i_value, f.J as mag_j_value, f.K as mag_k_value
        """#which parameters to query from simbad and what alias to give them
    #,f.I as mag_i_value, f.J as mag_j_value
    tables="""
    FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    LEFT JOIN allfluxes AS f ON b.oid=f.oidref
    
    """
    #JOIN allfluxes AS f ON b.oid=f.oidref
    
    adql_query=[
        select+
        tables+
        'WHERE b.plx_value >='+str(plx_cut)]
    #creating one table out of parameters from multiple ones and
    #keeping only objects with parallax bigger than ... mas

    upload_query=[
        #query for systems without parallax data but
        #children (in TAP_UPLOAD.t1 table) with parallax bigger than 50mas
        select+
        tables+
        """JOIN TAP_UPLOAD.t1 ON b.oid=t1.parent_oid
        WHERE (b.plx_value IS NULL) AND (otype='**..')""",
        #query for planets without parallax data but
        #host star (in TAP_UPLOAD.t1 table) with parallax bigger than 50mas
        select+
        tables+
        """JOIN TAP_UPLOAD.t1 ON b.oid=t1.oid
        WHERE (b.plx_value IS NULL) AND (otype='Pl..')""",
        #query all distance measurements for objects in TAP_UPLOAD.t1 table
        """SELECT oid, dist AS dist_st_value, plus_err, qual AS dist_st_qual,
        bibcode AS dist_st_ref,minus_err,dist_prec AS dist_st_prec
        FROM mesDistance
        JOIN TAP_UPLOAD.t1 ON oidref=t1.oid""",
        #query all identifiers for objects in TAP_UPLOAD.t1 table
        """SELECT id, t1.*
        FROM ident
        JOIN TAP_UPLOAD.t1 ON oidref=t1.oid"""]
    #------------------querrying----------------------------------------
    print('Creating ',sim_provider['provider_name'][0],' tables ...')
    #perform query for objects with in distance given
    simbad=query(sim_provider['provider_url'][0],adql_query[0])
    #querries parent and children objects with no parallax value
    parents_without_plx=query(sim_provider['provider_url'][0],
                                upload_query[0],[simbad])
    children_without_plx=query(sim_provider['provider_url'][0],
                                upload_query[1],[simbad])
    
    test_objects=np.array(test_objects)
    if len(test_objects)>0:
        print('in sim through plx query', 
                  test_objects[np.where(np.in1d(test_objects,
                                                simbad['main_id']))])
        print('in sim through child plx query', 
                test_objects[np.where(np.in1d(test_objects,
                                            parents_without_plx['main_id']))])
        print('in sim through parent plx query', 
                  test_objects[np.where(np.in1d(test_objects,
                                            children_without_plx['main_id']))])
    
    #adding of no_parallax objects to rest of simbad query objects
    
    simbad=ap.table.vstack([simbad,parents_without_plx])
    simbad=ap.table.vstack([simbad,children_without_plx])
    
    print(' sorting object types...')

    #sorting from object type into star, system and planet type
    simbad['type']=['None' for i in range(len(simbad))]
    simbad['binary_flag']=['False' for i in range(len(simbad))]
    to_remove_list=[]
    removed_otypes=[]
    for i in range(len(simbad)):
        #planets
        if "Pl" in simbad['otypes'][i]:
            simbad['type'][i]='pl'
        #stars
        elif "*" in simbad['otypes'][i]:
            #system containing multiple stars
            if "**" in simbad['otypes'][i]:
                simbad['type'][i]='sy'
                simbad['binary_flag'][i]='True'
            #individual stars
            else:
                simbad['type'][i]='st'
        else:
            removed_otypes.append(simbad['otypes'][i])
            #most likely single brown dwarfs
            #storing information for later removal from table called simbad
            to_remove_list.append(i)
    #removing any objects that are neither planet, star nor system in type
    if to_remove_list!=[]:
        print('removing',len(removed_otypes),' objects that had object types:',
              list(set(removed_otypes)))
        print('example object of them:', simbad['main_id'][to_remove_list[0]])
        simbad.remove_rows(to_remove_list)
        
    if len(test_objects)>0:
        print('in sim through otype criteria', 
                  test_objects[np.where(np.in1d(test_objects,
                                                simbad['main_id']))])

    #creating helpter table stars
    temp_stars=simbad[np.where(simbad['type']!='pl')]
    #removing double objects (in there due to multiple parents)
    stars=ap.table.Table(ap.table.unique(temp_stars,keys='main_id'),copy=True)
    
    print(' creating output tables...')
    #-----------------creating output table sim_ident-------------------
    #issue if I want to replace this here with sdc is that I have main_id column but no id_source_idref
    #sdc_simbad.table('ident').remove_column('id_source_idref')
    sim_ident=query(sim_provider['provider_url'][0],upload_query[3],
                    [simbad['oid','main_id'][:].copy()]) #adds column id
    sim_ident['id_ref']=[sim_provider['provider_bibcode'][0] \
                            for j in range(len(sim_ident))]
    sim_ident.remove_column('oid')
    
    #--------------creating output table sim_h_link --------------------
    sim_h_link=simbad['main_id','parent_oid','h_link_ref','membership']
    #sim_h_link=nullvalues(sim_h_link,'parent_oid',0,verbose=False)
    ###sim_h_link=nullvalues(sim_h_link,'membership',-1,verbose=False)
    
    # removing entries in h_link where parent objects are clusters or 
    # associations as we are 
    #only interested in hierarchical multiples. 
    sim_h_link=sim_h_link[np.where(np.in1d(sim_h_link['parent_oid'],
                                            stars['oid']))]
    
    
    sim_h_link=fetch_main_id(sim_h_link,'parent_oid','parent_main_id')
    sim_h_link.remove_column('parent_oid')
    #typeconversion needed as smallint fill value != int null value
    sim_h_link['membership']=sim_h_link['membership'].astype(int)
    sim_h_link=nullvalues(sim_h_link,'membership',999999)
    sim_h_link=replace_value(sim_h_link,'h_link_ref','',
                             sim_provider['provider_bibcode'][0])
    sim_h_link=ap.table.unique(sim_h_link)
                
    #--------------------creating helper table sim_stars----------------
    #updating multiplicity object type
    #no children and sptype does not contain + -> type needs to be st

    #all objects in stars table: stars['main_id','type']
    stars[np.where(stars['type']=='sy')]=stars_in_multiple_system(
            stars[np.where(stars['type']=='sy')],sim_h_link[:],
            simbad['main_id','type'][:])    
    
    # binary_flag 'True' for all stars with parents
    # meaning stars[main_id] in sim_h_link[child_main_id] 
    #-> stars[binary_flag]=='True'    
    stars['binary_flag'][np.where(np.in1d(stars['main_id'],
                                        sim_h_link['main_id']))]=\
                    ['True' for j in range(len(stars[np.where(
                    np.in1d(stars['main_id'],sim_h_link['main_id']))]))]   
                
    #change null value of plx_qual
    stars['plx_qual']=stars['plx_qual'].astype(object)
    stars=replace_value(stars,'plx_qual','',stars['plx_qual'].fill_value)
    
    for band in ['i','j','k']:
        #initiate some of the ref columns
        stars[f'mag_{band}_ref']=ap.table.MaskedColumn(dtype=object,
                                    length=len(stars),
                                    mask=[True for j in range(len(stars))])
        #add simbad reference where no other is given
        stars[f'mag_{band}_ref'][np.where(
                stars[f'mag_{band}_value'].mask==False)]=[
                sim_provider['provider_bibcode'][0] for j in range(len(
                stars[f'mag_{band}_ref'][np.where(
                stars[f'mag_{band}_value'].mask==False)]))]
        
    stars=replace_value(stars,'plx_ref','',sim_provider['provider_bibcode'][0])
    stars=replace_value(stars,'sptype_ref','',
            sim_provider['provider_bibcode'][0])
    stars=replace_value(stars,'coo_ref','',sim_provider['provider_bibcode'][0])
        
    stars['binary_ref']=[sim_provider['provider_bibcode'][0] for j in range(
            len(stars))]
    stars['binary_qual']=['D' for j in range(len(stars))]

    #-----------------creating output table sim_planets-----------------
    temp_sim_planets=simbad['main_id','ids',
                            'type'][np.where(simbad['type']=='pl')]
    sim_planets=ap.table.Table(ap.table.unique(
                    temp_sim_planets,keys='main_id'),copy=True)
    #-----------------creating output table sim_objects-----------------
    sim_objects=ap.table.vstack([sim_planets['main_id','ids','type'],
                             stars['main_id','ids','type']])
    sim_objects['ids']=sim_objects['ids'].astype(object)
    #tbd: add identifier simbad main_id without leading * and whitespaces
    #--------------creating output table sim_sources -------------------
    sim_sources=ap.table.Table()
    tables=[sim_provider,stars, sim_h_link,sim_ident]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['coo_ref','plx_ref','mag_i_ref',
                    'mag_j_ref','mag_k_ref','binary_ref','sptype_ref'],
                    ['h_link_ref'],['id_ref']]
    for cat,ref in zip(tables,ref_columns):
        sim_sources=sources_table(cat,ref,sim_provider['provider_name'][0],
                                sim_sources)
    #------------------------creating output table sim_star_basic-------
    sim_star_basic=stars['main_id','coo_ra','coo_dec','coo_err_angle',
                         'coo_err_maj','coo_err_min','coo_qual','coo_ref',
                         'mag_i_value','mag_i_ref','mag_j_value','mag_j_ref',
                         'mag_k_value','mag_k_ref',
                         'sptype_string','sptype_qual','sptype_ref',
                         'plx_value','plx_err','plx_qual','plx_ref']
    #-----------creating mes_binary table-------------------------------
    sim_mes_binary=stars['main_id','binary_flag','binary_qual','binary_ref']
    #-------------changing type from object to string for later join functions
    sim_star_basic['sptype_string']=sim_star_basic['sptype_string'].astype(str)
    sim_star_basic['sptype_qual']=sim_star_basic['sptype_qual'].astype(str)
    sim_star_basic['sptype_ref']=sim_star_basic['sptype_ref'].astype(str)
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': sim_list_of_tables[i]=sim_sources
        if table_names[i]=='provider': sim_list_of_tables[i]=sim_provider
        if table_names[i]=='objects': sim_list_of_tables[i]=sim_objects
        if table_names[i]=='ident': sim_list_of_tables[i]=sim_ident
        if table_names[i]=='h_link': sim_list_of_tables[i]=sim_h_link
        if table_names[i]=='star_basic': sim_list_of_tables[i]=sim_star_basic
        if table_names[i]=='mes_binary': sim_list_of_tables[i]=sim_mes_binary
        save([sim_list_of_tables[i][:]],['sim_'+table_names[i]])
    return sim_list_of_tables
