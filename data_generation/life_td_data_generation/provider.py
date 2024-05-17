import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
from datetime import datetime
import importlib #reloading external functions after modification

#self created modules
import helperfunctions as hf
importlib.reload(hf)#reload module after changing it
import sdc
importlib.reload(sdc)#reload module after changing it


#------------------------------provider helper functions----------------
def query(link,query,catalogs=[]):
    """
    Performs a query via TAP on the service given in the link parameter.
    
    If a list of tables is given in the catalogs parameter,
    those are uploaded to the service beforehand.
    
    :param str link: Service access URL.
    :param str query: Query to be asked of the external database service
         in ADQL.
    :param catalogs: List of astropy tables to be uploaded to the 
        service.
    :type catalogs: list(astropy.table.table.Table)
    :returns: Result of the query.
    :rtype: astropy.table.table.Table
    """
    
    #defining the vo service using the given link
    service = vo.dal.TAPService(link)
    #without upload tables
    if catalogs==[]:
        result=service.run_async(query.format(**locals()), maxrec=1600000)
    #with upload tables
    else:
        tables={}
        for i in range(len(catalogs)):
            tables.update({f"t{i+1}":catalogs[i]})
        result = service.run_async(query,uploads=tables,timeout=None,
                                   maxrec=1600000)
    cat=result.to_table()
    return cat

def sources_table(cat,ref_columns,provider,old_sources=ap.table.Table()):
    """
    Creates or updates the source table out of the given references.
    
    The entries are unique and the columns consist out of the
    reference and provider_name.
    
    :param cat: Table on which the references should be gathered.
    :type cat: astropy.table.table.Table
    :param ref_columns: Header of the columns containing reference 
        information.
    :type ref_columns:
    :param str provider: Provider name.
    :param old_sources: Previously created reference table.
    :type old_sources:
    :return: Table containing references and provider information.
    :rtype: astropy.table.table.Table
    """
    
    if len(cat)>0:
        # table initialization to prevent error messages when assigning 
        # columns
        cat_sources=ap.table.Table() 
        #initialization of list to store reference information
        cat_reflist=[] 
        #for all the columns given add reference information 
        for k in range(len(ref_columns)):
            #In case the column has elements that are masked skip those
            if type(cat[ref_columns[k]])==ap.table.column.MaskedColumn:
                cat_reflist.extend(
                    cat[ref_columns[k]][np.where(
                            cat[ref_columns[k]].mask==False)])
            else:
                cat_reflist.extend(cat[ref_columns[k]])
        # add list of collected references to the table and call the 
        # column ref
        cat_sources['ref']=cat_reflist
        cat_sources=ap.table.unique(cat_sources)
        #attaches service information
        cat_sources['provider_name']=[provider for j in range(
                len(cat_sources))]
        #combine old and new sources into one table
        sources=ap.table.vstack([old_sources,cat_sources])
        sources=ap.table.unique(sources) #remove double entries
    else:
        sources=old_sources
    return sources

def fetch_main_id(cat,colname='oid',name='main_id',oid=True):
    """
    Joins main_id from simbad to the column colname. 
    
    Returns the whole table cat but without any rows where no simbad 
    main_id was found.
    
    :param cat: Astropy table containing column colname.
    :type cat: astropy.table.table.Table
    :param str colname: Column header of the identifiers that should be 
        searched in SIMBAD.
    :param str name: Column header for the SIMBAD main identifiers, 
        default is main_id.
    :param bool oid: Specifies wether colname is a SIMBAD oid or normal
         identifier.
    :return: Table with all main SIMBAD identifiers that could be found 
        in column "name".
    :rtype: astropy.table.table.Table
    """
    
    #improvement idea to be performed at one point
    # tbd option to match on position instead of main_id or oid
    #SIMBAD TAP service
    TAP_service="http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    #creating oid query
    if oid:
        main_id_query='SELECT b.main_id AS '+name+""",t1.*
                    FROM basic AS b
                    JOIN TAP_UPLOAD.t1 ON b.oid=t1."""+colname
    #creating identifier query
    else:
        main_id_query='SELECT b.main_id AS '+name+""",t1.*
                    FROM basic AS b
                    JOIN ident ON ident.oidref=b.oid
                        JOIN TAP_UPLOAD.t1 ON ident.id=t1."""+colname
    #performing query using external function
    cat=query(TAP_service,main_id_query,[cat])
    return cat

def distance_cut(cat,colname,main_id=True):
    """
    Sorts out objects not within the provider_simbad distance cut. 
    
    :param cat: Astropy table to be matched against sim_objects table.
    :type cat: astropy.table.table.Table
    :param str colname: Name of the column to use for the match.
    :return: Table like cat without any objects not found in 
        sim_objects.
    :rtype: astropy.table.table.Table
    """
    
    if main_id:
        [sim]=hf.load(['sim_objects'])
        sim.rename_columns(['main_id','ids'],['temp1','temp2'])
        cat=ap.table.join(cat,sim['temp1','temp2'],
                      keys_left=colname,keys_right='temp1')
        cat.remove_columns(['temp1','temp2'])
    else:
        [sim]=hf.load(['sim_ident'])
        sim.rename_columns(['id'],['temp1'])
        cat=ap.table.join(cat,sim['temp1','main_id'],
                      keys_left=colname,keys_right='temp1')
        cat.remove_columns(['temp1'])
    return cat

def nullvalues(cat,colname,nullvalue,verbose=False):
    """
    This function fills masked entries specified column. 
    
    :param cat: Astropy table containing the column colname.
    :type cat: astropy.table.table.Table
    :param str colname: Name of a column.
    :param nullvalue: Value to be placed instead of masked elements.
    :type nullvalue: str or float or bool
    :param verbose: If True prints message if the column is not an 
        astropy masked column, defaults to False.
    :type verbose: bool, optional
    :return: Astropy table with masked elements of colname replaced
        by nullvalue.
    :rtype: astropy.table.table.Table
    """
    
    if type(cat[colname])==ap.table.column.MaskedColumn:
                cat[colname].fill_value=nullvalue
                cat[colname]=cat[colname].filled()
    elif verbose:
        print(colname,'is no masked column')
    return cat

def replace_value(cat,column,value,replace_by):
    """
    This function replaces values.
    
    :param cat: Table containing column specified as colname.
    :type cat: astropy.table.table.Table
    :param str column: Designates column in which to replace the 
        entries.
    :param value: Entry to be replaced.
    :type value: str or float or bool
    :param replace_by: Entry to be put in place of param value.
    :type replace_by: str or float or bool
    :return: Table with replaced entries.
    :rtype: astropy.table.table.Table
    """
    
    cat[column][np.where(cat[column]==value)]= \
            [replace_by for i in range(
        len(cat[column][np.where(cat[column]==value)]))]
    return cat

def ids_from_ident(ident,objects):
    """
    Concatenates identifier entries of same main_id object.
    
    This function extracts the identifiers of common main_id objects in 
    column id of table ident using the delimiter | and stores the result
    in the column ids of the table objects.
    
    :param ident: Table containing the rows main_id and id
    :type ident: astropy.table.table.Table
    :param objects: Table containing the columns main_id and ids
    :type objects: astropy.table.table.Table
    :returns: Filled out table objects.
    :rtype: astropy.table.table.Table
    """
    
    grouped_ident=ident.group_by('main_id')
    ind=grouped_ident.groups.indices
    for i in range(len(ind)-1):
    # -1 is needed because else ind[i+1] is out of bonds
        ids=[]
        for j in range(ind[i],ind[i+1]):
            ids.append(grouped_ident['id'][j])
        ids="|".join(ids)
        objects.add_row([grouped_ident['main_id'][ind[i]],ids])
    return objects

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
    sdc_simbad.table('provider')['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
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
                  test_objects[np.where(np.in1d(test_objects,parents_without_plx['main_id']))])
        print('in sim through parent plx query', 
                  test_objects[np.where(np.in1d(test_objects,children_without_plx['main_id']))])
    
    #adding of no_parallax objects to rest of simbad query objects
    simbad=ap.table.vstack([simbad,parents_without_plx])
    simbad=ap.table.vstack([simbad,children_without_plx])
    #----------------------sorting object types-------------------------
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
    
    #-----------------creating output table sim_ident-------------------
    #issue if I want to replace this here with sdc is that I have main_id column but no id_source_idref
    #sdc_simbad.table('ident').remove_column('id_source_idref')
    sim_ident=query(sim_provider['provider_url'][0],upload_query[3],
                    [simbad['oid','main_id'][:].copy()]) #adds column id
    sim_ident['id_ref']=[sim_provider['provider_bibcode'][0] for j in range(len(sim_ident))]
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
    print('tbd: add identifier simbad main_id without leading * and \
            whitespaces')
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
        hf.save([sim_list_of_tables[i][:]],['sim_'+table_names[i]])
    return sim_list_of_tables


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
    gk_provider=ap.table.Table()
    gk_provider['provider_name']=['Grant Kennedy Disks']
    gk_provider['provider_url']=['http://drgmk.com/sdb/']
    gk_provider['provider_bibcode']=['priv. comm.']
    gk_provider['provider_access']=['2022-12-06']
    
    print('Creating ',gk_provider['provider_name'][0],' tables ...')
    #loading table obtained via direct communication from Grant Kennedy
    gk_disks=ap.io.votable.parse_single_table(
        "../../data/additional_data/Grant_absil_2013_.xml").to_table()
    #transforming from string type into object to have variable length
    gk_disks=hf.stringtoobject(gk_disks,212)
    #removing objects with plx_value=='None' or masked entries
    gk_disks['plx_value']=gk_disks['plx_value'].filled('None')
    gk_disks=gk_disks[np.where(gk_disks['plx_value']!='None')]
    gk_disks['plx_value']=gk_disks['plx_value'].astype(float)
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
    gk_sources=ap.table.Table()
    tables=[gk_provider,gk_disks]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['disks_ref']]
    for cat,ref in zip(tables,ref_columns):
        gk_sources=sources_table(cat,ref,gk_provider['provider_name'][0],gk_sources)
    #--------------creating output table gk_disk_basic------------------
    gk_disk_basic=gk_disks['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    for column in ['rdisk_bb','e_rdisk_bb']:
        #replacing 'None' with 'nan' as the first one is not float convertible
        gk_disk_basic=replace_value(gk_disk_basic,column,'None','nan')
        gk_disk_basic[column].fill_value='nan' #because defeault is None and not float convertible
        #though this poses the issue that the float default float fill_value is 1e20
        gk_disk_basic[column].filled()
        gk_disk_basic[column]=gk_disk_basic[column].astype(float)
    gk_disk_basic.rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    gk_disk_basic=gk_disk_basic[np.where(np.isfinite(gk_disk_basic['rad_value']))]
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': gk_list_of_tables[i]=gk_sources
        if table_names[i]=='provider': gk_list_of_tables[i]=gk_provider
        if table_names[i]=='objects': gk_list_of_tables[i]=gk_objects
        if table_names[i]=='ident': gk_list_of_tables[i]=gk_ident
        if table_names[i]=='h_link': gk_list_of_tables[i]=gk_h_link
        if table_names[i]=='disk_basic': gk_list_of_tables[i]=gk_disk_basic
        hf.save([gk_list_of_tables[i][:]],['gk_'+table_names[i]])
    return gk_list_of_tables



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
    exo_provider=ap.table.Table()
    exo_provider['provider_name']=['Exo-MerCat']
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    exo_provider['provider_bibcode']=['2020A&C....3100370A']
    
    
    print('Creating ',exo_provider['provider_name'][0],' tables ...')
    #---------------define query----------------------------------------
    adql_query="""SELECT *
                  FROM exomercat.exomercat"""
    #---------------obtain data-----------------------------------------
    if temp:
        exomercat=ap.io.ascii.read("../../data/additional_data/exo-mercat05-02-2023_v2.0.csv")
        exomercat=stringtoobject(exomercat,3000)
        exo_provider['provider_access']=['2023-02-05']

    else:
        exomercat=query(exo_provider['provider_url'][0],adql_query)
        exo_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    #----------------putting object main identifiers together-----------
    
    # initializing column
    exomercat['planet_main_id']=ap.table.Column(dtype=object,
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
                             colname='planet_main_id',name='sim_planet_main_id',oid=False)

    
    notinsimbad=exomercat['planet_main_id'][np.where(np.in1d(exomercat['planet_main_id'],
                                                             exomercat2['planet_main_id'],invert=True))]
    #I use a left join as otherwise I would loose some objects that are not in simbad
    exomercat=ap.table.join(exomercat,exomercat2['sim_planet_main_id','planet_main_id'],
                            keys='planet_main_id',join_type='left')

    #show which elements from exomercat were not found in sim_objects
    exo['name']=exo['name'].astype(object)
    removed_objects=ap.table.setdiff(exo,exomercat,keys=['name'])
    hf.save([removed_objects],['exomercat_removed_objects'])

    #-------------exo_ident---------------
    exo_ident=ap.table.Table(names=['main_id','id'],dtype=[object,object])
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
    exo_objects=ap.table.Table(names=['main_id','ids'],dtype=[object,object])
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
    exomercat['mass_pl_err']=ap.table.Column(dtype=float,length=len(exomercat))
    exomercat['mass_pl_rel']=ap.table.Column(dtype=object,length=len(exomercat))
    exomercat['mass_pl_qual']=ap.table.MaskedColumn(dtype=object,length=len(exomercat))
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
                print('tbd: check if relation is correct in case of maximum error on a lower limit value')
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
    print('tbd: include masssini measurements from exomercat')

    #-------------exo_h_link---------------
    exo_h_link=exomercat['planet_main_id', 'host_main_id']
    exo_h_link.rename_columns(['planet_main_id','host_main_id'],
                              ['main_id','parent_main_id'])
    exo_h_link['h_link_ref']=[exo_provider['provider_bibcode'][0] for j in range(len(exo_h_link))]
    #-------------exo_sources---------------
    ref_columns=[['provider_bibcode'],['h_link_ref'],['id_ref'],['mass_pl_ref']]
    exo_sources=ap.table.Table()
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
        hf.save([exo_list_of_tables[i][:]],['exo_'+table_names[i]])
    return exo_list_of_tables

def sptype_string_to_class(temp,ref):
    """
    Extracts stellar parameters from spectral type string one.

    This function extracts the temperature class, temperature class number
    and luminocity class information from the spectral type string (e.g. 
    M5V to M, 5 and V). It stores that information in the for this purpose
    generated new columns. Only objects of temperature class O, B, A, F,
    G, K, and M are processed. Only objects of luminocity class IV, V and VI
    are processed.

    :param temp: Table containing spectral type information in
        the column sptype_string.
    :type temp: astropy.table.table.Table
    :param str ref: Designates origin of data.
    :returns: Table like temp with additional columns class_temp,
        class_temp_nr, class_lum and class_ref.
    :rtype: astropy.table.table.Table
    """

    temp['class_temp']=ap.table.MaskedColumn(dtype=object,length=len(temp))
    temp['class_temp_nr']=ap.table.MaskedColumn(dtype=object,length=len(temp))
    temp['class_lum']=ap.table.MaskedColumn(dtype=object,length=len(temp))
    temp['class_ref']=ap.table.MaskedColumn(dtype=object,length=len(temp))
    #tbd: rewrite code using recoursive function

    for i in range(len(temp)):
        #sorting out objects like M5V+K7V
        #strip d for spectral types starting with small d because it is an old annotation for dwarf star
        sptype=temp['sptype_string'][i].strip('d')
        if (len(sptype.split('+'))==1 and
        #sorting out entries like '', DA2.9, T1V
                len(sptype)>0 and 
                sptype[0] in ['O','B','A','F','G','K','M']):
            temp['class_temp'][i]=sptype[0]
            temp['class_ref'][i]=ref
            #sorting out objects like DA2.9
            if len(sptype)>1 and sptype[1] in ['0','1','2','3','4','5','6','7','8','9']:
                temp['class_temp_nr'][i]=sptype[1]
                #distinguishing between objects like K5V and K5.5V
                if len(sptype)>2 and sptype[2]=='.':
                    temp['class_temp_nr'][i]=sptype[1:4]
                    if len(sptype)>4 and sptype[4] in ['I','V']:
                        temp['class_lum'][i]=sptype[4]
                        if len(sptype)>5 and sptype[5] in ['I','V']:
                            temp['class_lum'][i]=sptype[4:6]
                            if len(sptype)>6 and sptype[6] in ['I','V']:
                                temp['class_lum'][i]=sptype[4:7]
                    else:
                        temp['class_lum'][i]='?'
                elif len(sptype)>2 and sptype[2] in ['I','V']:
                    temp['class_lum'][i]=sptype[2]
                    if len(sptype)>3 and sptype[3] in ['I','V']:
                        temp['class_lum'][i]=sptype[2:4]
                        if len(sptype)>4 and sptype[4] in ['I','V']:
                            temp['class_lum'][i]=sptype[2:5]
                else:
                    temp['class_lum'][i]='?'
            else:
                temp['class_lum'][i]='?'
        else:
            temp['class_temp'][i]='?'
            temp['class_temp_nr'][i]='?'
            temp['class_lum'][i]='?'
            temp['class_ref'][i]='?'
        if len(temp['sptype_string'][i])>0:
            if temp['sptype_string'][i][0]=='d':
                temp['class_lum'][i]='V'
    return temp

def realspectype(cat):
    """
    Removes rows not containing main sequence stars.

    Removes rows of cat where elements in column named 'sim_sptype' are
    either '', 'nan' or start with an other letter than the main sequence
    spectral type classification.

    :param cat: Table containing 'sim_sptype' column
    :type cat: astropy.table.table.Table
    :returns: Table, param cat with undesired rows removed
    :rtype: astropy.table.table.Table
    """

    index=[]
    for j in range(len(cat['sptype_string'])):
        if cat['sptype_string'][j] in ['','nan']:
            index.append(j)
        elif cat['sptype_string'][j][0] not in ['O','B','A','F','G','K','M']:
            index.append(j)
    cat.remove_rows(index)
    return cat

def model_param():
    """
    Loads and cleans up model file.

    Loads the table of Eric E. Mamajek containing stellar parameters 
    modeled from spectral types. Cleans up the columns for spectral 
    type, effective temperature radius and mass.

    :returns: Table of the 4 parameters as columns
    :rtype: astropy.table.table.Table
    """

    EEM_table=ap.io.ascii.read("../../data/additional_data/Mamajek2022-04-16.csv")['SpT','Teff','R_Rsun','Msun']
    EEM_table.rename_columns(['R_Rsun','Msun'],['Radius','Mass'])
    EEM_table=replace_value(EEM_table,'Radius',' ...','nan')
    EEM_table=replace_value(EEM_table,'Mass',' ...','nan')
    EEM_table=replace_value(EEM_table,'Mass',' ....','nan')
    EEM_table['Teff'].unit=ap.units.K
    EEM_table['Radius'].unit=ap.units.R_sun
    EEM_table['Mass'].unit=ap.units.M_sun       
    ap.io.votable.writeto(ap.io.votable.from_table(EEM_table), \
                          f'../../data/model_param.xml')#saving votable
    return EEM_table

def match_sptype(cat,model_param,sptypestring='sim_sptype',teffstring='mod_Teff',\
                 rstring='mod_R',mstring='mod_M'):
    """
    Assigns modeled parameter values.

    Matches the spectral types with the ones in Mamajek's table and 
    includes the modeled effective Temperature,
    stellar radius and stellar mass into the catalog.

    :param cat: astropy table containing spectral type information
    :type cat: astropy.table.table.Table
    :param str sptypestring: Column name where the spectral 
        type information is located
    :param str teffstring: Column name for the stellar effective 
        temperature column
    :param str rstring: Column name for the stellar radius column
    :param str mstring: Column name for the stellar mass column
    :returns: Table cat with added new columns for 
        effective temperature, radius and mass filled with model values
    :rtype: astropy.table.table.Table
    """

    #initiating columns with right units

    arr=np.zeros(len(cat))
    cat[teffstring]=arr*np.nan*ap.units.K
    cat[teffstring]=ap.table.MaskedColumn(mask=np.full(len(cat),True), \
                                         length=len(cat),unit=ap.units.K)
    cat[rstring]=arr*np.nan*ap.units.R_sun
    cat[mstring]=arr*np.nan*ap.units.M_sun
    #go through all spectral types in cat
    for j in range(len(cat[sptypestring])): 
        # for all the entries that are not empty
        if cat[sptypestring][j]!='':
            #go through the model spectral types of Mamajek 
            for i in range(len(model_param['SpT'])): 
                #match first two letters
                if model_param['SpT'][i][:2]==cat[sptypestring][j][:2]: 
                        cat[teffstring][j]=model_param['Teff'][i]
                        cat[rstring][j]=model_param['Radius'][i]
                        cat[mstring][j]=model_param['Mass'][i]
            #as the model does not cover all spectral types on .5 accuracy, check those separately
            if cat[sptypestring][j][2:4]=='.5':
                for i in range(len(model_param['SpT'])):
                    # match first four letters
                    if model_param['SpT'][i][:4]==cat[sptypestring][j][:4]:
                        cat[teffstring][j]=model_param['Teff'][i]
                        cat[rstring][j]=model_param['Radius'][i]
                        cat[mstring][j]=model_param['Mass'][i] 
        else:
            cat[sptypestring][j]='None' 
    return cat

def spec(cat):
    """
    Runs the spectral type related functions realspectype and match_sptype. 

    It also removes all empty columns of the effective temperature, removes 
    rows that are not main sequence, removes rows with binary subtype and 
    non unique simbad name.

    :param cat: astropy table containing columns named 
        'sim_sptype','sim_name' and 'sim_otypes'
    :type cat: astropy.table.table.Table
    :returns: Catalog of mainsequence stars with unique 
        simbad names, no binary subtypes and modeled parameters.
    :rtype: astropy.table.table.Table
    """   

    #Do I even need realspectype function? I can just take cat where class_temp not empty
    cat=realspectype(cat)
    #model_param=ap.io.votable.parse_single_table(\
        #f"catalogs/model_param.xml").to_table()
    mp=model_param()#create model table as votable
    cat=match_sptype(cat,mp,sptypestring='sptype_string')
    cat.remove_rows([np.where(cat['mod_Teff'].mask==True)])
    cat.remove_rows([np.where(np.isnan(cat['mod_Teff']))])
    cat=ap.table.unique(cat, keys='main_id')
    return cat

def provider_life(table_names,life_list_of_tables):
    """
    Loads SIMBAD data and postprocesses it. 
    
    Postprocessing enables to provide more useful information. It uses a model
    from Eric E. Mamajek to predict temperature, mass and radius from the simbad 
    spectral type data.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param life_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :type life_list_of_tables: list(astropy.table.table.Table)
    :returns: List of astropy table containing
        reference data, provider data, basic stellar data, stellar effective
        temperature data, stellar radius data and stellar mass data.
    :rtype: list(astropy.table.table.Table)
    """
    
    #---------------define provider-------------------------------------
    life_provider=ap.table.Table()
    life_provider['provider_name']=['LIFE']
    life_provider['provider_url']=['www.life-space-mission.com']
    life_provider['provider_bibcode']=['2022A&A...664A..21Q']
    life_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
    print('Creating ',life_provider['provider_name'][0],' tables ...')
    #---------------------star_basic----------------
    #galactic coordinates:  transformed from simbad ircs coordinates using astropy
    [life_star_basic]=hf.load(['sim_star_basic'])
    ircs_coord=ap.coordinates.SkyCoord(\
            ra=life_star_basic['coo_ra'],dec=life_star_basic['coo_dec'],frame='icrs')
    gal_coord=ircs_coord.galactic
    life_star_basic['coo_gal_l']=gal_coord.l.deg*ap.units.degree
    life_star_basic['coo_gal_b']=gal_coord.b.deg*ap.units.degree
    life_star_basic['dist_st_value']=1000./life_star_basic['plx_value'] 
    life_star_basic['dist_st_value']=np.round(life_star_basic['dist_st_value'],2)
    #null value treatment: plx_value has masked entries therefore distance_values too
    #ref:
    life_star_basic['dist_st_ref']=ap.table.MaskedColumn(dtype=object,length=len(life_star_basic),
                                    mask=[True for j in range(len(life_star_basic))])
    life_star_basic['dist_st_ref'][np.where(life_star_basic['dist_st_value'].mask==False)]= \
            [life_provider['provider_name'][0] for j in range(len(
                life_star_basic['dist_st_ref'][np.where(life_star_basic['dist_st_value'].mask==False)]))]
    # can I do the same transformation with the errors? -> try on some examples and compare to simbad ones
    life_star_basic['coo_gal_err_angle']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_err_maj']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_err_min']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_qual']=['?'
                        for j in range(len(life_star_basic))]
    life_star_basic['main_id']=life_star_basic['main_id'].astype(str)
    # source
    # transformed from simbad ircs coordinates using astropy
    life_star_basic['coo_gal_ref']=ap.table.Column(dtype=object,length=len(life_star_basic))
    life_star_basic['coo_gal_ref']=life_provider['provider_name'][0] 
    #for all entries since coo_gal column not masked column
             
    life_star_basic['coo_gal_ref']=life_star_basic['coo_gal_ref'].astype(str)
    life_star_basic=life_star_basic['main_id','coo_gal_l','coo_gal_b','coo_gal_err_angle',
                                   'coo_gal_err_maj','coo_gal_err_min','coo_gal_qual',
                                   'coo_gal_ref','dist_st_value','dist_st_ref','sptype_string']
    

    life_star_basic=sptype_string_to_class(life_star_basic,life_provider['provider_name'][0])
    
    #-----------measurement tables -------------------------------------
    #applying model from E. E. Mamajek on SIMBAD spectral type

    [sim_objects]=hf.load(['sim_objects'],stringtoobjects=False)
    
    stars=sim_objects[np.where(sim_objects['type']=='st')]
    cat=ap.table.join(stars,life_star_basic)
    cat=spec(cat['main_id','sptype_string'])
    #if I take only st objects from sim_star_basic I don't loose objects during realspectype
    life_mes_teff_st=cat['main_id','mod_Teff']
    life_mes_teff_st.rename_column('mod_Teff','teff_st_value')
    life_mes_teff_st['teff_st_qual']=['D' for i in range(len(life_mes_teff_st))]
    life_mes_teff_st['teff_st_ref']=['2013ApJS..208....9P' for i in range(len(life_mes_teff_st))]
    
    life_mes_radius_st=cat['main_id','mod_R']
    life_mes_radius_st.rename_column('mod_R','radius_st_value')
    life_mes_radius_st['radius_st_qual']=['D' for i in range(len(life_mes_radius_st))]
    life_mes_radius_st['radius_st_ref']=['2013ApJS..208....9P' for i in range(len(life_mes_radius_st))]
    
    life_mes_mass_st=cat['main_id','mod_M']
    life_mes_mass_st.rename_column('mod_M','mass_st_value')
    life_mes_mass_st['mass_st_qual']=['D' for i in range(len(life_mes_mass_st))]
    life_mes_mass_st['mass_st_ref']=['2013ApJS..208....9P' for i in range(len(life_mes_mass_st))]
    
    #specifying stars cocerning multiplicity
    #main sequence simbad object type: MS*, MS? -> luminocity class
    #Interacting binaries and close CPM systems: **, **?
    
    #-----------------sources table-------------------------------------
    life_sources=ap.table.Table()
    tables=[life_provider,life_star_basic,life_mes_teff_st,life_mes_radius_st,life_mes_mass_st]
    ref_columns=[['provider_bibcode'],['coo_gal_ref'],['teff_st_ref'],['radius_st_ref'],['mass_st_ref']]
    for cat,ref in zip(tables,ref_columns):
        life_sources=sources_table(cat,ref,life_provider['provider_name'][0],life_sources)
    
    #removing this column because I had to adapt it where there was a leadin d entry but change not useful for db just for 
    #life parameter creation
    life_star_basic.remove_column('sptype_string')
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': life_list_of_tables[i]=life_sources
        if table_names[i]=='provider': life_list_of_tables[i]=life_provider
        if table_names[i]=='star_basic': life_list_of_tables[i]=life_star_basic
        if table_names[i]=='mes_teff_st': life_list_of_tables[i]=life_mes_teff_st
        if table_names[i]=='mes_radius_st': life_list_of_tables[i]=life_mes_radius_st
        if table_names[i]=='mes_mass_st': life_list_of_tables[i]=life_mes_mass_st
        hf.save([life_list_of_tables[i][:]],['life_'+table_names[i]])
    return life_list_of_tables



def provider_gaia(table_names,gaia_list_of_tables,distance_cut_in_pc,temp=True):
    """
    Obtains and arranges gaia data.
    
    Currently there is a provlem in obtaining the data through pyvo.
    A temporary method to ingest old gaia data was implemented and can be
    accessed by setting temp=True as argument.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param gaia_list_of_tables:
    :type gaia_list_of_tables: list(astropy.table.table.Table)
    :param distance_cut_in_pc:
    :param bool temp: 
    :returnss: List of astropy tables containing
        reference data, provider data, object data, identifier data,  
        stellar effective temperature, radius, mass and binarity data.
    :rtype: list(astropy.table.table.Table)
    """
    
    plx_in_mas_cut=1000./distance_cut_in_pc
    #making cut a bit bigger for correct treatment of objects on boundary
    plx_cut=plx_in_mas_cut-plx_in_mas_cut/10.
    #---------------define provider-------------------------------------
    gaia_provider=ap.table.Table()
    gaia_provider['provider_name']=['Gaia']
    gaia_provider['provider_url']=["https://gea.esac.esa.int/tap-server/tap"]
    gaia_provider['provider_bibcode']=['2016A&A...595A...1G']
    gaia_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
    print('Creating ',gaia_provider['provider_name'][0],' tables ...')
    
    #query
    adql_query="""
    SELECT s.source_id ,p.mass_flame, p.radius_flame,
        p.teff_gspphot, p.teff_gspspec, m.nss_solution_type 
    FROM gaiadr3.gaia_source as s
        JOIN gaiadr3.astrophysical_parameters as p ON s.source_id=p.source_id
            LEFT JOIN gaiadr3.nss_two_body_orbit as m ON s.source_id=m.source_id
    WHERE s.parallax >="""+str(plx_in_mas_cut)
    
    if temp: #because of bug in gaia server where async not working currently
        service = vo.dal.TAPService(gaia_provider['provider_url'][0])
        result=service.run_sync(adql_query.format(**locals()), maxrec=160000)
        gaia=result.to_table()

    else:
        gaia=query(gaia_provider['provider_url'][0],adql_query) 
        
    gaia.rename_columns(['mass_flame','radius_flame'],
                        ['mass_st_value','radius_st_value'])
    gaia['gaia_id']=['Gaia DR3 '+str(gaia['source_id'][j]) for j in range(len(gaia))]
    gaia['ref']=['2022arXiv220800211G' for j in range(len(gaia))]#dr3 paper
    
    #---------------gaia_ident-----------------------
    gaia_sim_idmatch=fetch_main_id(gaia['gaia_id','ref'],colname='gaia_id',
                           oid=False) 
    #should be gaia_id, main_id, ref minus 40 objects that have only gaia_id
    gaia_ident=gaia_sim_idmatch.copy()
    gaia_ident.rename_columns(['gaia_id','ref'],['id','id_ref'])
    gaia_ident['id_ref']=gaia_ident['id_ref'].astype(str)
    #creating simbad main_id ident rows
    sim_main_id_ident=ap.table.Table()
    sim_main_id_ident['main_id']=gaia_ident['main_id']
    sim_main_id_ident['id']=gaia_ident['main_id']
    sim_main_id_ident['id_ref']=['2000A&AS..143....9W' for j in range(len(gaia_ident))]
    gaia_ident=ap.table.vstack([gaia_ident,sim_main_id_ident])
    #now need to add the 40 objects that have only gaia_identifiers
    #for setdiff need both columns to be same type
    for col in gaia_sim_idmatch.colnames:
        gaia_sim_idmatch[col]=gaia_sim_idmatch[col].astype(str)
    gaia_only_id=ap.table.setdiff(gaia['gaia_id','ref'],gaia_sim_idmatch['gaia_id','ref'])
    gaia_only_id['main_id']=gaia_only_id['gaia_id']
    gaia_only_id.rename_columns(['gaia_id','ref'],['id','id_ref'])
    #for vstack need both columns to be same type
    for col in gaia_ident.colnames:
        gaia_ident[col]=gaia_ident[col].astype(str)
    gaia_ident=ap.table.vstack([gaia_ident,gaia_only_id])
    #add main_id to gaia table
    gaia=ap.table.join(gaia_ident['main_id','id'],gaia,
                       keys_left='id', keys_right='gaia_id')
    gaia.remove_column('id')

    #-----------------gaia_objects------------------
    gaia_objects=ap.table.Table(names=['main_id','ids'],dtype=[object,object])
    gaia_objects=ids_from_ident(gaia_ident['main_id','id'],gaia_objects)
    #grouped_gaia_ident=gaia_ident.group_by('main_id')
    #ind=grouped_gaia_ident.groups.indices
    #for i in range(len(ind)-1):
    # -1 is needed because else ind[i+1] is out of bonds
     #   ids=[]
      #  for j in range(ind[i],ind[i+1]):
       #     ids.append(grouped_gaia_ident['id'][j])
       # ids="|".join(ids)
        #gaia_objects.add_row([grouped_gaia_ident['main_id'][ind[i]],ids])
    gaia_objects['type']=['st' for j in range(len(gaia_objects))]
    gaia_objects['main_id']=gaia_objects['main_id'].astype(str)
    gaia_objects=ap.table.join(gaia_objects,gaia['main_id','nss_solution_type'],join_type='left')
    gaia_objects['type'][np.where(gaia_objects['nss_solution_type']!='')]=['sy' for j in range(len(
            gaia_objects['type'][np.where(gaia_objects['nss_solution_type']!='')]))]
    gaia_objects.remove_column('nss_solution_type')

    #gaia_mes_binary
    gaia_mes_binary=gaia_objects['main_id','type']
    print('tbd add binary flag True to children of system objects once I get h_link info from gaia')
    gaia_mes_binary.rename_column('type','binary_flag')
    gaia_mes_binary['binary_flag']=gaia_mes_binary['binary_flag'].astype(object)
    gaia_mes_binary=replace_value(gaia_mes_binary,'binary_flag','sy','True')
    gaia_mes_binary=replace_value(gaia_mes_binary,'binary_flag','st','False')
    gaia_mes_binary['binary_ref']=['2016A&A...595A...1G' for j in range(len(gaia_mes_binary))]
    gaia_mes_binary['binary_qual']=['B' if gaia_mes_binary['binary_flag'][j]=='True' \
                                    else 'E' for j in range(len(gaia_mes_binary))]
    #if necessary lower binary_qual for binary_flag = False to level of simbad.
    
    #gaia_mes_teff
    gaia_mes_teff_st=gaia['main_id','teff_gspphot']
    gaia_mes_teff_st['ref']=[gaia['ref'][j]+ ' GSP-Phot'
                         for j in range(len(gaia))]
    
    #remove masked rows
    gaia_mes_teff_st.remove_rows(gaia_mes_teff_st['teff_gspphot'].mask.nonzero()[0])
    gaia_mes_teff_st.rename_columns(['teff_gspphot','ref'],['teff_st_value','teff_st_ref'])
    
    temp=gaia['main_id','teff_gspspec']
    temp['ref']=[gaia['ref'][j]+ ' GSP-Spec'
                for j in range(len(gaia))]
    temp.remove_rows(temp['teff_gspspec'].mask.nonzero()[0])
    temp.rename_columns(['teff_gspspec','ref'],['teff_st_value','teff_st_ref'])
    
    gaia_mes_teff_st=ap.table.vstack([gaia_mes_teff_st,temp])
    gaia_mes_teff_st['teff_st_qual']=['B' for j in range(len(gaia_mes_teff_st))]
    gaia_mes_teff_st=gaia_mes_teff_st['main_id','teff_st_value',
                                      'teff_st_qual','teff_st_ref']
    
    #------------------gaia_mes_radius---------------
    gaia_mes_radius_st=gaia['main_id','radius_st_value','ref']
    gaia_mes_radius_st.remove_rows(gaia_mes_radius_st['radius_st_value'].mask.nonzero()[0])
    gaia_mes_radius_st['radius_st_qual']=['B' for j in range(len(gaia_mes_radius_st))]
    gaia_mes_radius_st['radius_st_ref']=[gaia_mes_radius_st['ref'][j] + ' FLAME'
                                   for j in range(len(gaia_mes_radius_st))]
    gaia_mes_radius_st.remove_column('ref')
    
    #------------------gaia_mes_mass---------------
    gaia_mes_mass_st=gaia['main_id','mass_st_value','ref']
    gaia_mes_mass_st.remove_rows(gaia_mes_mass_st['mass_st_value'].mask.nonzero()[0])
    gaia_mes_mass_st['mass_st_qual']=['B' for j in range(len(gaia_mes_mass_st))]
    gaia_mes_mass_st['mass_st_ref']=[gaia_mes_mass_st['ref'][j] + ' FLAME'
                                   for j in range(len(gaia_mes_mass_st))]
    gaia_mes_mass_st.remove_column('ref')
    
    #sources table
    gaia_sources=ap.table.Table()
    tables=[gaia_provider,gaia_ident,gaia_mes_teff_st,gaia_mes_radius_st,
            gaia_mes_mass_st,gaia_mes_binary]
    ref_columns=[['provider_bibcode'],['id_ref'],['teff_st_ref'],['radius_st_ref'],
                 ['mass_st_ref'],['binary_ref']]
    for cat,ref in zip(tables,ref_columns):
        gaia_sources=sources_table(cat,ref,gaia_provider['provider_name'][0],gaia_sources)
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': gaia_list_of_tables[i]=gaia_sources
        if table_names[i]=='provider': gaia_list_of_tables[i]=gaia_provider
        if table_names[i]=='objects': gaia_list_of_tables[i]=gaia_objects
        if table_names[i]=='ident': gaia_list_of_tables[i]=gaia_ident
        if table_names[i]=='mes_teff_st': gaia_list_of_tables[i]=gaia_mes_teff_st
        if table_names[i]=='mes_radius_st': gaia_list_of_tables[i]=gaia_mes_radius_st
        if table_names[i]=='mes_mass_st': gaia_list_of_tables[i]=gaia_mes_mass_st
        if table_names[i]=='mes_binary': gaia_list_of_tables[i]=gaia_mes_binary
        hf.save([gaia_list_of_tables[i][:]],['gaia_'+table_names[i]])
    return gaia_list_of_tables


def provider_wds(table_names,wds_list_of_tables,temp=False,test_objects=[]):
    """
    This function obtains and arranges wds data.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param wds_provider_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :type wds_provider_list_of_tables: list(astropy.table.table.Table)
    :param bool temp: Defaults to False. Used for debugging. saves querrying time.
    :returns: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data, basic stellar data and binarity data.
    :rtype:  list(astropy.table.table.Table)
    """
    
    #---------------define provider--------------------------------------------
    wds_provider=ap.table.Table()
    wds_provider['provider_name']=['WDS']
    wds_provider['provider_url']=["http://tapvizier.u-strasbg.fr/TAPVizieR/tap"]
    wds_provider['provider_bibcode']=['2001AJ....122.3466M']
    wds_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    #---------------define queries---------------------------------------------
    adql_query=["""SELECT
                    wds.WDS as wds_name, wds.Comp as wds_comp,
                    wds.sep1 as wds_sep1, wds.sep2 as wds_sep2, 
                    wds.Obs1 as wds_obs1, wds.Obs2 as wds_obs2
                    FROM "B/wds/wds" as wds """]
    
    #------------------querrying-----------------------------------------------
    print('Creating ',wds_provider['provider_name'][0],' tables ...')
    #perform query for objects with parallax >50mas
    test_objects=np.array(test_objects)
    if temp:
        print(' loading...')
        [wds]=hf.load(['wds'])
        #currently temp=True not giving same result because wds['system_main_id'][j] are '' and not masked
        for col in ['system_main_id','primary_main_id','secondary_main_id']:
            wds[col][np.where(wds[col]=='')]=np.ma.masked
        print('tbd: add provider_access of last query')
    else:
        print(' querying VizieR for WDS...')
        wds=query(wds_provider['provider_url'][0],adql_query[0])
        print('length query',len(wds))
        
        # I need to match the wds objects with the simbad ones to inforce the
        # distance cut since wds does not have distance information.
        
        # initializing and setting type for object comparison in later join 
        for col in ['sim_wds_id','system_name','primary','secondary']:
            wds[col]=wds['wds_name'].astype(object)
        
        # assigning correct name of system, primary and secondary for each wds object
        for j in range(len(wds)):
            if wds['wds_comp'][j]=='':#trivial binaries
                wds['system_name'][j]='WDS J'+wds['wds_name'][j]+'AB'
                #AB added since apparently simbad calls trivial binary system AB too
                wds['primary'][j]='WDS J'+wds['wds_name'][j]+'A'
                wds['secondary'][j]='WDS J'+wds['wds_name'][j]+'B'
            else:#higer order multiples
                wds['system_name'][j]='WDS J'+wds['wds_name'][j]+wds['wds_comp'][j]
                if len(wds['wds_comp'][j])==2:
                    wds['primary'][j]='WDS J'+wds['wds_name'][j]+wds['wds_comp'][j][0]
                    wds['secondary'][j]='WDS J'+wds['wds_name'][j]+wds['wds_comp'][j][1]
                else:
                    components=wds['wds_comp'][j].split(',')
                    wds['primary'][j]='WDS J'+wds['wds_name'][j]+components[0]
                    wds['secondary'][j]='WDS J'+wds['wds_name'][j]+components[1]
        print('number of trivial binary systems:',len(wds[np.where(wds['wds_comp']=='')]))
                
        if len(test_objects)>0:
            print('in wds as system_name', test_objects[np.where(np.in1d(test_objects,wds['system_name']))])
            print('in wds as primary',test_objects[np.where(np.in1d(test_objects,wds['primary']))])
            print('in wds as secondary', test_objects[np.where(np.in1d(test_objects,wds['secondary']))])


    # an alternative would be to query simbad for the main id and then cut by distance
    # this however takes way longer as it joins 150'000 elements
    #    wds=fetch_main_id(wds,colname='wds_full_name',name='main_id',oid=False)
    #    wds=distance_cut(wds,colname='wds_full_name',main_id=True)
        print(' performing distance cut...')
        
        #assigning main_id for system using sim_hlink and cutting on the system or the components
        wds_system_cut=distance_cut(wds,colname='system_name',main_id=False)
        wds_system_cut.rename_column('main_id','system_main_id')
        
        wds_primary_cut=distance_cut(wds,colname='primary',main_id=False)
        
        wds_secondary_cut=distance_cut(wds,colname='secondary',main_id=False)
        [sim_h_link]=hf.load(['sim_h_link'])
        #joining parent object
        wds_primary_cut=ap.table.join(wds_primary_cut,sim_h_link['main_id','parent_main_id'],
                                  keys='main_id',join_type='left')
        wds_primary_cut.rename_columns(['main_id','parent_main_id'],['primary_main_id','system_main_id'])
        
        wds_secondary_cut=ap.table.join(wds_secondary_cut,sim_h_link['main_id','parent_main_id'],
                                  keys='main_id',join_type='left')
        wds_secondary_cut.rename_columns(['main_id','parent_main_id'],['secondary_main_id','system_main_id'])
        #here some empty ones when child is known in simbad but no parent. in this case would I want to assign system_name in system main_id? do it later
        
        wds=ap.table.vstack([wds_system_cut,wds_primary_cut])
        wds=ap.table.vstack([wds,wds_secondary_cut])
        print('lenwds',len(wds))
                        
        if len(test_objects)>0:
            print(wds['system_main_id','primary_main_id','secondary_main_id'])
            print('in wds as system_main_id',test_objects[np.where(np.in1d(test_objects,wds['system_main_id']))])
            print('in wds as primary_main_id',test_objects[np.where(np.in1d(test_objects,wds['primary_main_id']))])
            print('in wds as secondary_main_id',test_objects[np.where(np.in1d(test_objects,wds['secondary_main_id']))])
        
        hf.save([wds],['wds'])
    
    wds['system_main_id']=wds['system_main_id'].astype(object)
    wds['system_name']=wds['system_name'].astype(object)

    
    #-----------------creating output table wds_ident and wds_h_link------------------------
    wds_ident=ap.table.Table(names=['main_id','id'],dtype=[object,object],masked=True)
    # create wds_h_link (for systems)
    wds_h_link=ap.table.Table(names=['main_id','parent_main_id'],dtype=[object,object])
    #add all relevant invormation
    #about identifiers
    table_main=['system_name','system_main_id','system_main_id',
               'primary','primary_main_id','primary_main_id',
               'secondary','secondary_main_id','secondary_main_id']
    table_id=['system_name','system_main_id','system_name',
             'primary','primary_main_id','primary',
             'secondary','secondary_main_id','secondary']
    empty=ap.table.Table(names=['main_id'],dtype=[object],masked=True)
    for id1,id2 in zip(table_main,table_id):
        temp=empty.copy()
        temp['main_id']=wds[id1].astype(object)
        temp['id']=wds[id2].astype(object)
        wds_ident=ap.table.vstack([wds_ident,temp])
    
    #about relations of objects
    table_main_id=['primary','primary','primary_main_id','primary_main_id',
                   'secondary','secondary','secondary_main_id','secondary_main_id']
    table_parent=['system_name','system_main_id','system_name','system_main_id',
                  'system_name','system_main_id','system_name','system_main_id']
    for id1,id2 in zip(table_main_id,table_parent):
        temp=empty.copy()
        temp['main_id']=wds[id1].astype(object)
        temp['parent_main_id']=wds[id2].astype(object)
        wds_h_link=ap.table.vstack([wds_h_link,temp])
    #delete all rows containing masked entries
    wds_ident.remove_rows(wds_ident['main_id'].mask.nonzero()[0])
    wds_ident.remove_rows(wds_ident['id'].mask.nonzero()[0])
    wds_h_link.remove_rows(wds_h_link['main_id'].mask.nonzero()[0])
    wds_h_link.remove_rows(wds_h_link['parent_main_id'].mask.nonzero()[0])

    #uniqueness
    wds_ident=ap.table.unique(wds_ident)
    wds_h_link=ap.table.unique(wds_h_link)

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

    wds_h_link=ap.table.unique(wds_h_link) 

    #refs
    wds_ident['id_ref']=[wds_provider['provider_bibcode'][0] for j in range(len(wds_ident))]
    wds_h_link['h_link_ref']=[wds_provider['provider_bibcode'][0] for j in range(len(wds_h_link))]

    #create objects table
    #create ids
    wds_objects=ap.table.Table(names=['main_id','ids'],dtype=[object,object])
    wds_objects=ids_from_ident(wds_ident['main_id','id'],wds_objects)
    #if it has children, it is type system
    #if it has no children it can either be star or close in system
    wds_objects['type']=['sy' for j in range(len(wds_objects))]
    #change to st for those that have no children
    wds_objects['type'][np.invert(np.in1d(wds_objects['main_id'],wds_h_link['parent_main_id']))]=['st' for j in range(len(
            [np.invert(np.in1d(wds_objects['main_id'],wds_h_link['parent_main_id']))]))] 
    
    if len(test_objects)>0:
        print('number of test objects that are in objects main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_objects['main_id']))])
    
    #-----------------creating output table wds_mes_binary------------------------
    wds_mes_binary=wds_objects['main_id','type']#[np.where(wds_objects['type']=='sy')]
    wds_mes_binary.rename_column('type','binary_flag')
    wds_mes_binary['binary_flag']=wds_mes_binary['binary_flag'].astype(object)
    wds_mes_binary['binary_flag']=['True' for j in range(len(wds_mes_binary))]
    wds_mes_binary['binary_ref']=[wds_provider['provider_bibcode'][0] for j in range(len(wds_mes_binary))]
    wds_mes_binary['binary_qual']=['C' for j in range(len(wds_mes_binary))]
    
    if len(test_objects)>0:
        print('number of test objects that are in mes_binary main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_mes_binary['main_id']))])
    
    #-----------------creating output table wds_mes_sep_ang------------------------
    #better join them
    wds_mes_sep_ang0=ap.table.join(wds['system_name','wds_sep1','wds_obs1','wds_sep2','wds_obs2'],
                                  wds_ident['main_id','id'],keys_left='system_name', keys_right='id')
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
    wds_mes_sep_ang=ap.table.vstack([wds_mes_sep_ang1,wds_mes_sep_ang2])
    #add a quality to sep1 which is better than sep2. because newer measurements should be better.
    wds_mes_sep_ang['sep_ang_ref']=[wds_provider['provider_bibcode'][0] for j in range(len(wds_mes_sep_ang))]
    #wds_mes_sep_ang.rename_column('system_main_id','main_id')
    #remove columns where sep_ang_value is masked
    wds_mes_sep_ang.remove_columns(wds_mes_sep_ang['sep_ang_value'].mask.nonzero()[0])
    #uniqueness where obs date not known 
    if len(wds_mes_sep_ang['sep_ang_obs_date'].mask.nonzero()[0])>0:
        unique_unknown_obs_date=ap.table.unique(wds_mes_sep_ang[np.where(
                wds_mes_sep_ang['sep_ang_obs_date'].mask.nonzero()[0])],keys=['main_id','sep_ang_value'])
        unique_known_obs_date=ap.table.unique(wds_mes_sep_ang[np.where(
                np.invert(wds_mes_sep_ang['sep_ang_obs_date'].mask.nonzero()[0]))],
                keys=['main_id','sep_ang_value','sep_ang_obs_date'])   
        wds_mes_sep_ang=ap.table.vstack([unique_unknown_obs_date,unique_known_obs_date])
    else:
        wds_mes_sep_ang=ap.table.unique(wds_mes_sep_ang)
        
    if len(test_objects)>0:
        print('number of test objects that are in mes_sep_ang main_id \n', \
              test_objects[np.where(np.in1d(test_objects,wds_mes_sep_ang['main_id']))])
    #--------------creating output table wds_sources --------------------------
    wds_sources=ap.table.Table()
    tables=[wds_provider,wds_ident]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['id_ref']]
    for cat,ref in zip(tables,ref_columns):
        wds_sources=sources_table(cat,ref,wds_provider['provider_name'][0],wds_sources)

    #saving tables
    for i in range(len(table_names)):
        if table_names[i]=='sources': wds_list_of_tables[i]=wds_sources
        if table_names[i]=='provider': wds_list_of_tables[i]=wds_provider  
        if table_names[i]=='objects': wds_list_of_tables[i]=wds_objects
        if table_names[i]=='ident': wds_list_of_tables[i]=wds_ident
        if table_names[i]=='h_link': wds_list_of_tables[i]=wds_h_link
        if table_names[i]=='mes_sep_ang': wds_list_of_tables[i]=wds_mes_sep_ang
        if table_names[i]=='mes_binary': wds_list_of_tables[i]=wds_mes_binary 
        hf.save([wds_list_of_tables[i][:]],['wds_'+table_names[i]])
    return wds_list_of_tables
