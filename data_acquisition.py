"""
Authors Note
This code was written by Franziska Menti in 2023.
It creates the tables for the LIFE target database.
For more information on LIFE please visit www.life-space-mission.com.
Beta warning: As a prototype this code will change. In particular
the addition of more providers, parameters as well as a possible
overhaul of the planet data is planned.
The document is structured as follows:
- Authors note
- Definition of functions
- Main code
"""
###############################################################################
#-------------------------definition of functions------------------------------
###############################################################################
#general
#importing libraries
import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
from datetime import datetime

#-------------------global helper functions------------------------------------
def save(cats,names):
    """
    This functions saves the tables given as python list in the cats parameter.
    The saving location is 'data/{name}.xml' where path is given in the paths
    parameter.
    :param cats: Python list of astropy table to be saved.
    :param names: Python list of strings containing names for saving location
        of tables in cats.
    """
    #go through all the elements in both lists
    for cat,path in zip(cats,names):
        #for each column header
        for i in cat.colnames:
            #if the type of the column is object (=adaptable length string)
            if cat[i].dtype == object:
                #transform the type into string
                cat[i] = cat[i].astype(str)
        #save the table
        ap.io.votable.writeto(
        	    ap.io.votable.from_table(cat), f'data/{path}.xml')
    return

def stringtoobject(cat,number=100):
    """
    This function changes from string to object format.
    The later has the advantace of allowing strings of varying length.
    Without it truncation of entries is a risk.
    :param cat: Astropy table.
    :param number: Length of longest string type element in the table.
        Default is 100.
    :return cat: Astropy table with all string columns transformed into
        object type ones.
    """
    #defining string types as calling them string does not work and instead
    #the type name <U3 is needed for a string of length 3
    stringtypes=[np.dtype(f'<U{j}') for j in range(1,number)]
    #for each column header
    for i in cat.colnames:
        #if the type of the column is string
        if cat[i].dtype in stringtypes:
            #transform the type into object
            cat[i] = cat[i].astype(object)
    return cat

def load(paths,stringtoobjects=True):
    """
    This function loads the tables saved in XML format at saving locations
    specified in paths. If stringtoobject is True the function 
    stringtoobjects is invoked.
    :param paths: Python list of saving locations.
    :return cats: Python list of loaded astropy tables.
    """
    #initialize return parameter as list
    cats=[]
    #go through all the elements in the paths list
    for path in paths:
        #read the saved data into the cats lists as astropy votable element
        to_append=ap.io.votable.parse_single_table(f'data/{path}.xml')
        cats.append(to_append.to_table())
    #go through all the tables in the cats list
    if stringtoobjects:
        for cat in cats:
            cat=stringtoobject(cat,3000)
    return cats

#-------------------initialization function------------------------------------
def initialize_database_tables(table_names,list_of_tables):
    """
    This function initializes the database tables with column name and data
    type specified but no actual data in them.
    :return list_of_tables: List of astropy tables in the order sources,
                objects, provider, ident (identifiers), h_link (relation between
                objects),star_basic,planet_basic, disk_basic, mes_mass_pl 
                (planetary mass measurements), mes_teff_st (stellar 
                effective temperature measurements), mes_radius_st (stellar
                radius measurements), mes_mass_st (stellar mass measurements),
                mes_binary (binarity information) and mes_sep_phys (stellar
                physical separation between binaries measurement).
    """
    #explanation of abbreviations: id stands for identifier, idref for
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
        #object idref, id, source idref for the id parameter, id reference
        names=['object_idref','id','id_source_idref','id_ref'],
        dtype=[int,object,int,object])

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
               'sep_phys_value','sep_phys_err','sep_phys_qual',
               'sep_phys_source_idref','sep_phys_ref'],
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
               float,float,object,#sep_phys
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
               'rad_source_iderf','rad_ref'],
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
    
    mes_sep_phys=ap.table.Table(
        names=['object_idref',
               'sep_phys_value','sep_phys_err','sep_phys_qual',
               'sep_phys_source_idref','sep_phys_ref'],
        dtype=[int,
               float,float,object,#sep_phys
               int,object])

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
        if table_names[i]=='mes_sep_phys': list_of_tables[i]=mes_sep_phys
        save([list_of_tables[i]],['empty_'+table_names[i]])
    return list_of_tables

#------------------------------provider helper functions-----------------------
def query(link,query,catalogs=[]):
    """
    Performs a query via TAP on the service given in the link parameter.
    If a list of tables is given in the catalogs parameter,
    those are uploaded to the service beforehand.
    :param link: Service access URL.
    :param query: Query to be asked of the external database service in ADQL.
    :param catalogs: List of astropy tables to be uploaded to the service.
    :return cat: Astropy table containing the result of the query.
    """
    #defining the vo service using the given link
    service = vo.dal.TAPService(link)
    #without upload tables
    if catalogs==[]:
        result=service.run_async(query.format(**locals()), maxrec=160000)
    #with upload tables
    else:
        tables={}
        for i in range(len(catalogs)):
            tables.update({f"t{i+1}":catalogs[i]})
        result = service.run_async(query,uploads=tables,timeout=None,
                                   maxrec=160000)
    cat=result.to_table()
    return cat

def sources_table(cat,ref_columns,provider,old_sources=ap.table.Table()):
    """
    This function creates or updates the source table out of the given
    references. The entries are unique and the columns consist out of the
    reference and provider_name.
    :param cat: Astropy table on which the references should be gathered.
    :param ref_columns: Header of the columns containing reference information.
    :param provider: String for provider name.
    :param old_sources: Previously created reference table.
    :return sources: Astropy table containing references and provider
        information.
    """
    cat_sources=ap.table.Table() #table initialization
    cat_reflist=[] #initialization of list to store reference information
    #for all the columns given add reference information to the cat_reflist
    for k in range(len(ref_columns)):
        #In case the column has elements that are masked skip those
        if type(cat[ref_columns[k]])==ap.table.column.MaskedColumn:
            cat_reflist.extend(
            	cat[ref_columns[k]][np.where(cat[ref_columns[k]].mask==False)])
        else:
            cat_reflist.extend(cat[ref_columns[k]])
    #add list of collected references to the table and call the column ref
    cat_sources['ref']=cat_reflist
    cat_sources=ap.table.unique(cat_sources)#keeps only unique values
    #attaches service information
    cat_sources['provider_name']=[provider for j in range(len(cat_sources))]
    #combine old and new sources into one table
    sources=ap.table.vstack([old_sources,cat_sources])
    sources=ap.table.unique(sources) #remove double entries
    return sources

def fetch_main_id(cat,colname='oid',name='main_id',oid=True):
    """
    Joins main_id from simbad to the column colname. Returns the whole
    table cat but without any rows where no simbad main_id was found.
    :param cat: Astropy table containing column colname.
    :param colname: Column header of the identifiers that should be searched
        in SIMBAD.
    :param name: Column header for the SIMBAD main identifiers, default is
        main_id.
    :param oid: Specifies wether colname is a SIMBAD oid or normal identifier.
    :return cat: Astropy table with all main SIMBAD identifiers that could 
        be found in column "name".
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

def nullvalues(cat,colname,nullvalue,verbose=False):
    """
    This function transforms all masked entries of the column colname of
    the table cat into the values given in nullvalue.
    :param cat: Astropy table containing the column colname.
    :param colname: String of the name of an astropy table column.
    :param nullvalue: Value to be placed instead of masked elements.
    :param verboxe: Bool default False, in case of True prints message
        if the column is not an astropy masked column.
    :return cat: Astropy table with masked elements of colname replaced
        by nullvalue.
    """
    if type(cat[colname])==ap.table.column.MaskedColumn:
                cat[colname].fill_value=nullvalue
                cat[colname]=cat[colname].filled()
    elif verbose:
        print(colname,'is no masked column')
    return cat

def lowerquality(cat,colname):
    """
    This function lowers the quality (A highest, E lowest) entry of the 
    in colname specified column of the astropy table cat.
    :param cat: Astropy table containing the column colname and only 
        rows that should be altered.
    :param colname: String of the name of an astropy table column.
    :return cat: Astropy table with lowered quality entries in the column
        colname.
    """
    for i in range(len(cat)):
        if cat[colname][i]=='A':
            cat[colname][i]='B'
        elif cat[colname][i]=='B':
            cat[colname][i]='C'
        elif cat[colname][i]=='C':
            cat[colname][i]='D'
        elif cat[colname][i]=='D':
            cat[colname][i]='E'
    return cat

def replace_value(cat,column,value,replace_by):
    """
    This function replaces the in the parameter value
    specified entries of the column colname in the table
    cat with the in replaced_by specified entry.
    :param cat: Astropy table with column colname.
    :param column: String designating column in which to replace
        the entries.
    :param value: Entry to be replaced.
    :param replace_by: Entry to be put in place of param value.
    :return cat: Astropy table with replaced entries.
    """
    cat[column][np.where(cat[column]==value)]= \
            [replace_by for i in range(
        len(cat[column][np.where(cat[column]==value)]))]
    return cat

#-----------------------------provider data ingestion--------------------------
def provider_simbad(table_names,sim_list_of_tables):
    """
    This function obtains the SIMBAD data and arranges it in a way
    easy to ingest into the database.
    :param table_names: List of strings containing the names for the 
        output tables.
    :param sim_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :return simbad_list_of_tables: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data, basic stellar data, distance measurement data
        and binarity data.
    """
    #---------------define provider--------------------------------------------
    sim_provider=ap.table.Table()
    sim_provider['provider_name']=['SIMBAD']
    sim_provider['provider_url']=["http://simbad.u-strasbg.fr:80/simbad/sim-tap"]
    sim_provider['provider_bibcode']=['2000A&AS..143....9W']
    sim_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    #---------------define queries---------------------------------------------
    select="""SELECT b.main_id,b.ra AS coo_ra,b.dec AS coo_dec,
        b.coo_err_angle, b.coo_err_maj, b.coo_err_min,b.oid,
        b.coo_bibcode AS coo_ref, b.coo_qual,b.sp_type AS sptype_string,
        b.sp_qual AS sptype_qual, b.sp_bibcode AS sptype_ref,
        b.plx_err, b.plx_value, b.plx_bibcode AS plx_ref,b.plx_qual,
        h_link.membership, h_link.parent AS parent_oid,
        h_link.link_bibcode AS h_link_ref, a.otypes,ids.ids,
        f.I as mag_i_value, f.J as mag_j_value
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
    #------------------querrying-----------------------------------------------
    print('Creating ',sim_provider['provider_name'][0],' tables ...')
    #perform query for objects with parallax >50mas
    simbad=query(sim_provider['provider_url'][0],adql_query[0])
    #querries parent and children objects with no parallax value
    parents_without_plx=query(sim_provider['provider_url'][0],upload_query[0],[simbad])
    children_without_plx=query(sim_provider['provider_url'][0],upload_query[1],[simbad])
    #adding of no_parallax objects to rest of simbad query objects
    simbad=ap.table.vstack([simbad,parents_without_plx])
    simbad=ap.table.vstack([simbad,children_without_plx])
    #----------------------sorting object types--------------------------------
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
    #removing any objects that are neither planet, star or system in type
    if to_remove_list!=[]:
        simbad.remove_rows(to_remove_list)
        print('removed',len(removed_otypes),' objects that had object types:',
              list(set(removed_otypes)))

    #creating helpter table stars
    temp_stars=simbad[np.where(simbad['type']!='pl')]
    #removing double objects (in there due to multiple parents)
    stars=ap.table.Table(ap.table.unique(temp_stars,keys='main_id'),copy=True)

    #-----------------creating output table sim_ident--------------------------
    sim_ident=query(sim_provider['provider_url'][0],upload_query[3],
                    [simbad['oid','main_id'][:].copy()])
    sim_ident['id_ref']=[sim_provider['provider_bibcode'][0] for j in range(len(sim_ident))]
    sim_ident.remove_column('oid')
    
    #--------------creating output table sim_h_link ---------------------------
    sim_h_link=simbad['main_id','parent_oid','h_link_ref','membership']
    #sim_h_link=nullvalues(sim_h_link,'parent_oid',0,verbose=False)
    ###sim_h_link=nullvalues(sim_h_link,'membership',-1,verbose=False)
    sim_h_link=fetch_main_id(sim_h_link,'parent_oid','parent_main_id')
    sim_h_link.remove_column('parent_oid')
    #typeconversion needed as smallint fill value != int null value
    sim_h_link['membership']=sim_h_link['membership'].astype(int)
    sim_h_link=nullvalues(sim_h_link,'membership',999999)

    # binary_flag 'True' for all stars with parents
    # meaning stars[main_id] in sim_h_link[child_main_id] -> stars[binary_flag]=='True'
    #could do this via two for loops but maybe easier way? maybe join. 
    for i_star in range(len(stars['main_id'])):
        for i_child in range(len(sim_h_link['main_id'])):
            if stars['main_id'][i_star]==sim_h_link['main_id'][i_child]:
                stars['binary_flag'][i_star]=='True'
                
    #--------------------creating helper table sim_stars-----------------------
    #change null value of plx_qual
    stars['plx_qual']=stars['plx_qual'].astype(object)
    stars=replace_value(stars,'plx_qual','',stars['plx_qual'].fill_value)
    #initiate some of the ref columns
    stars['mag_i_ref']=ap.table.MaskedColumn(dtype=object,length=len(stars),
                                    mask=[True for j in range(len(stars))])
    stars['mag_j_ref']=ap.table.MaskedColumn(dtype=object,length=len(stars),
                                    mask=[True for j in range(len(stars))])
    #add simbad reference where no other is given
    stars['mag_i_ref'][np.where(stars['mag_i_value'].mask==False)]=[
            sim_provider['provider_bibcode'][0] for j in range(len(
            stars['mag_i_ref'][np.where(stars['mag_i_value'].mask==False)]))]
    stars['mag_j_ref'][np.where(stars['mag_j_value'].mask==False)]=[
            sim_provider['provider_bibcode'][0] for j in range(len(
            stars['mag_j_ref'][np.where(stars['mag_j_value'].mask==False)]))]
    stars=replace_value(stars,'plx_ref','',sim_provider['provider_bibcode'][0])
    sim_h_link=replace_value(sim_h_link,'h_link_ref','',
                             sim_provider['provider_bibcode'][0])
        
    stars['binary_ref']=[sim_provider['provider_bibcode'][0] for j in range(len(stars))]
    stars['binary_qual']=['C' for j in range(len(stars))]
    #-----------------creating output table sim_planets------------------------
    temp_sim_planets=simbad['main_id','ids',
                            'type'][np.where(simbad['type']=='pl')]
    sim_planets=ap.table.Table(ap.table.unique(
                    temp_sim_planets,keys='main_id'),copy=True)
    #-----------------creating output table sim_objects------------------------
    sim_objects=ap.table.vstack([sim_planets['main_id','ids','type'],
                             stars['main_id','ids','type']])
    sim_objects['ids']=sim_objects['ids'].astype(object)
    #--------------creating output table sim_sources --------------------------
    sim_sources=ap.table.Table()
    tables=[sim_provider,stars, sim_h_link,sim_ident]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['coo_ref','plx_ref','mag_i_ref',
                    'mag_j_ref','binary_ref'],['h_link_ref'],
                    ['id_ref']]
    for cat,ref in zip(tables,ref_columns):
        sim_sources=sources_table(cat,ref,sim_provider['provider_name'][0],sim_sources)
    #------------------------creating output table sim_star_basic--------------
    sim_star_basic=stars['main_id','coo_ra','coo_dec','coo_err_angle',
                         'coo_err_maj','coo_err_min','coo_qual','coo_ref',
                         'mag_i_value','mag_i_ref','mag_j_value','mag_j_ref',
                         'sptype_string','sptype_qual','sptype_ref',
                         'plx_value','plx_err','plx_qual','plx_ref']
    #-----------creating mes_binary table------------
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
        save([sim_list_of_tables[i]],['sim_'+table_names[i]])
        
    return sim_list_of_tables

def provider_gk(table_names,gk_list_of_tables):
    """
    This function obtains the disk data and arranges it in a way
    easy to ingest into the database.
    :param table_names: List of strings containing the names for the 
        output tables.
    :param gk_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :return gk_list_of_tables: List of astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data and basic disk data.
    """
    #---------------define provider--------------------------------------------
    gk_provider=ap.table.Table()
    gk_provider['provider_name']=['Grant Kennedy Disks']
    gk_provider['provider_url']=['http://drgmk.com/sdb/']
    gk_provider['provider_bibcode']=['priv. comm.']
    gk_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
    print('Creating ',gk_provider['provider_name'][0],' tables ...')
    #loading table obtained via direct communication from Grant Kennedy
    gk_disks=ap.io.votable.parse_single_table(
        "data/additional_data/Grant_absil_2013_.xml").to_table()
    #transforming from string type into object to have variable length
    gk_disks=stringtoobject(gk_disks,212)
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
    #--------------creating output table gk_h_link ----------------------------
    gk_h_link=gk_disks['id','main_id','disks_ref']
    gk_h_link.rename_columns(['main_id','disks_ref'],
                             ['parent_main_id','h_link_ref'])
    gk_h_link.rename_column('id','main_id')
    #--------------creating output table gk_objects ---------------------------
    gk_disks['ids']=gk_disks['id']#because only single id per source given
    gk_objects=gk_disks['id','ids','type']
    gk_objects.rename_column('id','main_id')
    #--------------creating output table gk_ident -----------------------------
    gk_ident=gk_disks['ids','id','disks_ref']
    # would prefer to use id instad of ids paremeter but this raises an error
    # so I use ids which has the same content as id
    gk_ident.rename_columns(['ids','disks_ref'],['main_id','id_ref'])
    #--------------creating output table gk_sources ---------------------------
    gk_sources=ap.table.Table()
    tables=[gk_provider,gk_disks]
    #define header name of columns containing references data
    ref_columns=[['provider_bibcode'],['disks_ref']]
    for cat,ref in zip(tables,ref_columns):
        gk_sources=sources_table(cat,ref,gk_provider['provider_name'][0],gk_sources)
    #--------------creating output table gk_disk_basic-------------------------
    gk_disk_basic=gk_disks['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    for column in ['rdisk_bb','e_rdisk_bb']:
        #replacing 'None' with 'nan' as the first one is not float convertible
        temp_length=len(gk_disk_basic[column][np.where(
                        gk_disk_basic[column]=='None')])
        gk_disk_basic=replace_value(gk_disk_basic,column,'None','nan')
        gk_disk_basic[column].fill_value='nan' #because defeault is None and not float convertible
        #though this poses the issue that the float default float fill_value is 1e20
        gk_disk_basic[column].filled()
        gk_disk_basic[column]=gk_disk_basic[column].astype(float)
    gk_disk_basic.rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    gk_disk_basic=gk_disk_basic[np.where(gk_disk_basic['rad_value']!='nan')]
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': gk_list_of_tables[i]=gk_sources
        if table_names[i]=='provider': gk_list_of_tables[i]=gk_provider
        if table_names[i]=='objects': gk_list_of_tables[i]=gk_objects
        if table_names[i]=='ident': gk_list_of_tables[i]=gk_ident
        if table_names[i]=='h_link': gk_list_of_tables[i]=gk_h_link
        if table_names[i]=='disk_basic': gk_list_of_tables[i]=gk_disk_basic
        save([gk_list_of_tables[i]],['gk_'+table_names[i]])
    return gk_list_of_tables

def provider_exo(table_names,exo_list_of_tables,temp=True):
    """
    This function obtains the exomercat data and arranges it in a way easy to
    ingest into the database. Currently the exomercat server is not online.
    A temporary method to ingest old exomercat data was implemented and can be
    accessed by setting temp=True as argument.
    :param table_names: List of strings containing the names for the 
        output tables.
    :param exo_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :param temp: Bool with default value True determining if the exomercat
        data gets queried (False) or loaded from an old version (True).
    :return exo_list_of_tables: List of astropy table containing
        reference data, object data, identifier data, object to object
        relation data, basic planetary data and planetary mass measurement 
        data.
    """
    #---------------define provider--------------------------------------------
    exo_provider=ap.table.Table()
    exo_provider['provider_name']=['Exo-MerCat']
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    exo_provider['provider_bibcode']=['2020A&C....3100370A']
    
    
    print('Creating ',exo_provider['provider_name'][0],' tables ...')
    #---------------define query-----------------------------------------------
    adql_query="""SELECT *
                  FROM exomercat.exomercat"""
    #---------------obtain data------------------------------------------------
    if temp:
        exomercat=ap.io.ascii.read("data/additional_data/exo-mercat05-02-2023_v2.0.csv")
        exomercat=stringtoobject(exomercat,3000)
        exo_provider['provider_access']=['2023-02-05']

    else:
        exomercat=query(exo_provider['provider_url'][0],adql_query)
        exo_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    #----------------putting object main identifiers together-------------------

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

    def sort_out_20pc(cat,colname):
        """
        This Function sorts out objects not within the provider_simbad
        distance cut (previously 20pc hence the name). 
        :param cat: Astropy table to be matched against sim_objects table.
        :param colname: Name of the column to use for the match.
        :return cat: Table like cat without any objects not found in sim_objects.
        """
        [sim_objects]=load(['sim_objects'])
        sim_objects.rename_column('main_id','temp')
        cat=ap.table.join(cat,sim_objects['temp','ids'],
                          keys_left=colname,keys_right='temp')
        cat.remove_columns(['temp','ids'])
        return cat

    exo=exomercat
    exomercat=sort_out_20pc(exomercat,'host_main_id')

    # removing whitespace in front of main_id and name.
    # done after sort_out_20pc function to prevent missing values error
    for i in range(len(exomercat)):
        exomercat['planet_main_id'][i]=exomercat['planet_main_id'][i].strip()
        exomercat['main_id'][i]=exomercat['main_id'][i].strip()
        exomercat['name'][i]=exomercat['name'][i].strip()

    #show which elements from exomercat were not found in sim_objects
    exo['name']=exo['name'].astype(object)
    removed_objects=ap.table.setdiff(exo,exomercat,keys=['name'])
    save([removed_objects],['exomercat_removed_objects'])

    #-------------exo_ident---------------
    exo_ident=exomercat['planet_main_id','name']
    exo_ident.rename_columns(['planet_main_id','name'],['main_id','id'])
    for i in range(len(exomercat)):
        if exomercat['planet_main_id'][i]!=exomercat['name'][i]:
            exo_ident.add_row([exomercat['planet_main_id'][i],
                               exomercat['planet_main_id'][i]])
    exo_ident['id_ref']=[exo_provider['provider_bibcode'][0] for j in range(len(exo_ident))]
    # TBD: I have a wrong double object
    print("""TBD: I have a wrong double object because of different amount of white
          spaces between catalog and number""")
    print(exo_ident[np.where(exo_ident['main_id']=='Wolf  940 b')])

    #-------------exo_objects---------------
    # tbd at one point: I think I want to add hosts to object
    exo_objects=ap.table.Table(names=['main_id','ids'],dtype=[object,object])
    grouped_exo_ident=exo_ident.group_by('main_id')
    ind=grouped_exo_ident.groups.indices
    for i in range(len(ind)-1):
    # -1 is needed because else ind[i+1] is out of bonds
        ids=[]
        for j in range(ind[i],ind[i+1]):
            ids.append(grouped_exo_ident['id'][j])
        ids="|".join(ids)
        exo_objects.add_row([grouped_exo_ident['main_id'][ind[i]],ids])
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
        else:
            if type(exomercat['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_min'][i]==np.inf:
                exomercat['mass_pl_rel'][i]='>'
                exomercat['mass_pl_err'][i]=exomercat['mass_max'][i]
            else:
                exomercat['mass_pl_rel'][i]='='
                exomercat['mass_pl_err'][i]=max(exomercat['mass_max'][i],
                                        exomercat['mass_min'][i])
    exo_mes_mass_pl=exomercat['planet_main_id','mass','mass_pl_err','mass_url',
                            'mass_pl_rel','mass_pl_qual']
    exo_mes_mass_pl.rename_columns(['planet_main_id','mass','mass_url'],
                                    ['main_id','mass_pl_value','mass_pl_ref'])
    #remove masked rows
    exo_mes_mass_pl.remove_rows(exo_mes_mass_pl['mass_pl_value'].mask.nonzero()[0])
    #remove null values
    exo_mes_mass_pl=exo_mes_mass_pl[np.where(exo_mes_mass_pl['mass_pl_value']!=1e+20)]


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
        save([exo_list_of_tables[i]],['exo_'+table_names[i]])
    return exo_list_of_tables

def provider_life(table_names,life_list_of_tables):
    """
    This function loads the SIMBAD data obtained by the function provider_simbad
    and postprocesses it to provide more useful information. It uses a model
    from Eric E. Mamajek to predict temperature, mass and radius from the simbad 
    spectral type data.
    :param table_names: List of strings containing the names for the 
        output tables.
    :param life_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :return life_list_of_tables: List of astropy table containing
        reference data, provider data, basic stellar data, stellar effective
        temperature data, stellar radius data and stellar mass data.
    """
    #---------------define provider--------------------------------------------
    life_provider=ap.table.Table()
    life_provider['provider_name']=['LIFE']
    life_provider['provider_url']=['www.life-space-mission.com']
    life_provider['provider_bibcode']=['2022A&A...664A..21Q']
    life_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
    print('Creating ',life_provider['provider_name'][0],' tables ...')
    #---------------------star_basic----------------
    #galactic coordinates:  transformed from simbad ircs coordinates using astropy
    [life_star_basic]=load(['sim_star_basic'])
    ircs_coord=ap.coordinates.SkyCoord(\
            ra=life_star_basic['coo_ra'],dec=life_star_basic['coo_dec'],frame='icrs')
    gal_coord=ircs_coord.galactic
    life_star_basic['coo_gal_l']=gal_coord.l.deg*ap.units.degree
    life_star_basic['coo_gal_b']=gal_coord.b.deg*ap.units.degree
    life_star_basic['dist_st_value']=1000./life_star_basic['plx_value'] 
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
    
    def sptype_string_to_class(temp,ref):
        """
        This function extracts the temperature class, temperature class number
        and luminocity class information from the spectral type string (e.g. 
        M5V to M, 5 and V). It stores that information in the for this purpose
        generated new columns. Only objects of temperature class O, B, A, F,
        G, K, and M are processed. Only objects of luminocity class IV, V and VI
        are processed.
        :param temp: Astropy table containing spectral type information in
            the column sptype_string.
        :param ref: String designating origin of data.
        :return temp: Astropy table like temp with additional columns class_temp,
            class_temp_nr, class_lum and class_ref.
        """
        temp['class_temp']=ap.table.MaskedColumn(dtype=object,length=len(temp))
        temp['class_temp_nr']=ap.table.MaskedColumn(dtype=object,length=len(temp))
        temp['class_lum']=ap.table.MaskedColumn(dtype=object,length=len(temp))
        temp['class_ref']=ap.table.MaskedColumn(dtype=object,length=len(temp))
        for i in range(len(temp)):
            #sorting out objects like M5V+K7V
            if (len(temp['sptype_string'][i].split('+'))==1 and
            #sorting out entries like '', DA2.9, T1V
                    len(temp['sptype_string'][i])>0 and 
                    temp['sptype_string'][i][0] in ['O','B','A','F','G','K','M']):
                temp['class_temp'][i]=temp['sptype_string'][i][0]
                temp['class_ref'][i]=ref
                #sorting out objects like DA2.9
                if len(temp['sptype_string'][i])>1 and temp['sptype_string'][i][1] in ['0','1','2','3','4','5','6','7','8','9']:
                    temp['class_temp_nr'][i]=temp['sptype_string'][i][1]
                    #distinguishing between objects like K5V and K5.5V
                    if len(temp['sptype_string'][i])>2 and temp['sptype_string'][i][2]=='.':
                        temp['class_temp_nr'][i]=temp['sptype_string'][i][1:4]
                        if len(temp['sptype_string'][i])>4 and temp['sptype_string'][i][4] in ['I','V']:
                            temp['class_lum'][i]=temp['sptype_string'][i][4]
                            if len(temp['sptype_string'][i])>5 and temp['sptype_string'][i][5] in ['I','V']:
                                temp['class_lum'][i]=temp['sptype_string'][i][4:6]
                                if len(temp['sptype_string'][i])>6 and temp['sptype_string'][i][6] in ['I','V']:
                                    temp['class_lum'][i]=temp['sptype_string'][i][4:7]
                        else:
                            temp['class_lum'][i]='?'
                    elif len(temp['sptype_string'][i])>2 and temp['sptype_string'][i][2] in ['I','V']:
                        temp['class_lum'][i]=temp['sptype_string'][i][2]
                        if len(temp['sptype_string'][i])>3 and temp['sptype_string'][i][3] in ['I','V']:
                            temp['class_lum'][i]=temp['sptype_string'][i][2:4]
                            if len(temp['sptype_string'][i])>4 and temp['sptype_string'][i][4] in ['I','V']:
                                temp['class_lum'][i]=temp['sptype_string'][i][2:5]
                    else:
                        temp['class_lum'][i]='?'
            else:
                temp['class_temp'][i]='?'
                temp['class_temp_nr'][i]='?'
                temp['class_lum'][i]='?'
                temp['class_ref'][i]='?'
        return temp
    life_star_basic=sptype_string_to_class(life_star_basic,life_provider['provider_name'][0])
    
    #-----------measurement tables -----------------
    #applying model from E. E. Mamajek on SIMBAD spectral type
    
             
    def realspectype(cat):
        """
        Removes rows of cat where elements in column named 'sim_sptype' are
        either '', 'nan' or start with an other letter than the main sequence
        spectral type classification.
        :param cat: astropy table containing 'sim_sptype' column
        :return cat: astropy table, param cat with undesired rows removed
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
        Loads the table of Eric E. Mamajek containing stellar parameters 
        modeled from spectral types. Cleans up the columns for spectral 
        type, effective temperature radius and mass.
        :return votable: astropy table of the 4 parameters as columns
        """
        EEM_table=ap.io.ascii.read("data/Mamajek2022-04-16.csv")['#SpT','Teff','R_Rsun','Msun']
        EEM_table.rename_columns(['R_Rsun','Msun'],['Radius','Mass'])
        EEM_table=replace_value(EEM_table,'Radius',' ...','nan')
        EEM_table=replace_value(EEM_table,'Mass',' ...','nan')
        EEM_table=replace_value(EEM_table,'Mass',' ....','nan')
        EEM_table['Teff'].unit=ap.units.K
        EEM_table['Radius'].unit=ap.units.R_sun
        EEM_table['Mass'].unit=ap.units.M_sun       
        ap.io.votable.writeto(ap.io.votable.from_table(EEM_table), \
                              f'data/model_param.xml')#saving votable
        return EEM_table

    def match_sptype(cat,model_param,sptypestring='sim_sptype',teffstring='mod_Teff',\
                     rstring='mod_R',mstring='mod_M'):
        """
        Matches the spectral types with the ones in Mamajek's table and 
        includes the modeled effective Temperature,
        stellar radius and stellar mass into the catalog.
        :param cat: astropy table containing spectral type information
        :param sptypestring: string of column name where the spectral 
            type information is located
        :param teffstring: string of column name for the stellar effective 
            temperature column
        :param rstring: string of column name for the stellar radius column
        :param mstring: string of column name for the stellar mass column
        :return cat: the astropy table cat with added new columns for 
            effective temperature, radius and mass filled with model values
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
                for i in range(len(model_param['#SpT'])): 
                    #match first two letters
                    if model_param['#SpT'][i][:2]==cat[sptypestring][j][:2]: 
                            cat[teffstring][j]=model_param['Teff'][i]
                            cat[rstring][j]=model_param['Radius'][i]
                            cat[mstring][j]=model_param['Mass'][i]
                #as the model does not cover all spectral types on .5 accuracy, check those separately
                if cat[sptypestring][j][2:4]=='.5':
                    for i in range(len(model_param['#SpT'])):
                        # match first four letters
                        if model_param['#SpT'][i][:4]==cat[sptypestring][j][:4]:
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
        :return cat: Catalog of mainsequence stars with unique 
            simbad names, no binary subtypes and modeled parameters.
        """    
        cat=realspectype(cat)
        #model_param=ap.io.votable.parse_single_table(\
            #f"catalogs/model_param.xml").to_table()
        mp=model_param()#create model table as votable
        cat=match_sptype(cat,mp,sptypestring='sptype_string')
        cat.remove_rows([np.where(cat['mod_Teff'].mask==True)])
        cat.remove_rows([np.where(np.isnan(cat['mod_Teff']))])
        cat=ap.table.unique(cat, keys='main_id')
        return cat

    [sim_objects]=load(['sim_objects'],stringtoobjects=False)
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
    
    #-----------------sources table----------------------
    life_sources=ap.table.Table()
    tables=[life_provider,life_star_basic,life_mes_teff_st,life_mes_radius_st,life_mes_mass_st]
    ref_columns=[['provider_bibcode'],['coo_gal_ref'],['teff_st_ref'],['radius_st_ref'],['mass_st_ref']]
    for cat,ref in zip(tables,ref_columns):
        life_sources=sources_table(cat,ref,life_provider['provider_name'][0],life_sources)
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': life_list_of_tables[i]=life_sources
        if table_names[i]=='provider': life_list_of_tables[i]=life_provider
        if table_names[i]=='star_basic': life_list_of_tables[i]=life_star_basic
        if table_names[i]=='mes_teff_st': life_list_of_tables[i]=life_mes_teff_st
        if table_names[i]=='mes_radius_st': life_list_of_tables[i]=life_mes_radius_st
        if table_names[i]=='mes_mass_st': life_list_of_tables[i]=life_mes_mass_st
        save([life_list_of_tables[i]],['life_'+table_names[i]])
    return life_list_of_tables

def provider_gaia(table_names,gaia_list_of_tables,temp=True):
    """
    This function obtains the gaia data and arranges it in a way
    easy to ingest into the database. Currently there is a provlem 
    in obtaining the data through pyvo.
    A temporary method to ingest old gaia data was implemented and can be
    accessed by setting temp=True as argument.
    :return gaia_list_of_tables: List of astropy tables containing
        reference data, object data and
        basic stellar data.
    """
    #---------------define provider--------------------------------------------
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
    
    if temp:
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
    grouped_gaia_ident=gaia_ident.group_by('main_id')
    ind=grouped_gaia_ident.groups.indices
    for i in range(len(ind)-1):
    # -1 is needed because else ind[i+1] is out of bonds
        ids=[]
        for j in range(ind[i],ind[i+1]):
            ids.append(grouped_gaia_ident['id'][j])
        ids="|".join(ids)
        gaia_objects.add_row([grouped_gaia_ident['main_id'][ind[i]],ids])
    gaia_objects['type']=['None' for j in range(len(gaia_objects))]
    gaia_objects['main_id']=gaia_objects['main_id'].astype(str)
    gaia_objects=ap.table.join(gaia_objects,gaia['main_id','nss_solution_type'],join_type='left')
    gaia_objects['type'][np.where(gaia_objects['nss_solution_type']!='')]='sy'
    gaia_objects.remove_column('nss_solution_type')
    #there might be issue in building merging now

    #gaia_mes_binary
    gaia_mes_binary=gaia_objects['main_id','type'][np.where(gaia_objects['type']=='sy')]
    gaia_mes_binary.rename_column('type','binary_flag')
    gaia_mes_binary['binary_flag']=['True' for j in range(len(gaia_mes_binary))]
    gaia_mes_binary['binary_ref']=['2016A&A...595A...1G' for j in range(len(gaia_mes_binary))]
    gaia_mes_binary['binary_qual']=['B' for j in range(len(gaia_mes_binary))]
    #there might be issue in building merging now
    
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
        save([gaia_list_of_tables[i]],['gaia_'+table_names[i]])
    return gaia_list_of_tables

def provider_orb6(table_names,orb6_list_of_tables):
    
    orb6_provider=ap.table.Table()
    orb6_provider['provider_name']=['ORB6']
    orb6_provider['provider_url']=["http://tapvizier.cds.unistra.fr/TAPVizieR/tap"]
    orb6_provider['provider_bibcode']=['https://crf.usno.navy.mil/wds-orb6/']
    orb6_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
    print('Creating ',orb6_provider['provider_name'][0],' tables ...')
    
    #query
    adql_query="""
    SELECT "id1-DR3" as gaia_id, Axis, e_Axis, Name, WDS
    FROM "J/MNRAS/517/2925/tablea3" 
    WHERE "plx1-DR3" >="""+str(plx_in_mas_cut)
    
    orb6=query(orb6_provider['provider_url'][0],adql_query)
    orb6['gaia_id']=['Gaia DR3 '+str(orb6['gaia_id'][j]) for j in range(len(orb6))]
    orb6=fetch_main_id(orb6,colname='gaia_id',name='main_id',oid=False)
    orb6['ref']=['ORB6' for j in range(len(orb6))]
    orb6['ref']=orb6['ref'].astype(object)
    
    #------------ident-------------------
    orb6_ident=orb6['main_id','Name','ref']
    orb6_ident.rename_columns(['Name','ref'],['id','id_ref'])
    #add main_id main_id sim_ref for completeness
    sim_main_id=orb6_ident.copy()
    sim_main_id['id_ref']='2000A&AS..143....9W'
    #this is wrong, provider needs to be simbad or life not orb6 for simbad main_id
    #though I do have same issue with gaia stuff
    
    sim_main_id['id']=sim_main_id['main_id']
    orb6_ident=ap.table.vstack([orb6_ident,sim_main_id])
    orb6_ident=ap.table.unique(orb6_ident)
    #------------objects---------------------
    orb6_objects=ap.table.Table(names=['main_id','ids'],dtype=[object,object])
    grouped_orb6_ident=orb6_ident.group_by('main_id')
    ind=grouped_orb6_ident.groups.indices
    for i in range(len(ind)-1):
    # -1 is needed because else ind[i+1] is out of bonds
        ids=[]
        for j in range(ind[i],ind[i+1]):
            ids.append(grouped_orb6_ident['id'][j])
        ids="|".join(ids)
        orb6_objects.add_row([grouped_orb6_ident['main_id'][ind[i]],ids])
    orb6_objects['type']=['sy' for j in range(len(orb6_objects))]
        
    #---------------mes_binary-------------------------------------
    orb6_mes_binary=orb6['main_id','Axis','e_Axis','ref']
    orb6_mes_binary.rename_columns(['Axis','e_Axis','ref'],
                                   ['sep_phys_value','sep_phys_err','sep_phys_ref'])
    orb6_mes_binary['binary_flag']=['True' for j in range(len(orb6_mes_binary))]
    orb6_mes_binary['binary_ref']=orb6['ref']
    orb6_mes_binary['binary_qual']=['B' for j in range(len(orb6_mes_binary))]
    orb6_mes_binary['sep_phys_qual']=['B' for j in range(len(orb6_mes_binary))]
    orb6_mes_sep_phys=orb6_mes_binary['main_id','sep_phys_value','sep_phys_err',
                                      'sep_phys_qual','sep_phys_ref']
    orb6_mes_sep_phys=ap.table.unique(orb6_mes_sep_phys,silent=True)
    orb6_mes_sep_phys[orb6_mes_sep_phys['sep_phys_err'].mask.nonzero()[0]]=lowerquality(
            orb6_mes_sep_phys[orb6_mes_sep_phys['sep_phys_err'].mask.nonzero()[0]],'sep_phys_qual')
    
    orb6_mes_binary.remove_columns(['sep_phys_value','sep_phys_err',
                                      'sep_phys_qual','sep_phys_ref'])
    orb6_mes_binary=ap.table.unique(orb6_mes_binary,silent=True)
    #---------------sources---------------------------------------
    orb6_sources=ap.table.Table()
    tables=[orb6_provider,orb6_ident,orb6_mes_binary,orb6_mes_sep_phys]
    ref_columns=[['provider_bibcode'],['id_ref'],['binary_ref'],['sep_phys_ref']]
    for cat,ref in zip(tables,ref_columns):
        orb6_sources=sources_table(cat,ref,orb6_provider['provider_name'][0],orb6_sources)
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': orb6_list_of_tables[i]=orb6_sources
        if table_names[i]=='provider': orb6_list_of_tables[i]=orb6_provider
        if table_names[i]=='objects': orb6_list_of_tables[i]=orb6_objects
        if table_names[i]=='ident': orb6_list_of_tables[i]=orb6_ident
        if table_names[i]=='mes_binary': orb6_list_of_tables[i]=orb6_mes_binary
        if table_names[i]=='mes_sep_phys': orb6_list_of_tables[i]=orb6_mes_sep_phys
        save([orb6_list_of_tables[i]],['orb6_'+table_names[i]])
    return orb6_list_of_tables


#------------------------provider combining-----------------
def building(providers,table_names,list_of_tables):
    """
    This function builds from the input parameters the tables
    for the LIFE database.
    :param providers: List of astropy tables containing simbad, 
        grant kennedy, exomercat, gaia and orb6 data.
    :param table_names: List of strings corresponding to the names of the 
        astropy tables contained in providers and the return list.
    :param list_of_tables: List of empty astropy tables to be filled
        in and returned.
    :return list_of_tables: List of astropy table containing data 
        combined from the different providers.
    """
    #creates empty tables as needed for final database ingestion
    init=initialize_database_tables(table_names,list_of_tables)
    n_tables=len(init)

    cat=[ap.table.Table() for i in range(n_tables)]
    #for the sources and objects joins tables from different providers
    
    print(f'Building {table_names[0]} table ...')
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

    def idsjoin(cat,column_ids1,column_ids2):
        """
        This function merges the identifiers from two different columns
        into one.
        :param cat: Astropy table containing two identifer columns.
        :param column_ids1: column name for the first identifier column
        :param column_ids2: column name for the second identifier column
        :return cat: Astropy table like input cat but with only one identifier
            column containing the unique identifiers of both of the
            previous identifier columns.
        """
        #initializing column
        cat['ids']=ap.table.Column(dtype=object, length=len(cat))
        for column in [column_ids1,column_ids2]:
            cat=nullvalues(cat,column,'')
        for i in range(len(cat)):
            #splitting object into list of elements
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
        This function merges the data from the same physical object obtained
        from different providers into one entry.
        :param cat: Astropy table containing multiple entries for the same
            physical objects due to data from different providers.
        :return cat: Astropy table with unique object entries.
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
    
    print(f'Building {table_names[1]} table ...')
    
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
    def match(cat,sources,paras,provider):
        """
        This function joins the source identifiers to the in paras specified
        parameters of cat.
        :param cat: Astropy table with empty para_source_id columns.
        :param sources: Astropy table containing reference data.
        :param paras: String describing a parameter in cat.
        :param provider: Name of the data provider as string.
        :return cat: Astropy table containing para_source_id data.
        """
        #for all parameters specified
        for para in paras:
            #if they have reference columns
            if para+'_ref' in cat.colnames:
                #if those reference columns are masked
                cat=nullvalues(cat,para+'_ref','None')
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
    
    print(f'Building {table_names[2]} table ...')
    
    paras=[['id'],['h_link'],['coo','plx','dist_st','coo_gal','mag_i','mag_j','class'],
           ['mass_pl'],['rad'],['dist_st'],['mass_pl'],['teff_st'],
          ['radius_st'],['mass_st'],['binary'],['sep_phys']]
    
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
        
        for j in range(len(providers)):#for the different providers
            if len(providers[j][i])>0:
                #replacing ref with source_idref columns
                #getting source_idref to each ref
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
        
        #if resulting catalog not empty
        if len(cat[i])>0:
            #only keeping unique entries
            cat[i]=ap.table.unique(cat[i],silent=True)
            
            def best_para(para,mes_table):
                """
                This function creates a table containing only the highest
                quality measurement for each object.
                :param para: string describing parameter e.g. mass
                :param mes_table: Astropy table containing only columns
                    'main_id', para+'_value',para+'_err',para+'_qual' and para+'_ref'
                """
                if para=='binary':
                    columns=['main_id',para+'_flag',para+'_qual',para+'_source_idref']
                elif para=='mass_pl':
                    columns=['main_id',para+'_value',para+'_rel',para+'_err',para+'_qual',para+'_source_idref']
                else:
                    columns=['main_id',para+'_value',para+'_err',para+'_qual',para+'_source_idref']
                mes_table=mes_table[columns[0:]]
                best_para=mes_table[columns[0:]][:0].copy()
                #group mes_table by object (=main_id)
                grouped_mes_table=mes_table.group_by('main_id')
                #take highest quality
                for j in range(len(grouped_mes_table.groups.keys)):#go through all objects
                    for qual in ['A','B','C','D','E','?']:
                        for i in range(len(grouped_mes_table.groups[j])):
                            if grouped_mes_table.groups[j][i][para+'_qual']==qual:
                                best_para.add_row(grouped_mes_table.groups[j][i])
                                break
                        else:
                            continue  # only executed if the inner loop did NOT break
                        break  # only executed if the inner loop DID break
                    
                return best_para
                #I have quite some for loops here, will be slow
                #this function is ready for testing

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

            #for all the other tables add object_idref
            else:
                #first remove the object_idref we got from empty initialization
                #though I would prefer a more elegant way to do this
                cat[i].remove_column('object_idref') 
                cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],
                                     join_type='left')
                cat[i].rename_column('object_id','object_idref')
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
                temp.rename_columns(['object_id','main_id'],
                                    ['object_idref','temp'])
                cat[i]=ap.table.join(cat[i],temp,join_type='outer',
                                     keys='object_idref')
                cat[i].remove_column('temp')
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
            if table_names[i]=='mes_sep_phys':
                sep_phys_best_para=best_para('sep_phys',cat[i])
                cat[5].remove_columns(['sep_phys_value','sep_phys_err',
                                      'sep_phys_qual','sep_phys_source_idref',
                                      'sep_phys_ref'])
                cat[5]=ap.table.join(cat[5],sep_phys_best_para,join_type='left')
            cat[i]=cat[i].filled()
            cat[i]=ap.table.unique(cat[i])
        else:
            print('error: empty table',i,table_names[i])
    cat[5]=cat[5].filled()
    
    #unify null values (had 'N' and '?' because of ap default fill_value and type conversion string vs object)
    tables=[cat[table_names.index('star_basic')],cat[table_names.index('planet_basic')],
            cat[table_names.index('disk_basic')],cat[table_names.index('mes_mass_pl')],
            cat[table_names.index('mes_teff_st')],cat[table_names.index('mes_radius_st')],
            cat[table_names.index('mes_mass_st')],cat[table_names.index('mes_binary')]]
    columns=[['coo_qual','coo_gal_qual','plx_qual','dist_st_qual',
              'sep_phys_qual','teff_st_qual','radius_st_qual',
              'mass_st_qual','binary_qual'],
             ['mass_pl_qual','mass_pl_rel'],
             ['rad_qual','rad_rel'],['mass_pl_qual','mass_pl_rel'],
             ['teff_st_qual'],['radius_st_qual'],['mass_st_qual'],['binary_qual']]
    for i in range(len(tables)):
        for col in columns[i]:
            tables[i]=replace_value(tables[i],col,'N','?')
            tables[i]=replace_value(tables[i],col,'N/A','?')
    
    save(cat,table_names)
    return cat

###############################################################################
#-------------------------Main code--------------------------------------------
###############################################################################

#------------------------initialize empty database tables----------------------
table_names=['sources','objects','provider','ident','h_link','star_basic',
              'planet_basic','disk_basic','mes_mass_pl',
              'mes_teff_st','mes_radius_st','mes_mass_st','mes_binary','mes_sep_phys']

#distance cut
distance_cut_in_pc=5#25.

#transforming from pc distance cut into parallax in mas cut
plx_in_mas_cut=1000./distance_cut_in_pc
#making cut a bit bigger for correct treatment of objects on boundary
plx_cut=plx_in_mas_cut-plx_in_mas_cut/10.

#------------------------obtain data from external sources---------------------
empty_provider=[ap.table.Table() for i in range(len(table_names))]

sim=provider_simbad(table_names,empty_provider[:])

gk=provider_gk(table_names,empty_provider[:])

exo=provider_exo(table_names,empty_provider[:])

life=provider_life(table_names,empty_provider[:])

gaia=provider_gaia(table_names,empty_provider[:])

orb6=provider_orb6(table_names,empty_provider[:])

#------------------------combine data from external sources-------------------
database_tables=building([sim[:],gk[:],exo[:],life[:],gaia[:],orb6[:]],table_names,empty_provider[:])