"""
Authors Note
This code was written by Franziska Menti in 2023.
It creates the tables for the prototype of the LIFE target database.
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
import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables

#distance cut
distance_cut_in_pc=25.
#-------------------global helper functions------------------------------------
def save(cats,paths):
    """
    This functions saves the tables given as python list in the cats parameter.
    The saving location is 'data/{path}.xml' where path is given in the paths
    parameter.
    :param cats: Python list of astropy table to be saved.
    :param paths: Python list of paths to where to save the tables given in
    	cats.
    """
    #go through all the elements in both lists
    for cat,path in zip(cats,paths):
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
    Without it strings can get truncated.
    :param cat: Astropy table.
    :param number: Length of longest string type element in the table.
        Default is 100.
    :return cat: Astropy table.
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
    specified in paths.
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

def nullvalues(cat,colname,nullvalue,verbose=False):
    if type(cat[colname])==ap.table.column.MaskedColumn:
                cat[colname].fill_value=nullvalue
                cat[colname]=cat[colname].filled()
    elif verbose:
        print(colname,'is no masked column')
    return cat

#-------------------initialization function------------------------------------
def initialize_database_tables():
    """
    This function initializes the database tables with column name and data
    type specified but no actual data in them.
    :return list_of_tables: List of astropy tables in the order sources,
                objects, ident (identifiers), h_link (relation between
                objects),star_basic,planet_basic, disk_basic,
                mes_dist (distance measurements) and
                mes_mass_st (mass measurements).
    """
    #explanation of abbreviations: id stands for identifier, idref for
    # reference identifier and parameter_source_idref for the identifier in the
    # source table corresponding to the mentioned parameter

    list_of_tables=[]

    sources=ap.table.Table(
        #reference,...
        names=['ref','provider_name','provider_url','provider_bibcode',
               'source_id'],
        dtype=[object,object,object,object,int])
    list_of_tables.append(sources)

    objects=ap.table.Table(
        #object id, type of object, all identifiers, main id
        names=['object_id','type','ids','main_id'],
        dtype=[int,object,object,object])
    list_of_tables.append(objects)

    #identifier table
    ident=ap.table.Table(
        #object idref, id, source idref for the id parameter, id reference
        names=['object_idref','id','id_source_idref','id_ref'],
        dtype=[int,object,int,object])
    list_of_tables.append(ident)

    #hierarchical link table (which means relation between objects)
    h_link=ap.table.Table(
        #child object idref (e.g. planet X),
        #parent object idref (e.g. host star of planet X)
        #source idref of h_link parameter, h_link reference,
        #membership probability
        names=['child_object_idref','parent_object_idref',
               'h_link_source_idref','h_link_ref','membership'],
        dtype=[int,int,int,object,int])
    list_of_tables.append(h_link)

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
               'dist_value','dist_err','dist_qual','dist_source_idref',
               'dist_ref',
               'sptype_value','sptype_err','sptype_qual','sptype_source_idref',
               'sptype_ref',
               'teff_st_value','teff_st_err','teff_st_qual','teff_st_source_idref',
               'teff_st_ref',
               'radius_st_value','radius_st_err','radius_st_qual','radius_st_source_idref',
               'radius_st_ref',
               'mass_st_value','mass_st_err','mass_st_qual','mass_st_source_idref',
               'mass_st_ref',
               'binary_flag','binary_source_idref','binary_ref',
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
               float,float,object,int,#teff
               object,
               float,float,object,int,#rad
               object,
               float,float,object,int,#mass
               object,
               bool,int,object,#binary
               float,float,object,#sep_phys
               int,object])
    list_of_tables.append(star_basic)

    planet_basic=ap.table.Table(
        #object idref, mass value, mass error, mass realtion (min, max, equal),
        #mass quality, source idref of mass parameter, mass reference
        names=['object_idref','mass_pl_val','mass_pl_err','mass_pl_rel','mass_pl_qual',
               'mass_pl_source_idref','mass_pl_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(planet_basic)

    disk_basic=ap.table.Table(
        #object idref, black body radius value, bbr error,
        #bbr relation (min, max, equal), bbr quality,...
        names=['object_idref','rad_value','rad_err','rad_rel','rad_qual',
               'rad_source_iderf','rad_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(disk_basic)

    mes_dist=ap.table.Table(
        names=['object_idref','dist_value','dist_err','dist_qual',
               'dist_source_idref','dist_ref'],
        dtype=[int,float,float,object,
               int,object])
    list_of_tables.append(mes_dist)

    mes_mass_pl=ap.table.Table(
        names=['object_idref','mass_pl_value','mass_pl_err','mass_pl_rel','mass_pl_qual',
               'mass_pl_source_idref','mass_pl_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(mes_mass_pl)
    
    mes_teff_st=ap.table.Table(
        names=['object_idref','teff_st_value','teff_st_err','teff_st_rel','teff_st_qual',
               'teff_st_source_idref','teff_st_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(mes_teff_st)
    
    mes_radius_st=ap.table.Table(
        names=['object_idref','radius_st_value','radius_st_err','radius_st_rel','radius_st_qual',
               'radius_st_source_idref','radius_st_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(mes_radius_st)
    
    mes_mass_st=ap.table.Table(
        names=['object_idref','mass_st_value','mass_st_err','mass_st_rel','mass_st_qual',
               'mass_st_source_idref','mass_st_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(mes_mass_st)

    #save all tables
    save(list_of_tables,
         ['empty_sources','empty_objects','empty_ident','empty_h_link',
         'empty_star_basic','empty_planet_basic','empty_disk_basic',
         'empty_mes_dist','empty_mes_mass_pl','empty_mes_teff_st',
         'empty_mes_radius_st','empty_mes_mass_st'])
    return list_of_tables

#transforming from pc distance cut into parallax in mas cut
plx_in_mas_cut=1000./distance_cut_in_pc
#making cut a bit bigger for correct treatment of objects on boundary
plx_cut=plx_in_mas_cut-plx_in_mas_cut/10.

#------------------------------provider helper functions-----------------------
def query(link,query,catalogs=[]):
    """
    Performs a query via TAP on the service given in the link parameter.
    If a list of tables is given in the catalogs parameter,
    those are uploaded to the service.
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
    reference, provider_name, provider_url and provider_bibcode.
    :param cat: Astropy table on which the references should be gathered.
    :param ref_columns: Header of the columns containing reference information.
    :param provider: List containing name, url and bibcode of provider.
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
    cat_sources['provider_name']=[provider[0]]*len(cat_sources)
    cat_sources['provider_url']=[provider[1]]*len(cat_sources)
    cat_sources['provider_bibcode']=[provider[2]]*len(cat_sources)
    #combine old and new sources into one table
    sources=ap.table.vstack([old_sources,cat_sources])
    sources=ap.table.unique(sources) #remove double entries
    return sources

def fetch_main_id(cat,colname='oid',name='main_id',oid=True):
    """
    Joins main_id from simbad to the column colname. Returns the whole
    table cat but without any rows where no simbad main_id was found.
    :param cat: Astropy table.
    :param colname: Column header of the identifiers that should be searched
        in SIMBAD.
    :param name: Column header for the SIMBAD main identifiers, default is
        main_id.
    :param oid: Specifies wether colname is a SIMBAD oid or normal identifier.
    :return cat: Astropy table with all identifiers that could be found
        in SIMBAD and column contining ther main identifier.
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

#-----------------------------provider data ingestion--------------------------
def provider_simbad():
    """
    This function obtains the SIMBAD data and arranges it in a way
    easy to ingest into the database.
    :return simbad_list_of_tables: List of astropy tables containing
        reference data, object data, identifier data, object to object
        relation data, basic stellar data and distance measurement data.
    """
    #---------------define provider--------------------------------------------
    TAP_service="http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    provider_name='SIMBAD'
    provider_bibcode='2000A&AS..143....9W'
    #---------------define queries---------------------------------------------
    select="""SELECT b.main_id,b.ra AS coo_ra,b.dec AS coo_dec,
        b.coo_err_angle, b.coo_err_maj, b.coo_err_min,b.oid,
        b.coo_bibcode AS coo_ref, b.coo_qual,b.sp_type AS sptype_value,
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
        """SELECT oid, dist AS dist_value, plus_err, qual AS dist_qual,
        bibcode AS dist_ref,minus_err,dist_prec
        FROM mesDistance
        JOIN TAP_UPLOAD.t1 ON oidref=t1.oid""",
        #query all identifiers for objects in TAP_UPLOAD.t1 table
        """SELECT id, t1.*
        FROM ident
        JOIN TAP_UPLOAD.t1 ON oidref=t1.oid"""]
    #define header name of columns containing references data
    ref_columns=[['coo_ref','plx_ref','mag_i_ref','mag_j_ref'],['h_link_ref'],
                 ['dist_ref'],['id_ref']]
    #------------------querrying-----------------------------------------------
    print(f'Creating {provider_name} tables ...')
    #perform query for objects with parallax >50mas
    simbad=query(TAP_service,adql_query[0])
    #querries parent and children objects with no parallax value
    parents_without_plx=query(TAP_service,upload_query[0],[simbad])
    children_without_plx=query(TAP_service,upload_query[1],[simbad])
    #adding of no_parallax objects to rest of simbad query objects
    simbad=ap.table.vstack([simbad,parents_without_plx])
    simbad=ap.table.vstack([simbad,children_without_plx])
    #----------------------sorting object types--------------------------------
    #sorting from object type into star, system and planet type
    simbad['type']=['None' for i in range(len(simbad))]
    simbad['multiple']=[False for i in range(len(simbad))]
    to_remove_list=[]
    for i in range(len(simbad)):
        #planets
        if "Pl" in simbad['otypes'][i]:
            simbad['type'][i]='pl'
        #stars
        elif "*" in simbad['otypes'][i]:
            #system containing multiple stars
            if "**" in simbad['otypes'][i]:
                simbad['type'][i]='sy'
                simbad['multiple'][i]=True
            #individual stars
            else:
                simbad['type'][i]='st'
        else:
            print('Removed one object because its type was',
                  simbad['otypes'][i],
                  'which is neither planet, star nor system.')
            #most likely single brown dwarfs
            #storing information for later removal from table called simbad
            to_remove_list.append(i)
    #removing any objects that are neither planet, star or system in type
    if to_remove_list!=[]:
        simbad.remove_rows(to_remove_list)

    #creating helpter table stars
    temp_stars=simbad[np.where(simbad['type']!='pl')]
    #removing double objects (in there due to multiple parents)
    stars=ap.table.Table(ap.table.unique(temp_stars,keys='main_id'),copy=True)

    #-----------------creating output table sim_ident--------------------------
    sim_ident=query(TAP_service,upload_query[3],
                    [simbad['oid','main_id'][:].copy()])
    sim_ident['id_ref']=[provider_bibcode for j in range(len(sim_ident))]
    sim_ident.remove_column('oid')

    #-------------------creating output table sim_mes_dist---------------------
    sim_mes_dist=query(TAP_service,upload_query[2],[stars[:].copy()])
    sim_mes_dist=fetch_main_id(sim_mes_dist)
    sim_mes_dist['dist_err']=np.maximum(sim_mes_dist['plus_err'],
                                       -sim_mes_dist['minus_err'])
    sim_mes_dist.remove_rows(sim_mes_dist['dist_err'].mask.nonzero()[0])
    #change provider given quality null values all to same: 'N'
    sim_mes_dist['dist_qual'][np.where(sim_mes_dist['dist_qual']=='')]='N'
    sim_mes_dist['dist_qual'][np.where(sim_mes_dist['dist_qual']==':')]='N'
    #group by oid
    grouped_mes_dist=sim_mes_dist.group_by('main_id')
    best_mes_dist=sim_mes_dist['main_id','dist_value','plus_err',
                             'dist_qual','dist_ref'][:0]
    best_mes_dist.rename_column('plus_err','dist_err')
    for i in range(len(grouped_mes_dist.groups.keys)):
        #sort by quality
        row=grouped_mes_dist.groups[i][np.where(
                    grouped_mes_dist['dist_prec'].groups[i]
                    ==np.max(grouped_mes_dist['dist_prec'].groups[i]))][0]
        #take first and add to best_paras
        #which error to take when there are multiples...
        best_mes_dist.add_row([row['main_id'],row['dist_value'],
                              row['dist_err'],row['dist_qual'],
                              row['dist_ref']])
    #join with other multimes thingis
    best_paras=best_mes_dist
    # TBD when more multi measurement tables are implemented: vstack them here
    sim_mes_dist=sim_mes_dist['main_id','dist_value','dist_err',
                            'dist_qual','dist_ref']

    #--------------------creating helper table sim_stars-----------------------
    #add best para from multiple measurements tables
    stars=ap.table.join(stars,best_paras,keys='main_id',join_type='left')
    #change null value of plx_qual from '' to 'N'
    stars['plx_qual'][np.where(stars['plx_qual']=='')]='N'
    stars['dist_qual'][np.where(stars['dist_qual']=='')]='N'
    stars['dist_qual'][np.where(stars['dist_qual']==':')]='N'
    stars=nullvalues(stars,'mag_i_value',0)
    stars=nullvalues(stars,'mag_j_value',0)
    #first need to initiate those columns
    stars['mag_i_ref']=ap.table.MaskedColumn(dtype=object,length=len(stars),
                                    mask=[True for j in range(len(stars))])
    stars['mag_j_ref']=ap.table.MaskedColumn(dtype=object,length=len(stars),
                                    mask=[True for j in range(len(stars))])
    stars['mag_i_ref'][np.where(stars['mag_i_value']!=0)]=provider_bibcode
    stars['mag_j_ref'][np.where(stars['mag_j_value']!=0)]=provider_bibcode
    
    #--------------creating output table sim_h_link ---------------------------
    sim_h_link=simbad['main_id','parent_oid','h_link_ref','membership']
    # if you want to exclude objects with lower membership probability
    # use this line instead:
    # sim_h_link=simbad['main_id','parent_oid','h_link_ref',
    #                   'membership'][np.where(simbad['membership']>50)]
    # consequence is that you loose objects with no membership value given (~)
    # e.g. alf cen system
    sim_h_link=nullvalues(sim_h_link,'parent_oid',0,verbose=False)
    sim_h_link=nullvalues(sim_h_link,'membership',-1,verbose=False)
    sim_h_link=fetch_main_id(sim_h_link,'parent_oid','parent_main_id')
    sim_h_link.remove_column('parent_oid')
    #null values
    sim_h_link['membership'].fill_value=-1
    sim_h_link['membership']=sim_h_link['membership'].filled()
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
    tables=[stars, sim_h_link, sim_mes_dist,sim_ident]
    for cat,ref in zip(tables,ref_columns):
        sim_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],sim_sources)
    simbad_list_of_tables=[sim_sources,sim_objects, sim_ident, sim_h_link]
    #------------------------creating output table sim_star_basic--------------
    sim_star_basic=stars['main_id','coo_ra','coo_dec','coo_err_angle',
                         'coo_err_maj','coo_err_min','coo_qual','coo_ref',
                         'mag_i_value','mag_i_ref','mag_j_value','mag_j_ref',
                         'sptype_value','sptype_qual','sptype_ref',
                         'plx_value','plx_err','plx_qual','plx_ref',
                         'dist_value','dist_err','dist_qual','dist_ref']
    sim_star_basic['sptype_value']=sim_star_basic['sptype_value'].astype(str)
    sim_star_basic['sptype_qual']=sim_star_basic['sptype_qual'].astype(str)
    sim_star_basic['sptype_ref']=sim_star_basic['sptype_ref'].astype(str)
    
    simbad_list_of_tables.extend([sim_star_basic,sim_mes_dist])
    save(simbad_list_of_tables,['sim_sources','sim_objects','sim_ident',
                                'sim_h_link','sim_star_basic','sim_mes_dist'])
    return simbad_list_of_tables

def provider_gk():
    """
    This function obtains the disk data and arranges it in a way
    easy to ingest into the database.
    :return gk_list_of_tables: List of astropy tables containing
        reference data, object data, identifier data, object to object
        relation data and basic disk data.
    """
    #---------------define provider--------------------------------------------
    provider_name='priv. comm.'
    TAP_service='None'
    provider_bibcode='None'
    print(f'Creating {provider_name} tables ...')
    #loading table obtained via direct communication from Grant Kennedy
    gk_disks=ap.io.votable.parse_single_table(
        "data/Grant_absil_2013_.xml").to_table()
    #transforming from string type into object to have variable length
    gk_disks=stringtoobject(gk_disks,212)
    #removing objects with plx_value=='None'
    gk_disks=gk_disks[np.where(gk_disks['plx_value']!='None')]
    #converting masked plx_value into -99
    gk_disks['plx_value'].fill_value=-99
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
    gk_sources=sources_table(gk_disks,['disks_ref'],[provider_name,TAP_service,
                                           provider_bibcode])
    gk_list_of_tables=[gk_sources,gk_objects,gk_ident,gk_h_link]
    #--------------creating output table gk_disk_basic-------------------------
    gk_disk_basic=gk_disks['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    for column in ['rdisk_bb','e_rdisk_bb']:
        #replacing 'None' with 'nan' as the first one is not float convertible
        temp_length=len(gk_disk_basic[column][np.where(
                        gk_disk_basic[column]=='None')])
        gk_disk_basic[column][np.where(gk_disk_basic[column]=='None')]=[
                                        'nan' for i in range(temp_length)]
        gk_disk_basic[column].fill_value='nan'
        gk_disk_basic[column].filled()
        gk_disk_basic[column]=gk_disk_basic[column].astype(float)
    gk_disk_basic.rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    gk_list_of_tables.append(gk_disk_basic)
    save(gk_list_of_tables,
         ['gk_sources','gk_objects', 'gk_ident', 'gk_h_link','gk_disk_basic'])
    return gk_list_of_tables

def provider_exo(temp=True):
    """
    This function obtains the exomercat data and arranges it in a way easy to
    ingest into the database. Currently the exomercat server is not online.
    A temporary method to ingest old exomercat data was implemented and can be
    accessed by setting temp=True as argument.
    :return exo_list_of_tables: List of astropy table containing
        reference data, object data, identifier data, object to object
        relation data, basic planetary data and mass measurement data.
    """
    #---------------define provider--------------------------------------------
    TAP_service="http://archives.ia2.inaf.it/vo/tap/projects"
    provider_name='Exo-MerCat'
    provider_bibcode='2020A&C....3100370A'
    #---------------define query-----------------------------------------------
    adql_query="""SELECT *
                  FROM exomercat.exomercat"""
    #---------------obtain data------------------------------------------------
    print(f'Creating {provider_name} tables ...')
    if temp:
        exomercat=ap.io.ascii.read("data/exomercat_Sep2.csv")
        exomercat=stringtoobject(exomercat,3000)

    else:
        exomercat=query(TAP_service,adql_query)
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
        This Function sorts out objects not within 20pc. The value comes from
        the LIFE database distance cut.
        :param cat: Input table to be matched against sim_objects table.
        :param colname: Name of the column to use for the match.
        :return cat: Table like cat without any objects not found in sim_objects
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
    exo_ident['id_ref']=[provider_bibcode for j in range(len(exo_ident))]
    # TBD: I have a wrong double object
    # print(exo_ident[np.where(exo_ident['main_id']=='Wolf  940 b')])
    # in exo_ident because there are different amount of white spaces between
    # catalog and number.

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
    exomercat['mass_pl_err']=ap.table.Column(dtype=object,length=len(exomercat))
    exomercat['mass_pl_rel']=ap.table.Column(dtype=object,length=len(exomercat))
    #transforming mass errors from upper (mass_max) and lower (mass_min) error
    # into instead error (mass_error) as well as relation (mass_pl_rel)
    for i in range(len(exomercat)):
        if type(exomercat['mass_max'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_max'][i]==np.inf:
            if type(exomercat['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_min'][i]==np.inf:
                exomercat['mass_pl_rel'][i]=None
                exomercat['mass_pl_err'][i]=None
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
                            'mass_pl_rel']
    exo_mes_mass_pl.rename_columns(['planet_main_id','mass','mass_url'],
                                    ['main_id','mass_pl_value','mass_pl_ref'])
    #remove masked rows
    exo_mes_mass_pl.remove_rows(exo_mes_mass_pl['mass_pl_value'].mask.nonzero()[0])

    grouped_mes_mass_pl=exo_mes_mass_pl.group_by('main_id')
    #initialize best_mes_mass_pl as a table that hase same columns as exo_mes_mass_pl
    # but no rows
    best_mes_mass_pl=exo_mes_mass_pl['main_id','mass_pl_value','mass_pl_err','mass_pl_rel',
                            'mass_pl_ref'][:0]
    for i in range(len(grouped_mes_mass_pl.groups.keys)):
        #sort by quality
        row=grouped_mes_mass_pl.groups[i][np.where(
                grouped_mes_mass_pl['mass_pl_err'].groups[i]==np.min(
                grouped_mes_mass_pl['mass_pl_err'].groups[i]))][0]
        #take first and add to best_paras
        #which error to take when there are multiples...
        best_mes_mass_pl.add_row([row['main_id'],row['mass_pl_value'],
                              row['mass_pl_err'],row['mass_pl_rel'],row['mass_pl_ref']])
    # changing from string to float, if I do it earlier it raises an error.
    # change in future
    exo_mes_mass_pl['mass_pl_err']=exo_mes_mass_pl['mass_pl_err'].astype(float)
    best_mes_mass_pl['mass_pl_err']=best_mes_mass_pl['mass_pl_err'].astype(float)
    # vstack other multiple measurements tables (currently none)
    best_paras=best_mes_mass_pl

    #-------------exo_h_link---------------
    exo_h_link=exomercat['planet_main_id', 'host_main_id']
    exo_h_link.rename_columns(['planet_main_id','host_main_id'],
                              ['main_id','parent_main_id'])
    exo_h_link['h_link_ref']=[provider_bibcode for j in range(len(exo_h_link))]

    #-------------exo_planet_basic
    exo_planet_basic=best_mes_mass_pl

    #-------------exo_sources---------------
    ref_columns=[['h_link_ref'],['id_ref'],['mass_pl_ref']]
    exo_sources=ap.table.Table()
    tables=[exo_h_link,exo_ident,exo_mes_mass_pl]
    for cat,ref in zip(tables,ref_columns):
        exo_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],exo_sources)

    exo_list_of_tables=[exo_sources,exo_objects,exo_ident,
                        exo_h_link,exo_planet_basic,exo_mes_mass_pl]
    save(exo_list_of_tables,['exo_sources','exo_objects','exo_ident',
                         'exo_h_link','exo_planet_basic','exo_mes_mass_pl'])
    return exo_list_of_tables

def provider_life():
    """
    This function loads the SIMBAD data obtained by the function provider_simbad
    and manipulates it.
    :return life_list_of_tables: List of astropy table containing
        reference data and basic stellar data.
    """
    #---------------define provider--------------------------------------------
    TAP_service="None"
    provider_name='adapted data'
    provider_bibcode='None'#'2022A&A...664A..21Q'
    print(f'Creating {provider_name} tables ...')
    #---------------------star_basic----------------
    #galactic coordinates:  transformed from simbad ircs coordinates using astropy
    [life_star_basic]=load(['sim_star_basic'])
    ircs_coord=ap.coordinates.SkyCoord(\
            ra=life_star_basic['coo_ra'],dec=life_star_basic['coo_dec'],frame='icrs')
    gal_coord=ircs_coord.galactic
    life_star_basic['coo_gal_l']=gal_coord.l.deg*ap.units.degree
    life_star_basic['coo_gal_l'].fill_value=0 #null value
    life_star_basic['coo_gal_b']=gal_coord.b.deg*ap.units.degree
    life_star_basic['coo_gal_b'].fill_value=0 #null value
    # can I do the same transformation with the errors? -> try on some examples and compare to simbad ones
    life_star_basic['coo_gal_err_angle']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_err_maj']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_err_min']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_qual']=['N'
                        for j in range(len(life_star_basic))]
    life_star_basic['main_id']=life_star_basic['main_id'].astype(str)
    # source
    # transformed from simbad ircs coordinates using astropy
    life_star_basic['coo_gal_ref']=ap.table.MaskedColumn(dtype=object,length=len(life_star_basic),
                                    mask=[True for j in range(len(life_star_basic))])
    life_star_basic['coo_gal_ref'][np.where(life_star_basic['coo_gal_l']!=0)]='LIFE' 
    life_star_basic['coo_gal_ref']=life_star_basic['coo_gal_ref'].astype(str)
    life_star_basic=life_star_basic['main_id','coo_gal_l','coo_gal_b','coo_gal_err_angle',
                                   'coo_gal_err_maj','coo_gal_err_min','coo_gal_qual',
                                   'coo_gal_ref']
    #-----------measurement tables -----------------
    #applying model from E. E. Mamajek on SIMBAD spectral type
    
    def replace_value(cat,column,value,replace_by):
        cat[column][np.where(cat[column]==value)]= \
                [replace_by for i in range(
            len(cat[column][np.where(cat[column]==value)]))]
        return cat
             
    def realspectype(cat):
        """
        Removes rows of cat where elements in column named 'sim_sptype' are
        either '', 'nan' or start with an other letter than the main sequence
        spectral type classification.
        :param cat: astropy table containing 'sim_sptype' column
        :return cat: astropy table, param cat with undesired rows removed
        """
        index=[]
        for j in range(len(cat['sptype_value'])):
            if cat['sptype_value'][j] in ['','nan']:
                index.append(j)
            elif cat['sptype_value'][j][0] not in ['O','B','A','F','G','K','M']:
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
        EEM_table=ap.io.ascii.read("data/updatedEEM_dwarf_UBVIJHK_colors_Teff.csv")['SpT','Teff','R_Rsun','Msun']
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
        cat=match_sptype(cat,mp,sptypestring='sptype_value')
        cat.remove_rows([np.where(cat['mod_Teff'].mask==True)])
        cat.remove_rows([np.where(np.isnan(cat['mod_Teff']))])
        cat=ap.table.unique(cat, keys='main_id')
        return cat

    [sim_objects]=load(['sim_objects'],stringtoobjects=False)
    stars=sim_objects[np.where(sim_objects['type']=='st')]
    cat=ap.table.join(stars,sim_star_basic)
    cat=spec(cat['main_id','sptype_value'])
    #if I take only st objects from sim_star_basic I don't loose objects during realspectype
    life_mes_teff_st=cat['main_id','mod_Teff']
    life_mes_teff_st.rename_column('mod_Teff','teff_st_value')
    life_mes_teff_st['teff_st_qual']=['C' for i in range(len(life_mes_teff_st))]
    life_mes_teff_st['teff_st_ref']=['LIFE' for i in range(len(life_mes_teff_st))]
    
    life_mes_radius_st=cat['main_id','mod_R']
    life_mes_radius_st.rename_column('mod_R','radius_st_value')
    life_mes_radius_st['radius_st_qual']=['C' for i in range(len(life_mes_radius_st))]
    life_mes_radius_st['radius_st_ref']=['LIFE' for i in range(len(life_mes_radius_st))]
    
    life_mes_mass_st=cat['main_id','mod_M']
    life_mes_mass_st.rename_column('mod_M','mass_st_value')
    life_mes_mass_st['mass_st_qual']=['C' for i in range(len(life_mes_mass_st))]
    life_mes_mass_st['mass_st_ref']=['LIFE' for i in range(len(life_mes_mass_st))]
    
#ok so next step is to implement this into life provider and create tables mes_teff_st, mes_rad_st, mes_mass_st
#step after that one is to merge in building function different mes tables and get best para from it
    
    #-----------------sources table----------------------
    life_sources=ap.table.Table()
    tables=[life_star_basic,life_mes_teff_st,life_mes_radius_st,life_mes_mass_st]
    ref_columns=[['coo_gal_ref'],['teff_st_ref'],['radius_st_ref'],['mass_st_ref']]
    for cat,ref in zip(tables,ref_columns):
        life_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],life_sources)
    
    life_list_of_tables=[life_sources,life_star_basic,life_mes_teff_st,
                         life_mes_radius_st,life_mes_mass_st]
    save(life_list_of_tables,['life_sources','life_star_basic','life_mes_teff_st',
                         'life_mes_radius_st','life_mes_mass_st'])
    return life_list_of_tables

def provider_gaia(temp=True):
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
    TAP_service="https://gea.esac.esa.int/tap-server/tap" #might not be all of gaia
    provider_name='Gaia'
    provider_bibcode='2016A&A...595A...1G'
    #query
    adql_query="""
    SELECT s.source_id ,p.mass_flame, p.radius_flame,
        p.teff_gspphot, p.teff_gspspec 
    FROM gaiadr3.gaia_source as s
        JOIN gaiadr3.astrophysical_parameters as p ON s.source_id=p.source_id
    WHERE parallax >="""+str(plx_in_mas_cut)
    
    print(f'Creating {provider_name} tables ...')
    if temp:
        gaia=ap.io.ascii.read("data/gaia26mai23.csv")
        gaia=stringtoobject(gaia,3000)

    else:
        gaia=query(TAP_service,adql_query) 
        
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
    #gaia_objects['type']=['None' for j in range(len(gaia_objects))]

    
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
    tables=[gaia_ident,gaia_mes_teff_st,gaia_mes_radius_st,gaia_mes_mass_st]
    ref_columns=[['id_ref'],['teff_st_ref'],['radius_st_ref'],['mass_st_ref']]
    for cat,ref in zip(tables,ref_columns):
        gaia_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],gaia_sources)
        
    gaia_list_of_tables=[gaia_sources,gaia_objects,gaia_ident,
                         gaia_mes_teff_st,
                        gaia_mes_radius_st,gaia_mes_mass_st]
    save(gaia_list_of_tables,['gaia_sources','gaia_objects','gaia_ident',
                              'gaia_mes_teff_st',
                             'gaia_mes_radius_st','gaia_mes_mass_st'])
    return gaia_list_of_tables

#------------------------provider combining-----------------
def building(providers,column_names):
    """
    This function builds from the input parameters the tables
    for the LIFE database.
    :param sim: List of astropy table containing simbad data.
    :param gk: List of astropy table containing grant kennedy data.
    :param exo: List of astropy table containing exomercat data.
    :return cat: List of astropy table containing
        reference data, object data, identifier data, object to object
        relation data, basic stellar data, basic planetary data, basic disk
        data, distance measurement data and mass measurement data.
    """
    #creates empty tables as needed for final database ingestion
    init=initialize_database_tables()
    n_tables=len(init)

    cat=[ap.table.Table() for i in range(n_tables)]

    #for the sources and objects joins tables from different providers
    
    print(f'Building {column_names[0]} table ...')
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
                if type(cat['type_2'][i])==np.ma.core.MaskedConstant:
                    cat['type'][i]=cat['type_1'][i]
                else:
                    cat['type'][i]=cat['type_2'][i]
            cat.remove_columns(['type_1','type_2'])
        return cat
    
    print(f'Building {column_names[1]} table ...')
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
                cat=nullvalues(cat,para+'_ref','')
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
                            cat[f'{para}_source_idref'][i]=0
        return cat
    paras=[['id'],['h_link'],['coo','plx','dist','coo_gal'],
           ['mass_pl'],['rad'],['dist'],['mass_pl'],['teff_st'],
          ['radius_st'],['mass_st']]
    
    prov_ref=['SIMBAD','priv. comm.','Exo-MerCat','adapted data','Gaia']
    
    for i in range(2,n_tables):
        print(f'Building {column_names[i]} table ...')
        
        for j in range(len(providers)):
            if len(providers[j][i])>0:
                #replacing ref with source_idref columns
                #getting source_idref to each ref
                providers[j][i]=match(providers[j][i],cat[0],paras[i-2],prov_ref[j])
            if len(cat[i])>0:
                #joining data from different providers
                if len(providers[j][i])>0:
                    cat[i]=ap.table.join(cat[i],providers[j][i],join_type='outer')
            else:
                cat[i]=providers[j][i]
        
        #I do this to get those columns that are empty in the data
        cat[i]=ap.table.vstack([cat[i],init[i]])

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
                mes_table=mes_table['main_id',para+'_value',para+'_err',para+'_qual',para+'_source_idref']
                best_para=mes_table['main_id',para+'_value',para+'_err',para+'_qual',para+'_source_idref'][:0].copy()
                #group mes_table by object (=main_id)
                grouped_mes_table=mes_table.group_by('main_id')
                #take highest quality
                for j in range(len(grouped_mes_table.groups.keys)):#go through all objects
                    for qual in ['A','B','C','D','E']:
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

            if i==3:#--------------------h_link--------------------------
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
                cat[i]['membership'].fill_value=-1
                cat[i]['membership']=cat[i]['membership'].filled()

            #for all the other tables add object_idref
            else:
                #first remove the object_idref we got from empty initialization
                #though I would prefer a more elegant way to do this
                cat[i].remove_column('object_idref') 
                cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],
                                     join_type='left')
                cat[i].rename_column('object_id','object_idref')
            if i==4:#--------------------star_basic--------------------------
                #choosing all objects with type star or system
                #I am just adding main_id because I have not found out how
                # to do join with just one column table
                stars=cat[1]['object_id','main_id'][np.where(
                                cat[1]['type']=='st')]
                systems=cat[1]['object_id','main_id'][np.where(
                                cat[1]['type']=='sys')]
                temp=ap.table.vstack([stars,systems])
                temp.rename_columns(['object_id','main_id'],
                                    ['object_idref','temp'])
                cat[i]=ap.table.join(cat[i],temp,join_type='outer',
                                     keys='object_idref')
                cat[i].remove_column('temp')
            if i==5:#--------------------planet_basic--------------------------
                temp=cat[1]['object_id','main_id'][np.where(
                                cat[1]['type']=='pl')]
                temp.rename_columns(['object_id','main_id'],
                                    ['object_idref','temp'])
                cat[i]=ap.table.join(cat[i],temp,keys='object_idref',
                                     join_type='outer')
                cat[i].remove_column('temp')
            if i==9:##-------mes_teff_st
                teff_st_best_para=best_para('teff_st',cat[i])
                cat[4].remove_columns(['teff_st_value','teff_st_err',
                                       'teff_st_qual','teff_st_source_idref',
                                       'teff_st_ref'])
                cat[4]=ap.table.join(cat[4],teff_st_best_para)
            if i==10:##-------mes_radius_st
                radius_st_best_para=best_para('radius_st',cat[i])
                cat[4].remove_columns(['radius_st_value','radius_st_err',
                                       'radius_st_qual','radius_st_source_idref',
                                       'radius_st_ref'])
                cat[4]=ap.table.join(cat[4],radius_st_best_para)
            if i==11:##-------mes_mass_st
                mass_st_best_para=best_para('mass_st',cat[i])
                cat[4].remove_columns(['mass_st_value','mass_st_err',
                                       'mass_st_qual','mass_st_source_idref',
                                       'mass_st_ref'])
                cat[4]=ap.table.join(cat[4],mass_st_best_para)
                             
        else:
            print('error: empty table')

    save(cat,column_names)
    return cat

###############################################################################
#-------------------------Main code--------------------------------------------
###############################################################################

#------------------------initialize empty database tables----------------------
initialized_tables=initialize_database_tables()

column_names=['sources','objects','ident','h_link','star_basic',
              'planet_basic','disk_basic','mes_dist','mes_mass_pl',
              'mes_teff_st','mes_radius_st','mes_mass_st']
save(initialized_tables,column_names)

#------------------------obtain data from external sources---------------------
sim_sources,sim_objects,sim_ident,sim_h_link \
            ,sim_star_basic,sim_mes_dist=provider_simbad()

gk_sources,gk_objects, gk_ident, gk_h_link,gk_disk_basic=provider_gk()

exo_sources,exo_objects,exo_ident,exo_h_link \
            ,exo_planet_basic,exo_mes_mass_pl=provider_exo()

life_sources,life_star_basic,life_mes_teff_st,life_mes_radius_st,life_mes_mass_st=provider_life()

gaia_sources,gaia_objects,gaia_ident,gaia_mes_teff_st,gaia_mes_radius_st,gaia_mes_mass_st=provider_gaia()

#------------------------construct the database tables-------------------------
empty=ap.table.Table()
sim=[sim_sources,sim_objects,sim_ident,sim_h_link,sim_star_basic,empty[:],empty[:],
     sim_mes_dist,empty[:],empty[:],empty[:],empty[:]]
gk=[gk_sources,gk_objects, gk_ident, gk_h_link,empty[:],empty[:],gk_disk_basic,
    empty[:],empty[:],empty[:],empty[:],empty[:]]
exo=[exo_sources,exo_objects, exo_ident, exo_h_link,empty[:],exo_planet_basic,empty[:],
     empty[:],exo_mes_mass_pl,empty[:],empty[:],empty[:]]
life=[life_sources,empty[:],empty[:],empty[:],life_star_basic,empty[:],empty[:],
      empty[:],empty[:],empty[:],empty[:],empty[:]]

gaia=[gaia_sources,gaia_objects,gaia_ident,empty[:],empty[:],empty[:],empty[:],
      empty[:],empty[:],gaia_mes_teff_st,gaia_mes_radius_st,gaia_mes_mass_st]

sources,objects,ident,h_link,star_basic,planet_basic,disk_basic,mes_dist, \
            mes_mass_pl,mes_teff_st,mes_radius_st,mes_mass_st=building([sim,gk,exo,life,gaia],column_names)