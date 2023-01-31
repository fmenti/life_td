"""
Authors Note
This code was written by Franziska Menti in 2022. 
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

def load(paths):
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
    for cat in cats:
        cat=stringtoobject(cat,3000)
    return cats

#-------------------initialization function------------------------------------
def initialize_database_tables():
    """
    This function initializes the database tables with column name and data 
    type specified but no actual data in them.
    :return list_of_tables: List of astropy tables in the order sources, 
                objects, ident (identifiers), h_link (relation between 
                objects),star_basic,planet_basic, disk_basic,
                mesDist (distance measurements) and 
                mesMass (mass measurements). 
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
               'plx_value','plx_err','plx_qual','plx_source_idref',
               'plx_ref',
               'dist_value','dist_err','dist_qual','dist_source_idref',
               'dist_ref'],
        dtype=[int,float,float,float,
               float,float,object,
               int,object,
               float,float,object,int,
               object,
               float,float,object,int,
               object])
    list_of_tables.append(star_basic)
    
    planet_basic=ap.table.Table(
        #object idref, mass value, mass error, mass realtion (min, max, equal),
        #mass quality, source idref of mass parameter, mass reference
        names=['object_idref','mass_val','mass_err','mass_rel','mass_qual',
               'mass_source_idref','mass_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(planet_basic)
    
    disk_basic=ap.table.Table(
        #object idref, black body radius value, bbr error, 
        #bbr relation (min, max, equal), bbr quality,...
        names=['object_idref','rad_value','rad_err','rad_rel','rad_qual',
               'rad_source_iderf','rad_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(disk_basic)
    
    mesDist=ap.table.Table(
        names=['object_idref','dist_value','dist_err','dist_qual',
               'dist_source_idref','dist_ref'],
        dtype=[int,float,float,object,
               int,object])
    list_of_tables.append(mesDist)
    
    mesMass=ap.table.Table(
        names=['object_idref','mass_val','mass_err','mass_rel','mass_qual',
               'mass_source_idref','mass_ref'],
        dtype=[int,float,float,object,object,int,object])
    list_of_tables.append(mesMass)
    
    #save all tables
    save(list_of_tables,
         ['empty_sources','empty_objects','empty_ident','empty_h_link',
         'empty_star_basic','empty_planet_basic','empty_disk_basic',
         'empty_mesDist','empty_mesMass'])
    return list_of_tables

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
        b.coo_bibcode AS coo_ref, b.coo_qual,
        b.plx_err, b.plx_value, b.plx_bibcode AS plx_ref,b.plx_qual,
        h_link.membership, h_link.parent AS parent_oid, 
        h_link.link_bibcode AS h_link_ref, a.otypes,ids.ids
        """#which parameters to query from simbad and what alias to give them
    adql_query=[
        select+"""
        FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
        WHERE b.plx_value >=50."""]
    #creating one table out of parameters from multiple ones and
    #keeping only objects with parallax bigger than 50mas
    
    upload_query=[
        #query for systems without parallax data but 
        #children (in TAP_UPLOAD.t1 table) with parallax bigger than 50mas
        select+"""
        FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    JOIN TAP_UPLOAD.t1 ON b.oid=t1.parent_oid
        WHERE (b.plx_value IS NULL) AND (otype='**..')""",
        #query for planets without parallax data but 
        #host star (in TAP_UPLOAD.t1 table) with parallax bigger than 50mas
        select+"""
        FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    JOIN TAP_UPLOAD.t1 ON b.oid=t1.oid
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
    ref_columns=[['coo_ref','plx_ref'],['h_link_ref'],['dist_ref'],['id_ref']]
    #------------------querrying-----------------------------------------------
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
    
    #-------------------creating output table sim_mesDist---------------------
    sim_mesDist=query(TAP_service,upload_query[2],[stars[:].copy()])
    sim_mesDist=fetch_main_id(sim_mesDist)
    sim_mesDist['dist_err']=np.maximum(sim_mesDist['plus_err'],
                                       -sim_mesDist['minus_err'])
    sim_mesDist.remove_rows(sim_mesDist['dist_err'].mask.nonzero()[0])
    #group by oid
    grouped_mesDist=sim_mesDist.group_by('main_id')
    best_mesDist=sim_mesDist['main_id','dist_value','plus_err',
                             'dist_qual','dist_ref'][:0]
    best_mesDist.rename_column('plus_err','dist_err')
    for i in range(len(grouped_mesDist.groups.keys)):
        #sort by quality
        row=grouped_mesDist.groups[i][np.where(
                    grouped_mesDist['dist_prec'].groups[i]
                    ==np.max(grouped_mesDist['dist_prec'].groups[i]))][0]
        #take first and add to best_paras
        #which error to take when there are multiples...
        best_mesDist.add_row([row['main_id'],row['dist_value'], 
                              row['dist_err'],row['dist_qual'], 
                              row['dist_ref']])
    #join with other multimes thingis
    best_paras=best_mesDist
    # TBD when more multi measurement tables are implemented: vstack them here
    sim_mesDist=sim_mesDist['main_id','dist_value','dist_err',
                            'dist_qual','dist_ref']
    
    #--------------------creating helper table sim_stars-----------------------
    #add best para from multiple measurements tables
    stars=ap.table.join(stars,best_paras,keys='main_id',join_type='left')
    #change null value of plx_qual from '' to 'N'
    stars['plx_qual'][np.where(stars['plx_qual']=='')]='N'
    stars['dist_qual'][np.where(stars['dist_qual']=='')]='N'
    stars['dist_qual'][np.where(stars['dist_qual']==':')]='N'
    
    #--------------creating output table sim_h_link ---------------------------
    sim_h_link=simbad['main_id','parent_oid','h_link_ref','membership']
    # if you want to exclude objects with lower membership probability 
    # use this line instead:
    # sim_h_link=simbad['main_id','parent_oid','h_link_ref',
    #                   'membership'][np.where(simbad['membership']>50)]
    # consequence is that you loose objects with no membership value given (~) 
    # e.g. alf cen system
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
    tables=[stars, sim_h_link, sim_mesDist,sim_ident]
    for cat,ref in zip(tables,ref_columns):
        sim_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],sim_sources)
    simbad_list_of_tables=[sim_sources,sim_objects, sim_ident, sim_h_link]
    #------------------------creating output table sim_star_basic--------------
    sim_star_basic=stars['main_id','coo_ra','coo_dec','coo_err_angle',
                         'coo_err_maj','coo_err_min','coo_qual','coo_ref',
                         'plx_value','plx_err','plx_qual','plx_ref',
                         'dist_value','dist_err','dist_qual','dist_ref']
    simbad_list_of_tables.extend([sim_star_basic,sim_mesDist])
    save(simbad_list_of_tables,['sim_sources','sim_objects','sim_ident',
                                'sim_h_link','sim_star_basic','sim_mesDist'])
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
    #sorting out everything with plx_value <50mas (corresponding to >20pc)
    gk_disks=gk_disks[np.where(gk_disks['plx_value']>50.)]
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

    #-------------------exo_mesMass---------------------
    #initialize columns exomercat['mass_rel'] and exomercat['mass_err']
    exomercat['mass_err']=ap.table.Column(dtype=object,length=len(exomercat))
    exomercat['mass_rel']=ap.table.Column(dtype=object,length=len(exomercat))
    #transforming mass errors from upper (mass_max) and lower (mass_min) error 
    # into instead error (mass_error) as well as relation (mass_rel)
    for i in range(len(exomercat)):
        if type(exomercat['mass_max'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_max'][i]==np.inf:
            if type(exomercat['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_min'][i]==np.inf:
                exomercat['mass_rel'][i]=None
                exomercat['mass_err'][i]=None
            else:
                exomercat['mass_rel'][i]='<'
                exomercat['mass_err'][i]=exomercat['mass_min'][i]
        else:
            if type(exomercat['mass_min'][i])==np.ma.core.MaskedConstant or \
                  exomercat['mass_min'][i]==np.inf:
                exomercat['mass_rel'][i]='>'
                exomercat['mass_err'][i]=exomercat['mass_max'][i]
            else:
                exomercat['mass_rel'][i]='='
                exomercat['mass_err'][i]=max(exomercat['mass_max'][i],
                                        exomercat['mass_min'][i])
    exo_mesMass=exomercat['planet_main_id','mass','mass_err','mass_url',
                            'mass_rel']
    exo_mesMass.rename_columns(['planet_main_id','mass','mass_url'],
                                    ['main_id','mass_val','mass_ref'])
    #remove masked rows
    exo_mesMass.remove_rows(exo_mesMass['mass_val'].mask.nonzero()[0])
    
    grouped_mesMass=exo_mesMass.group_by('main_id')
    #initialize best_mesMass as a table that hase same columns as exo_mesMass 
    # but no rows
    best_mesMass=exo_mesMass['main_id','mass_val','mass_err','mass_rel',
                            'mass_ref'][:0]
    for i in range(len(grouped_mesMass.groups.keys)):
        #sort by quality
        row=grouped_mesMass.groups[i][np.where(
                grouped_mesMass['mass_err'].groups[i]==np.min(
                grouped_mesMass['mass_err'].groups[i]))][0]
        #take first and add to best_paras
        #which error to take when there are multiples...
        best_mesMass.add_row([row['main_id'],row['mass_val'],
                              row['mass_err'],row['mass_rel'],row['mass_ref']])
    # changing from string to float, if I do it earlier it raises an error. 
    # change in future
    exo_mesMass['mass_err']=exo_mesMass['mass_err'].astype(float)
    best_mesMass['mass_err']=best_mesMass['mass_err'].astype(float)
    # vstack other multiple measurements tables (currently none) 
    best_paras=best_mesMass
    
    #-------------exo_h_link---------------
    exo_h_link=exomercat['planet_main_id', 'host_main_id']
    exo_h_link.rename_columns(['planet_main_id','host_main_id'],
                              ['main_id','parent_main_id'])
    exo_h_link['h_link_ref']=[provider_bibcode for j in range(len(exo_h_link))]

    #-------------exo_planet_basic
    exo_planet_basic=best_mesMass
    
    #-------------exo_sources---------------
    ref_columns=[['mass_url'],['h_link_ref'],['id_ref']]
    exo_sources=ap.table.Table()
    tables=[exomercat, exo_h_link,exo_ident]
    for cat,ref in zip(tables,ref_columns):
        exo_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],exo_sources) 
        
    exo_list_of_tables=[exo_sources,exo_objects,exo_ident,
                        exo_h_link,exo_planet_basic,exo_mesMass]
    save(exo_list_of_tables,['exo_sources','exo_objects','exo_ident',
                         'exo_h_link','exo_planet_basic','exo_mesMass'])
    return exo_list_of_tables

#------------------------provider combining-----------------
def building(sim,gk,exo):
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
    
    #initializes 8 table objects
    #corresponding to 'sources','objects','ident','h_link','star_basic',
    #                 'planet_basic','disk_basic','mesDist', 'mesMass'
    cat=[ap.table.Table() for i in range(9)]
    
    #for the sources and objects stacks tables from different providers 
    # keeping only unique values then create identifiers for those tables
    
    cat[0]=ap.table.vstack([init[0],sim[0]])
    cat[0]=ap.table.vstack([cat[0],gk[0]])
    cat[0]=ap.table.vstack([cat[0],exo[0]])
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
            if type(cat[column])==ap.table.column.MaskedColumn:
                cat[column].fill_value=''
                cat[column]=cat[column].filled()
        for i in range(len(cat)):
            #splitting object into list of elements
            ids1=cat[column_ids1][i].split('|')
            ids2=cat[column_ids2][i].split('|')
            if ids2==['']:
                cat['ids'][i]=cat[column_ids1][i]
            if ids1==['']:
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
        return cat
    #removed vstack with init to not have object_id as is empty anyways
    #cat[1]=ap.table.vstack([init[1],sim[1]]) 
    cat[1]=sim[1]
    cat[1]=ap.table.join(cat[1],gk[1],keys='main_id',join_type='outer')

    # could not solve the warning message about type merging in type part but 
    # the results seem to be okay

    def objectmerging(cat):
        """
        This function merges the data from the same physical object obtained 
        from different providers into one entry. 
        :param cat: Astropy table containing multiple entries for the same 
            physical objects due to data from different providers.
        :return cat: Astropy table with unique object entries.
        """
        cat=idsjoin(cat,'ids_1','ids_2')
        #merging types
        #initializing column
        cat['type']=ap.table.Column(dtype=object, length=len(cat))
        cat['type_1']=cat['type_1'].astype(object)
        cat['type_2']=cat['type_2'].astype(object)
        for i in range(len(cat)):
            if type(cat['type_2'][i])==np.ma.core.MaskedConstant:
                cat['type'][i]=cat['type_1'][i]
            else:
                cat['type'][i]=cat['type_2'][i]
        cat.remove_columns(['ids_1','ids_2','type_1','type_2'])
        return cat

    cat[1]=objectmerging(cat[1])
    cat[1]=ap.table.join(cat[1],exo[1],keys='main_id',join_type='outer')
    cat[1]=objectmerging(cat[1])
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
                if type(cat[para+'_ref'])==ap.table.column.MaskedColumn:
                    cat[para+'_ref'].fill_value=''
                    cat[para+'_ref']=cat[para+'_ref'].filled()
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
    paras=[['id'],['h_link'],['coo','plx','dist'],
           ['mass'],['rad'],['dist'],['mass']]
    
    for i in range(2,9):
        #replacing ref with source_idref columns
        #getting source_idref to each ref
        sim[i]=match(sim[i],cat[0],paras[i-2],'SIMBAD')
        gk[i]=match(gk[i],cat[0],paras[i-2],'priv. comm.')
        exo[i]=match(exo[i],cat[0],paras[i-2],'Exo-MerCat') 
        #joining data from different providers
        #I do this to get those columns that are empty in the data
        cat[i]=ap.table.vstack([init[i],sim[i]])
        cat[i]=ap.table.vstack([cat[i],gk[i]])
        cat[i]=ap.table.vstack([cat[i],exo[i]])
        
        #if resulting catalog not empty
        if len(cat[i])>0:
            #only keeping unique entries
            cat[i]=ap.table.unique(cat[i],silent=True)
            #i==2#--------------------ident--------------------------
            #no special action required

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
        else:
            print('error: empty table')
            
    save(cat,['sources','objects','ident','h_link','star_basic',
              'planet_basic','disk_basic','mesDist','mesMass'])
    return cat
###############################################################################
#-------------------------Main code--------------------------------------------
###############################################################################

#------------------------initialize empty database tables----------------------
initialized_tables=initialize_database_tables()

save(initialized_tables,['sources','objects','ident','h_link','star_basic',
                         'planet_basic','disk_basic','mesDist','mesMass'])

#------------------------obtain data from external sources---------------------
sim_sources,sim_objects,sim_ident,sim_h_link \
            ,sim_star_basic,sim_mesDist=provider_simbad()
gk_sources,gk_objects, gk_ident, gk_h_link,gk_disk_basic=provider_gk()
exo_sources,exo_objects,exo_ident,exo_h_link \
            ,exo_planet_basic,exo_mesMass=provider_exo()

#------------------------construct the database tables-------------------------
empty=ap.table.Table()
sim=[sim_sources,sim_objects,sim_ident,sim_h_link,sim_star_basic,
     empty[:],empty[:],sim_mesDist,empty[:]]
gk=[gk_sources,gk_objects, gk_ident, gk_h_link,empty[:],empty[:],gk_disk_basic,empty[:],empty[:]]
exo=[exo_sources,exo_objects, exo_ident, exo_h_link,empty[:],exo_planet_basic,empty[:],empty[:],exo_mesMass]

sources,objects,ident,h_link,star_basic,planet_basic \
            ,disk_basic,mesDist,mesMass=building(sim,gk,exo)
