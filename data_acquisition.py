import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #table handling

#-------------------global helper functions-----------------------
def save(cats,paths):
    """
    Saves cat in the location path as .xml file.
    :param cat: astropy table to be saved
    :param path: path to where to save the cat table
    """
    for cat,path in zip(cats,paths):
        for i in cat.colnames:
            if cat[i].dtype == object:
                cat[i] = cat[i].astype(str)
        ap.io.votable.writeto(ap.io.votable.from_table(cat), f'data/{path}.xml')
    return
def load(paths):
    """
    Loads tables saved in .xml format at locations specified in paths.
    :param path: path to file location.
    :return cats: list of loaded tables
    """
    cats=[]
    for path in paths:
        cats.append(ap.io.votable.parse_single_table(f'data/{path}.xml').to_table())
    stringtypes=[np.dtype(f'<U{j}') for j in range(1,3000)]
    for cat in cats:
        #making sure that objecttype is string which can change number of characters
        for i in cat.colnames:
            if cat[i].dtype in stringtypes:#no nicer way found so far, ==str does not work.
                cat[i] = cat[i].astype(object)
    return cats

#-------------------initialization function----------------------------
def initialize_database_tables():
    """
    This function initializes the database tables with no data in them.
    """
    #initialize tables with no data but column names and data type specified.
    objects=ap.table.Table(
        names=['object_id','type','ids','main_id'],
        dtype=[int,object,object,object])
    ident=ap.table.Table(
        names=['object_idref','id','id_source_idref','id_ref'],
        dtype=[int,object,int,object])
    h_link=ap.table.Table(
        names=['child_object_idref','parent_object_idref',
               'h_link_source_idref','h_link_ref','membership'],
        dtype=[int,int,int,object,int])
    star_basic=ap.table.Table(
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
    planet_basic=ap.table.Table(
        names=['object_idref','mass_val','mass_err','mass_rel','mass_qual',
               'mass_source_idref','mass_ref'],
        dtype=[int,float,float,object,object,int,object])
    disk_basic=ap.table.Table(
        names=['object_idref','rad_value','rad_err','rad_rel','rad_qual',
               'rad_source_iderf','rad_ref'],
        dtype=[int,float,float,object,object,int,object])
    sources=ap.table.Table(
        names=['ref','provider_name','provider_url','provider_bibcode','source_id'],
        dtype=[object,object,object,object,int])
    mesDist=ap.table.Table(
        names=['object_idref','dist_value','dist_err','dist_qual',
               'dist_source_idref','dist_ref'],
        dtype=[int,float,float,object,
               int,object])
    mesMass=ap.table.Table(
        names=['object_idref','mass_val','mass_err','mass_rel','mass_qual',
               'mass_source_idref','mass_ref'],
        dtype=[int,float,float,object,object,int,object])
    
    #save all tables
    save([sources,objects,ident,h_link,star_basic,planet_basic,disk_basic,mesDist,mesMass],
         ['empty_sources','empty_objects','empty_ident','empty_h_link','empty_star_basic',
          'empty_planet_basic','empty_disk_basic','empty_mesDist','empty_mesMass'])
    return [sources,objects,ident,h_link,star_basic,planet_basic,disk_basic,mesDist,mesMass]

#------------------------------provider helper functions--------------------------
def query(link,query,catalogs=[]):
    """Performs a query via TAP on the service given in the link parameter. 
    If catalogs is specified those tables are uploaded to the service. 
    :param link: service URL
    :param query: query in ADQL 
    :param catalogs: list of tables to be uploaded to the service
    :return result.to_table(): astropy table containing the result of the query"""
    service = vo.dal.TAPService(link)
    if catalogs==[]:
        result=service.run_async(query.format(**locals()), maxrec=160000)
    else:
        tables={}
        for i in range(len(catalogs)):
            tables.update({f"t{i+1}":catalogs[i]})
        result = service.run_async(query,uploads=tables,timeout=None, maxrec=160000)
    cat=result.to_table()
    #for i in cat.colnames:
     #   if cat[i].dtype == object:
      #      cat[i] = cat[i].astype(str)
    return cat

def sources_table(cat,ref_columns,provider,old_sources=ap.table.Table()): #put this into source function
    """
    This function collects all the references in the ref_colums, 
    keeps only unique entries and adds the columns provider_name, 
    provider_url and provider_bibcode.
    :param provider: List containing name, url and bibcode of provider.
    """
    cat_sources=ap.table.Table()
    cat_reflist=[]
    for k in range(len(ref_columns)):
        if type(cat[ref_columns[k]])==ap.table.column.MaskedColumn:
            cat_reflist.extend(cat[ref_columns[k]][np.where(cat[ref_columns[k]].mask==False)])
        else:
            cat_reflist.extend(cat[ref_columns[k]])
    cat_sources['ref']=cat_reflist
    cat_sources=ap.table.unique(cat_sources)#keeps only unique values
    #attaches service information
    cat_sources['provider_name']=[provider[0]]*len(cat_sources)
    cat_sources['provider_url']=[provider[1]]*len(cat_sources)
    cat_sources['provider_bibcode']=[provider[2]]*len(cat_sources)
    sources=ap.table.vstack([old_sources,cat_sources])
    sources=ap.table.unique(sources)
    return sources

def fetch_main_id(cat,colname='oid',name='main_id',oid=True):
    """
    Joins main_id from simbad to the column colname. Returns the whole
    table cat but without any rows where no simbad main_id was found.
    """
    print('tbd option to match on position instead of main_id or oid')
    TAP_service="http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    if oid:
        main_id_query='SELECT b.main_id AS '+name+""",t1.*
                    FROM basic AS b
                    JOIN TAP_UPLOAD.t1 ON b.oid=t1."""+colname
    else:
        main_id_query='SELECT b.main_id AS '+name+""",t1.*
                    FROM basic AS b
                    JOIN ident ON ident.oidref=b.oid
                        JOIN TAP_UPLOAD.t1 ON ident.id=t1."""+colname
    cat=query(TAP_service,main_id_query,[cat])
    return cat

def stringtoobject(cat,number=100):
    """
    This function changes from string to object format. 
    The later has the advantace of allowing strings of varying length.
    """
    stringtypes=[np.dtype(f'<U{j}') for j in range(1,number)]
    for i in cat.colnames:
        if cat[i].dtype in stringtypes:
            cat[i] = cat[i].astype(object)
    return cat
#-----------------------------provider data ingestion------------------
def provider_simbad():
    """
    This function obtains the simbad data and arranges it in a way 
    easy to ingest into the database.
    """
    #---------------define provider------------------------
    TAP_service="http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    provider_name='SIMBAD'
    provider_bibcode='2000A&AS..143....9W'
    #---------------define queries---------------------------
    select="""SELECT b.main_id,b.ra AS coo_ra,b.dec AS coo_dec,
        b.coo_err_angle, b.coo_err_maj, b.coo_err_min,b.oid, 
        b.coo_bibcode AS coo_ref, b.coo_qual,
        b.plx_err, b.plx_value, b.plx_bibcode AS plx_ref,b.plx_qual,
        h_link.membership, h_link.parent AS parent_oid, 
        h_link.link_bibcode AS h_link_ref, a.otypes,ids.ids
        """
    adql_query=[
        select+"""
        FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
        WHERE b.plx_value >=50."""]
    upload_query=[select+"""
        FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    JOIN TAP_UPLOAD.t1 ON b.oid=t1.parent_oid
        WHERE (b.plx_value IS NULL) AND (otype='**..')""",
        select+"""
        FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    JOIN TAP_UPLOAD.t1 ON b.oid=t1.oid
        WHERE (b.plx_value IS NULL) AND (otype='Pl..')""",
        """SELECT oid, dist AS dist_value, plus_err, qual AS dist_qual,
        bibcode AS dist_ref,minus_err,dist_prec
        FROM mesDistance
        JOIN TAP_UPLOAD.t1 ON oidref=t1.oid""",
        """SELECT id, t1.*
        FROM ident 
        JOIN TAP_UPLOAD.t1 ON oidref=t1.oid"""]
    #references in columns
    ref_columns=[['coo_ref','plx_ref'],['h_link_ref'],['dist_ref'],['id_ref']]
    
    simbad=query(TAP_service,adql_query[0])
    #adds parent and children objects with no parallax value
    parents_without_plx=query(TAP_service,upload_query[0],[simbad])
    children_without_plx=query(TAP_service,upload_query[1],[simbad])
    #adding of no_parallax objects to rest of simbad query objects
    simbad=ap.table.vstack([simbad,parents_without_plx])
    simbad=ap.table.vstack([simbad,children_without_plx])
    #sorting from object type into star, system and planet type
    simbad['type']=['None' for i in range(len(simbad))]
    simbad['multiple']=[False for i in range(len(simbad))]
    to_remove_list=[]
    for i in range(len(simbad)):
        if "Pl" in simbad['otypes'][i]:
            simbad['type'][i]='pl'
        elif "*" in simbad['otypes'][i]:
            if "**" in simbad['otypes'][i]:
                simbad['type'][i]='sy'
                simbad['multiple'][i]=True
            else:
                simbad['type'][i]='st'
        else:
            print('Removed object because type neither Pl,* or **:',simbad['otypes'][i])
            to_remove_list.append(i)
    if to_remove_list!=[]:
        simbad.remove_rows(to_remove_list)
            
    temp_stars=simbad[np.where(simbad['type']!='pl')]
    #removing double objects (in there due to multiple parents)
    stars=ap.table.Table(ap.table.unique(temp_stars,keys='main_id'),copy=True)
    
    #-----------------sim_ident
    #wait I actually need the oid here, or I need to rewrite the query -> use simbad instead of sim_objects.
    #retlased column type with main_id
    sim_ident=query(TAP_service,upload_query[3],[simbad['oid','main_id'][:].copy()])
    sim_ident['id_ref']=[provider_bibcode for j in range(len(sim_ident))]
    sim_ident.remove_column('oid')        
    
    #-------------------sim_mesDist---------------------
    ######wait, isn't this problematic? do I need another table as best_para. at the moment it works    
    sim_mesDist=query(TAP_service,upload_query[2],[stars[:].copy()])
    print(sim_mesDist['dist_value'].fill_value)
    print('insert here something to fill in null values that make more sense than like nan')
    print('so here it works but not further below, maybe do it in the building as maybe saving is issue')
    sim_mesDist=fetch_main_id(sim_mesDist)
    sim_mesDist['dist_err']=np.maximum(sim_mesDist['plus_err'],-sim_mesDist['minus_err'])
    sim_mesDist.remove_rows(sim_mesDist['dist_err'].mask.nonzero()[0])
    #group by oid
    grouped_mesDist=sim_mesDist.group_by('main_id')
    best_mesDist=sim_mesDist['main_id','dist_value','plus_err','dist_qual','dist_ref'][:0]
    best_mesDist.rename_column('plus_err','dist_err')
    for i in range(len(grouped_mesDist.groups.keys)):
        #sort by quality
        row=grouped_mesDist.groups[i][np.where(grouped_mesDist['dist_prec'].groups[i]==np.max(grouped_mesDist['dist_prec'].groups[i]))][0]
        #take first and add to best_paras
        #which error to take when there are multiples...
        best_mesDist.add_row([row['main_id'],row['dist_value'], row['dist_err'],row['dist_qual'], row['dist_ref']])
    #join with other multimes thingis
    best_paras=best_mesDist#vstack other multi meas tables
    sim_mesDist=sim_mesDist['main_id','dist_value','dist_err','dist_qual','dist_ref']
    
    #--------------------sim_stars---------------------
    #add best para from multiple measurements tables
    stars=ap.table.join(stars,best_paras,keys='main_id',join_type='left')
    
    #--------------sim_h_link ---------------
    sim_h_link=simbad['main_id','parent_oid','h_link_ref','membership']
    #sim_h_link=simbad['main_id','parent_oid','h_link_ref','membership'][np.where(simbad['membership']>50)]
    #downside of excluding lower membership probability objects is that you loose objects with ~ values e.g. alf cen system
    print('all membership values included, use commented out code to change')
    sim_h_link=fetch_main_id(sim_h_link,'parent_oid','parent_main_id')
    sim_h_link.remove_column('parent_oid')
    #null values
    sim_h_link['membership'].fill_value=-1
    sim_h_link['membership']=sim_h_link['membership'].filled()
    #when trying to fetch parent main id I am loosing many stars in simbad
    #ah, so some have missing values... that is the issue. so maybe do this only in h_link...
    #-----------------sim_planets
    temp_sim_planets=simbad['main_id','ids','type'][np.where(simbad['type']=='pl')]
    sim_planets=ap.table.Table(ap.table.unique(temp_sim_planets,keys='main_id'),copy=True)
    #-----------------sim_objects
    sim_objects=ap.table.vstack([sim_planets['main_id','ids','type'],
                             stars['main_id','ids','type']])
    sim_objects['ids']=sim_objects['ids'].astype(object)
    
    
    #--------------sim_sources ---------------
    sim_sources=ap.table.Table()
    tables=[stars, sim_h_link, sim_mesDist,sim_ident]
    for cat,ref in zip(tables,ref_columns):
        sim_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],sim_sources)
    
    #------------------------sim_star_basic------------------------------
    sim_star_basic=stars['main_id','coo_ra','coo_dec','coo_err_angle','coo_err_maj',
                         'coo_err_min','coo_qual','coo_ref',
                         'plx_value','plx_err','plx_qual','plx_ref',
                        'dist_value','dist_err','dist_qual','dist_ref']
    #sim_planet_basic=ap.table.Table()
    #sim_planet_basic['main_id']=sim_planets['main_id']
    #form data into wanted database tables
    save([sim_sources,sim_objects,sim_ident,sim_h_link,sim_star_basic,sim_mesDist],
         ['sim_sources','sim_objects','sim_ident','sim_h_link','sim_star_basic',
          'sim_mesDist'])
    return sim_sources,sim_objects,sim_ident,sim_h_link,sim_star_basic,sim_mesDist

def provider_gk():
    """
    This function obtains the disk data and arranges it in a way 
    easy to ingest into the database.
    """
    provider_name='priv. comm.'
    TAP_service='None'
    provider_bibcode='None'

    gk_disks=ap.io.votable.parse_single_table(
        "data/Grant_absil_2013_.xml").to_table()
    #transforming from string type into object to have variable length
    gk_disks=stringtoobject(gk_disks,212)

    #sorting out objects not within 20pc.
    #removing objects with plx_value=='None'
    gk_disks=gk_disks[np.where(gk_disks['plx_value']!='None')]
    #converting masked plx_value into -99
    gk_disks['plx_value'].fill_value=-99
    gk_disks['plx_value']=gk_disks['plx_value'].astype(float)
    #sorting out everything with plx_value not >50 (corresponding to 50 mas=20pc)
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

    #issue is that main id from gk is not in all cases same anylonger as in simbad
    gk_disks.rename_column('main_id','gk_host_main_id')
    gk_disks=fetch_main_id(gk_disks,colname='gk_host_main_id',name='main_id',oid=False)

    #--------------gk_h_link ---------------
    gk_h_link=gk_disks['id','main_id','disks_ref']
    gk_h_link.rename_columns(['main_id','disks_ref'],['parent_main_id','h_link_ref'])
    gk_h_link.rename_column('id','main_id')
    #--------------gk_objects ---------------
    gk_disks['ids']=gk_disks['id']#because only single id per source given
    gk_objects=gk_disks['id','ids','type']
    gk_objects.rename_column('id','main_id')
    #add stars?
    #--------------gk_ident ---------------
    gk_ident=gk_disks['ids','id','disks_ref']
    #actually would want to use id instad of ids but gets error and ids is same column as id
    gk_ident.rename_columns(['ids','disks_ref'],['main_id','id_ref'])
    #--------------gk_sources ---------------
    gk_sources=sources_table(gk_disks,['disks_ref'],[provider_name,TAP_service,
                                           provider_bibcode])
    #--------------gk_disk_basic----------------
    gk_disk_basic=gk_disks['id','rdisk_bb','e_rdisk_bb','disks_ref']
    #converting from string to float
    for column in ['rdisk_bb','e_rdisk_bb']:
        #replacing 'None' with 'nan' as the first one is not float convertible
        temp_length=len(gk_disk_basic[column][np.where(gk_disk_basic[column]=='None')])
        gk_disk_basic[column][np.where(gk_disk_basic[column]=='None')]=['nan' for i in range(temp_length)]
        gk_disk_basic[column].fill_value='nan'
        gk_disk_basic[column].filled()
        gk_disk_basic[column]=gk_disk_basic[column].astype(float)
    gk_disk_basic.rename_columns(['id','rdisk_bb','e_rdisk_bb','disks_ref'],
                                 ['main_id','rad_value','rad_err','rad_ref'])
    
    save([gk_sources,gk_objects, gk_ident, gk_h_link,gk_disk_basic],
         ['gk_sources','gk_objects', 'gk_ident', 'gk_h_link','gk_disk_basic'])
    
    return gk_sources,gk_objects, gk_ident, gk_h_link,gk_disk_basic

def provider_exo(temp=True):
    """
    This function obtains the exomercat data and arranges it in a way 
    easy to ingest into the database.
    """
    #---------------define provider------------------------
    TAP_service="http://archives.ia2.inaf.it/vo/tap/projects"
    provider_name='Exo-MerCat'
    provider_bibcode='2020A&C....3100370A'
    #---------------define queries---------------------------
    adql_query="""SELECT *
                  FROM exomercat.exomercat"""
    #getting data
    if temp:
        #exomercat=ap.io.votable.parse_single_table("data/raw_exomercat.xml").to_table()
        exomercat=ap.io.ascii.read("data/exomercat_Sep2.csv")
        exomercat=stringtoobject(exomercat,3000)

    else:
        exomercat=query(TAP_service,adql_query)

    #getting main id for planet and host for later comparison of objects from different sources    
    exomercat['planet_main_id']=ap.table.Column(dtype=object, length=len(exomercat))#initializing column
    #issue that string gets truncated in planet main id
    #need to assign more string digit space than 31 -> use object instead of string
    exomercat['host_main_id']=exomercat['main_id']
    for i in range(len(exomercat)):
        if type(exomercat['main_id'][i])!=np.ma.core.MaskedConstant:
            hostname=exomercat['main_id'][i]
        else:
            hostname=exomercat['host'][i]
        if type(exomercat['binary'][i])!=np.ma.core.MaskedConstant:
            exomercat['host_main_id'][i]=hostname+' '+exomercat['binary'][i]
        else:
            exomercat['host_main_id'][i]=hostname
        exomercat['planet_main_id'][i]=exomercat['host_main_id'][i]+' '+exomercat['letter'][i]

    #sorting out objects not within 20pc.
    def sort_out_20pc(cat,colnames):#this does not make sense...
        [sim_objects]=load(['sim_objects'])
        sim_objects.rename_column('main_id','temp')
        #print(cat)
        cat_old=ap.table.Table()
        for colname in colnames:
            cat=ap.table.join(cat,sim_objects['temp','ids'],keys_left=colname,keys_right='temp')
            cat.remove_columns(['temp','ids'])
            cat_old=ap.table.vstack([cat,cat_old])
        cat=ap.table.unique(cat_old,silent=True)
        return cat
        
    exo=exomercat
    exomercat=sort_out_20pc(exomercat,['host_main_id'])

    #removing whitespace in front of main_id and name. doing it down here as here I don't have missing values 
    #in masked column any longer because of sort20pc
    for i in range(len(exomercat)):
        exomercat['planet_main_id'][i]=exomercat['planet_main_id'][i].strip()
        exomercat['main_id'][i]=exomercat['main_id'][i].strip()
        exomercat['name'][i]=exomercat['name'][i].strip()
    #join exomercat on host_main_id and sim_objects main_id
    #using simbad instead of other cataloges to determine if a star is within 20 pc will mean I loose some of them.
    #that is, however, preferrable to having to do the work of checking the literature.
    #a compromise is to keep the list of objects I lost.
    
    #show which elements from exomercat were not found in sim_objects
    exo['name']=exo['name'].astype(object)
    removed_objects=ap.table.setdiff(exo,exomercat,keys=['name'])
    print('length exomercat before join',len(exo['host_main_id']))
    print('length exomercat after join on host_main_id',len(exomercat['host_main_id']))
    print('number of removed_objects',len(removed_objects['host_main_id']))###----------------hm 7000 there, so join didn't work----------
    save([removed_objects],['exomercat_removed_objects'])
    print('tbd: improve to not loose that many objects') 
    
    #-------------exo_ident---------------
    #['main_id','id','id_ref']
    exo_ident=exomercat['planet_main_id','name']
    exo_ident.rename_columns(['planet_main_id','name'],['main_id','id'])
    for i in range(len(exomercat)):
        if exomercat['planet_main_id'][i]!=exomercat['name'][i]:
            exo_ident.add_row([exomercat['planet_main_id'][i],exomercat['planet_main_id'][i]])
    exo_ident['id_ref']=[provider_bibcode for j in range(len(exo_ident))]
    print(exo_ident[np.where(exo_ident['main_id']=='Wolf  940 b')])
    print('I have a wrong double object in exo_ident because there are different amount of white spaces between catalog and number')
    #-------------exo_objects---------------
    print('tbd at one point: I think I want to add hosts to object')
    exo_objects=ap.table.Table(names=['main_id','ids'],dtype=[object,object])
    grouped_exo_ident=exo_ident.group_by('main_id')
    ind=grouped_exo_ident.groups.indices
    for i in range(len(ind)-1):#-1 is needed because else ind[i+1] is out of bonds
        ids=[]
        for j in range(ind[i],ind[i+1]):
            ids.append(grouped_exo_ident['id'][j])
        ids="|".join(ids)
        exo_objects.add_row([grouped_exo_ident['main_id'][ind[i]],ids])
    exo_objects['type']=['pl' for j in range(len(exo_objects))]

    #-------------------exo_mesMass---------------------
    exomercat['mass_max'].fill_value=999
    exomercat['mass_min'].fill_value=-999
    exomercat['mass_max']=exomercat['mass_max'].filled()
    exomercat['mass_min']=exomercat['mass_min'].filled()
    exomercat['mass_err']=np.maximum(exomercat['mass_max'],-exomercat['mass_min'])
    exo_mesMass=exomercat['planet_main_id','mass','mass_err','mass_url']
    exo_mesMass.rename_columns(['planet_main_id','mass','mass_url'],
                                    ['main_id','mass_val','mass_ref'])
    print('issue is that I have many null value ones (null expr is 1w+20 for val and 999 for err) in here but should only contain those with actual measurements')
    #exo_mesMass.remove_rows(exo_mesMass[np.where(mesMass['mass_val']==1e+20)])
    
    grouped_mesMass=exo_mesMass.group_by('main_id')
    #____ ahm here, dont forget to delete some stuff below
    best_mesMass=exo_mesMass['main_id','mass_val','mass_err','mass_ref'][:0]
    for i in range(len(grouped_mesMass.groups.keys)):
        #sort by quality
        row=grouped_mesMass.groups[i][np.where(grouped_mesMass['mass_err'].groups[i]==np.min(grouped_mesMass['mass_err'].groups[i]))][0]
        #take first and add to best_paras
        #which error to take when there are multiples...
        best_mesMass.add_row([row['main_id'],row['mass_val'], row['mass_err'],row['mass_ref']])
    #join with other multimes thingis
    best_paras=best_mesMass#vstack other multi meas tables
    
    #some of the stars are not in simbad
    #-------------exo_h_link---------------
    #['child_object_idref','parent_object_idref',
               #'h_link_source_idref','h_link_ref','membership']
    exo_h_link=exomercat['planet_main_id', 'host_main_id']
    exo_h_link.rename_columns(['planet_main_id','host_main_id'],
                              ['main_id','parent_main_id'])
    exo_h_link['h_link_ref']=[provider_bibcode for j in range(len(exo_h_link))]

    #I decided not to take catalog as h_link_ref, as it would introduce 16 different sources when it
    #just means providers of provider.
    #the correct way would be to add the four different databases as sources and have membership from
    #status string column
    #this, however, is more detail than needed for just prototype database
    #-------------exo_planet_basic
    #exchanged bestmass with mass because had issues of different references for same value and object
    #did not solve issue, need to create like for simbad special table for mass values and then only
    #add single ones to planet_basic
    exo_planet_basic=best_mesMass
    
    #-------------exo_sources---------------
    ref_columns=[['mass_url'],['h_link_ref'],['id_ref']]
    exo_sources=ap.table.Table()
    tables=[exomercat, exo_h_link,exo_ident]
    for cat,ref in zip(tables,ref_columns):
        exo_sources=sources_table(cat,ref,[provider_name,TAP_service,
                                           provider_bibcode],exo_sources) 
        
    #print('to be saved:',[exo_sources,exo_objects,exo_ident,exo_h_link,exo_planet_basic])

    save([exo_sources,exo_objects,exo_ident,exo_h_link,exo_planet_basic,exo_mesMass],
         ['exo_sources','exo_objects','exo_ident','exo_h_link','exo_planet_basic','exo_mesMass'])
        
    return exo_sources,exo_objects,exo_ident,exo_h_link,exo_planet_basic,exo_mesMass

#------------------------provider combining-----------------
def building(sim,gk,exo,temp=False):
    """
    This function builds from the input parameters the tables
    for the LIFE database.
    """
    #creates empty tables as needed for final database ingestion
    init=initialize_database_tables()
    
    #initializes 8 table objects
    #corresponding to 'sources','objects','ident','h_link','star_basic','planet_basic','disk_basic','mesDist'
    cat=[ap.table.Table() for i in range(9)]
    
    #for the sources and objects stacks tables from different providers keeping only unique values
    #then create identifiers for those tables
    
    cat[0]=ap.table.vstack([init[0],sim[0]])
    cat[0]=ap.table.vstack([cat[0],gk[0]])
    cat[0]=ap.table.vstack([cat[0],exo[0]])
    if len(cat[0])>0:
        cat[0]=ap.table.unique(cat[0],silent=True)
        cat[0]['source_id']=[j+1 for j in range(len(cat[0]))]
    
    def idsjoin(cat,column_ids1,column_ids2):
        """
        This function merges the identifiers from two different columns into one.
        """
        cat['ids']=ap.table.Column(dtype=object, length=len(cat))#initializing column
        for column in [column_ids1,column_ids2]:
            if type(cat[column])==ap.table.column.MaskedColumn:
                cat[column].fill_value=''
                cat[column]=cat[column].filled()
        for i in range(len(cat)):
            ids1=cat[column_ids1][i].split('|')#splitting object into list of elements
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
    #cat[1]=ap.table.vstack([init[1],sim[1]]) 
    cat[1]=sim[1]#removed vstack with init to not have object_id as is empty anyways
    cat[1]=ap.table.join(cat[1],gk[1],keys='main_id',join_type='outer')

    print('getting warning about column type mergeing of string types')
    print('ok could not solve the warning message about type merging in type part but stuff seems generally to work')

    def objectmerging(cat):
        cat=idsjoin(cat,'ids_1','ids_2')
        #merging types
        cat['type']=ap.table.Column(dtype=object, length=len(cat))#initializing column
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
    
    print('At one point I would like to be able to merge objects with main_id NAME Proxima Centauri b and Proxima Centauri b')
    print('should work with simbad upload')        
    #now get source_idref here
    def match(cat,sources,paras,provider):
        """
        This function joins the source identifiers to the in paras specified parameters of cat.
        """
        #for all parameters specified
        for para in paras:
            #if they have reference columns
            if para+'_ref' in cat.colnames:
                #if those reference columns are masked
                if type(cat[para+'_ref'])==ap.table.column.MaskedColumn:
                    cat[para+'_ref'].fill_value=''
                    cat=cat.filled()
                #join to each reference parameter its source_id
                cat=ap.table.join(cat,
                    sources['ref','source_id'][np.where(sources['provider_name']==provider)],
                    keys_left=para+'_ref',keys_right='ref',join_type='left')
                #renaming column to specify to which parameter the source_id correspond
                cat.rename_column('source_id',f'{para}_source_idref')
                #deleting double column containing reference information
                cat.remove_columns('ref')
                #in case the para_value entry is masked this if environment will put the source_id entry to null
                if para+'_value' in cat.colnames:
                    if type(cat[para+'_value'])==ap.table.column.MaskedColumn:
                        for i in cat[para+'_value'].mask.nonzero()[0]:
                            cat[f'{para}_source_idref'][i]=0
        return cat
    paras=[['id'],['h_link'],['coo','plx','dist'],['mass'],['rad'],['dist'],['mass']]
    
    for i in range(2,9):
        #replacing ref with source_idref columns
        #q for markus: do I need to do this or can dachs do that for me?
        #I mean creating an identifier for the refs
        
        #getting source_idref to each ref
        sim[i]=match(sim[i],cat[0],paras[i-2],'SIMBAD')
        gk[i]=match(gk[i],cat[0],paras[i-2],'priv. comm.')
        exo[i]=match(exo[i],cat[0],paras[i-2],'Exo-MerCat') 
        
        #joining data from different providers
        cat[i]=ap.table.vstack([init[i],sim[i]])#I do this to get those columns that are empty in the data
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
                #first remove the child_object_idref we got from empty initialization
                #yes there should be a more elegant way to do this but at least it works
                cat[i].remove_column('child_object_idref')
                cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],keys='main_id',
                   join_type='left')
                cat[i].rename_columns(['object_id','main_id'],
                                      ['child_object_idref','child_main_id'])
                
                #expanding from parent_main_id to parent_object_idref
                cat[i].remove_column('parent_object_idref')
                #kick out any h_link rows where parent_main_id not in objects (e.g. clusters)
                cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],
                   keys_left='parent_main_id',keys_right='main_id')
                cat[i].remove_column('main_id')#removing because same as parent_main_id
                cat[i].rename_column('object_id','parent_object_idref')
                #null values
                cat[i]['membership'].fill_value=-1
                cat[i]['membership']=cat[i]['membership'].filled()
            
            #for all the other tables add object_idref
            else:
                #first remove the object_idref we got from empty initialization
                #yes there should be a more elegant way to do this but at least it works
                cat[i].remove_column('object_idref')
                cat[i]=ap.table.join(cat[i],cat[1]['object_id','main_id'],join_type='left')
                cat[i].rename_column('object_id','object_idref')
            if i==4:#--------------------star_basic--------------------------
                #choosing all objects with type star or system
                #I am just adding main_id because I have not found out how to do join with just one column table
                stars=cat[1]['object_id','main_id'][np.where(cat[1]['type']=='st')]
                systems=cat[1]['object_id','main_id'][np.where(cat[1]['type']=='sys')]
                temp=ap.table.vstack([stars,systems])
                temp.rename_columns(['object_id','main_id'],['object_idref','temp'])
                cat[i]=ap.table.join(cat[i],temp,join_type='outer',keys='object_idref')
                cat[i].remove_column('temp')
            if i==5:#--------------------planet_basic--------------------------
                temp=cat[1]['object_id','main_id'][np.where(cat[1]['type']=='pl')]
                temp.rename_columns(['object_id','main_id'],['object_idref','temp'])
                cat[i]=ap.table.join(cat[i],temp,keys='object_idref',join_type='outer')
                cat[i].remove_column('temp')
        else:
            print('error: empty table')
    if temp:
        cat[5]=ap.table.unique(cat[5],keys='object_idref')
    #add missing data information to basic
    #pl=objects[np.where(objects['type']=='pl')]
    #di=objects[np.where(objects['type']=='di')]
    #st=objects-pl-di
    #if len(star_basic)< len(st):
        #for obj in st['main_id']:
            #if obj not in star_basic['main_id']:
                #star_basic.add_row()
    save(cat,
         ['sources','objects','ident','h_link','star_basic','planet_basic','disk_basic','mesDist','mesMass'])
    return cat

##############################################################################
#################initialize###################################################
sources,objects,ident,h_link,star_basic,planet_basic,disk_basic,mesDist,mesMass=initialize_database_tables()
save([sources,objects,ident,h_link,star_basic,planet_basic,disk_basic,mesDist,mesMass],
         ['sources','objects','ident','h_link','star_basic',
          'planet_basic','disk_basic','mesDist','mesMass'])
#################get source data###################################################
sim_sources,sim_objects,sim_ident,sim_h_link,sim_star_basic,sim_mesDist=provider_simbad()
gk_sources,gk_objects, gk_ident, gk_h_link,gk_disk_basic=provider_gk()
exo_sources,exo_objects,exo_ident,exo_h_link,exo_planet_basic,exo_mesMass=provider_exo()
#################building db###################################################
empty=ap.table.Table()
sim=[sim_sources,sim_objects,sim_ident,sim_h_link,sim_star_basic,empty[:],empty[:],sim_mesDist,empty[:]]
gk=[gk_sources,gk_objects, gk_ident, gk_h_link,empty[:],empty[:],gk_disk_basic,empty[:],empty[:]]
exo=[exo_sources,exo_objects, exo_ident, exo_h_link,empty[:],exo_planet_basic,empty[:],empty[:],exo_mesMass]

sources,objects,ident,h_link,star_basic,planet_basic,disk_basic,mesDist,mesMass=building(sim,gk,exo)

