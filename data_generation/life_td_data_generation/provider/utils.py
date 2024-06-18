""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from pyvo.dal import TAPService
from astropy.table import Table, column, unique, vstack, join
from datetime import datetime

#self created modules
from utils.utils import load

additional_data_path='../../additional_data/'

#------------------------------provider helper functions----------------
def query(link,adql_query,catalogs=[],no_description=True):
    """
    Performs a query via TAP on the service given in the link parameter.
    
    If a list of tables is given in the catalogs parameter,
    those are uploaded to the service beforehand.
    
    :param str link: Service access URL.
    :param str adql_query: Query to be asked of the external database service
         in ADQL.
    :param catalogs: List of astropy tables to be uploaded to the 
        service.
    :type catalogs: list(astropy.table.table.Table)
    :param bool no_description: Defaults to True, wether description gets removed
    :returns: Result of the query.
    :rtype: astropy.table.table.Table
    """
    
    #defining the vo service using the given link
    service = TAPService(link)
    #without upload tables
    if catalogs==[]:
        result=service.run_async(adql_query.format(**locals()), maxrec=1600000)
    #with upload tables
    else:
        tables={}
        for i in range(len(catalogs)):
            tables.update({f"t{i+1}":catalogs[i]})
        result = service.run_async(adql_query,uploads=tables,timeout=None,
                                   maxrec=1600000)
    cat=result.to_table()
    
    # removing descriptions because merging of data leaves wrong description 
    for col in cat.colnames:
        if no_description:
            cat[col].description=''
    
    return cat

def sources_table(cat,ref_columns,provider,old_sources=Table()):
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
        cat_sources=Table() 
        #initialization of list to store reference information
        cat_reflist=[] 
        #for all the columns given add reference information 
        for k in range(len(ref_columns)):
            #In case the column has elements that are masked skip those
            if type(cat[ref_columns[k]])==column.MaskedColumn:
                cat_reflist.extend(
                    cat[ref_columns[k]][np.where(
                            cat[ref_columns[k]].mask==False)])
            else:
                cat_reflist.extend(cat[ref_columns[k]])
        # add list of collected references to the table and call the 
        # column ref
        cat_sources['ref']=cat_reflist
        cat_sources=unique(cat_sources)
        #attaches service information
        cat_sources['provider_name']=[provider for j in range(
                len(cat_sources))]
        #combine old and new sources into one table
        sources=vstack([old_sources,cat_sources])
        sources=unique(sources) #remove double entries
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
        [sim]=load(['sim_objects'])
        sim.rename_columns(['main_id','ids'],['temp1','temp2'])
        cat=join(cat,sim['temp1','temp2'],
                      keys_left=colname,keys_right='temp1')
        cat.remove_columns(['temp1','temp2'])
    else:
        [sim]=load(['sim_ident'])
        sim.rename_columns(['id'],['temp1'])
        cat=join(cat,sim['temp1','main_id'],
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
    
    if type(cat[colname])==column.MaskedColumn:
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

