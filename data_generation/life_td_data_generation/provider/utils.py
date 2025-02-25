""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np  #arrays
from pyvo.dal import TAPService
from astropy.table import Table, column, unique, vstack, join, table, MaskedColumn, Column
from datetime import datetime
from typing import List

#self created modules
from utils.io import load


def initiate_columns(table, columns, types, maskc):
    for i in range(len(columns)):
        if maskc[i]:
            table[columns[i]] = MaskedColumn(dtype=types[i], length=len(table))
        else:
            table[columns[i]] = Column(dtype=types[i], length=len(table))
    return table


def create_provider_table(provider_name, provider_url, provider_bibcode,
                          provider_access=datetime.now().strftime('%Y-%m-%d')):
    print(f'Trying to create {provider_name} tables from {provider_access}...')
    provider_table = Table()
    provider_table['provider_name'] = [provider_name]
    provider_table['provider_url'] = [provider_url]
    provider_table['provider_bibcode'] = [provider_bibcode]
    provider_table['provider_access'] = [provider_access]
    return provider_table


def query(link: str, adql_query: str, upload_tables: List[table.Table] = []) -> table.Table:
    """
    Performs a query via TAP on the service given in the link parameter.
    
    If a list of tables is given in the catalogs parameter,
    those are uploaded to the service beforehand.
    
    :param str link: Service access URL.
    :param str adql_query: Query to be asked of the external database service
         in ADQL.
    :param upload_tables: List of astropy tables to be uploaded to the 
        service.
    :type upload_tables: list(astropy.table.table.Table)
    :returns: Result of the query.
    :rtype: astropy.table.table.Table
    """

    service = TAPService(link)
    if upload_tables == []:
        result = service.run_async(adql_query.format(**locals()), maxrec=1600000)
    else:
        tables = {}
        for i in range(len(upload_tables)):
            tables.update({f"t{i + 1}": upload_tables[i]})
        result = service.run_async(adql_query, uploads=tables, timeout=None,
                                   maxrec=1600000)
    return result.to_table()


def remove_catalog_description(cat: table.Table, no_description) -> table.Table:
    """
    Removes description meta data of catalog columns.
    
    This is useful if by merging different catalogs the description of the 
    resulting catalog is not correct any longer. 
    
    :param cat: Catalog
    :type cat: astropy.table.table.Table
    :returns: Catalog with no description of columns.
    :rtype: astropy.table.table.Table
    """

    for col in cat.colnames:
        if no_description:
            cat[col].description = ''
    return cat


def fill_sources_table(cat: table.Table, ref_columns: List[str], provider: str,
                       old_sources: table.Table = Table()) -> table.Table:
    """
    Creates or updates the source table out of the given references.
    
    The entries are unique and the columns consist out of the
    reference and provider_name.
    
    :param cat: Table on which the references should be gathered.
    :type cat: astropy.table.table.Table
    :param ref_columns: Header of the columns containing reference 
        information.
    :type ref_columns: list(str)
    :param str provider: Provider name.
    :param old_sources: Previously created reference table.
    :type old_sources: astropy.table.table.Table
    :return: Table containing references and provider information.
    :rtype: astropy.table.table.Table
    """

    if len(cat) > 0:
        # table initialization to prevent error messages when assigning 
        # columns
        cat_sources = Table()
        #initialization of list to store reference information
        cat_reflist = []
        #for all the columns given add reference information 
        for k in range(len(ref_columns)):
            #In case the column has elements that are masked skip those
            if type(cat[ref_columns[k]]) == column.MaskedColumn:
                cat_reflist.extend(
                    cat[ref_columns[k]][np.where(
                        cat[ref_columns[k]].mask == False)])
            else:
                cat_reflist.extend(cat[ref_columns[k]])
        # add list of collected references to the table and call the 
        # column ref
        cat_sources['ref'] = cat_reflist
        cat_sources = unique(cat_sources)
        #attaches service information
        cat_sources['provider_name'] = [provider for j in range(
            len(cat_sources))]
        #combine old and new sources into one table
        sources = vstack([old_sources, cat_sources])
        sources = unique(sources)  #remove double entries
    else:
        sources = old_sources
    return sources


def create_sources_table(tables, ref_columns, provider_name):
    #--------------creating output table sim_sources -------------------
    sources = Table()
    for cat, ref in zip(tables, ref_columns):
        sources = fill_sources_table(cat, ref, provider_name, sources)
    return sources


class OidCreator:
    """
    Create adql query for fetch_main_id function using oid column.
    """

    def __init__(self, name, colname):
        self.name = name
        self.colname = colname

    def create_main_id_query(self):
        return 'SELECT b.main_id AS ' + self.name + """,t1.*
                FROM basic AS b
                JOIN TAP_UPLOAD.t1 ON b.oid=t1.""" + self.colname


class IdentifierCreator:
    """
    Create adql query for fetch_main_id function using identifier column.
    """

    def __init__(self, name, colname):
        self.name = name
        self.colname = colname

    def create_main_id_query(self):
        return 'SELECT b.main_id AS ' + self.name + """,t1.*
                    FROM basic AS b
                    JOIN ident ON ident.oidref=b.oid
                        JOIN TAP_UPLOAD.t1 ON ident.id=t1.""" + self.colname


#looks better but don't think this will run. issue is that I pass variables to a class that doesn't take any

def fetch_main_id(cat: table.Table, id_creator=OidCreator(name='main_id', colname='oid')) -> table.Table:
    """
    Joins main_id from simbad to the column colname. 
    
    Returns the whole table cat but without any rows where no simbad 
    main_id was found.
    
    :param cat: Astropy table containing column colname.
    :type cat: astropy.table.table.Table
    :param id_creator: OidCreator or IdentifierCreator object.
    :return: Table with all main SIMBAD identifiers that could be found 
        in column "name".
    :rtype: astropy.table.table.Table
    """

    #improvement idea to be performed at one point
    # tbd option to match on position instead of main_id or oid
    #SIMBAD TAP service
    TAP_service = "http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    #performing query using external function
    main_id_query = id_creator.create_main_id_query()
    cat = query(TAP_service, main_id_query, [cat])
    return cat


def distance_cut(cat: table.Table, colname: str, main_id: bool = True):
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
        [sim] = load(['sim_objects'])
        sim.rename_columns(['main_id', 'ids'], ['temp1', 'temp2'])
        cat = join(cat, sim['temp1', 'temp2'],
                   keys_left=colname, keys_right='temp1')
        cat.remove_columns(['temp1', 'temp2'])
    else:
        [sim] = load(['sim_ident'])
        sim.rename_columns(['id'], ['temp1'])
        cat = join(cat, sim['temp1', 'main_id'],
                   keys_left=colname, keys_right='temp1')
        cat.remove_columns(['temp1'])
    return cat


def nullvalues(cat, colname, nullvalue, verbose=False):
    """
    This function fills masked entries of specified column. 
    
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

    if type(cat[colname]) == column.MaskedColumn:
        cat[colname].fill_value = nullvalue
        cat[colname] = cat[colname].filled()
    elif verbose:
        print(colname, 'is no masked column')
    return cat


def replace_value(cat, column, value, replace_by):
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

    cat[column][np.where(cat[column] == value)] = \
        [replace_by for i in range(
            len(cat[column][np.where(cat[column] == value)]))]
    return cat


def ids_from_ident(ident, objects):
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

    grouped_ident = ident.group_by('main_id')
    ind = grouped_ident.groups.indices
    for i in range(len(ind) - 1):
        # -1 is needed because else ind[i+1] is out of bonds
        ids = []
        for j in range(ind[i], ind[i + 1]):
            ids.append(grouped_ident['id'][j])
        ids = "|".join(ids)
        objects.add_row([grouped_ident['main_id'][ind[i]], ids])
    return objects


def lower_quality(qual):
    if qual == 'A':
        qual = 'B'
    elif qual == 'B':
        qual = 'C'
    elif qual == 'C':
        qual = 'D'
    elif qual == 'D':
        qual = 'E'
    return qual

def teff_st_spec_assign_quality(gaia_mes_teff_st_spec):
    interval = 41 * 9 / 5.
    gaia_mes_teff_st_spec['teff_st_qual'] = ['?' for j in range(len(gaia_mes_teff_st_spec))]
    for i, flag in enumerate(gaia_mes_teff_st_spec['flags_gspspec']):
        summed = 0
        for j in flag:
            summed += int(j)
        if summed in range(0, int(interval) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'A'
        elif summed in range(int(interval) + 1, int(interval * 2) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'B'
        elif summed in range(int(interval * 2) + 1, int(interval * 3) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'C'
        elif summed in range(int(interval * 3) + 1, int(interval * 4) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'D'
        elif summed in range(int(interval * 4) + 1, int(interval * 5) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'E'
    return gaia_mes_teff_st_spec

def assign_quality_elementwise(exo_helptab, para, i):
    qual = 'B'
    if exo_helptab[para + '_max'][i] == 1e+20:
        qual = lower_quality(qual)
    if exo_helptab[para + '_min'][i] == 1e+20:
        qual = lower_quality(qual)
    return qual

def exo_assign_quality(exo_helptab):
    for para in ['mass', 'msini']:
        exo_helptab[para + '_pl_qual'] = MaskedColumn(dtype=object, length=len(exo_helptab))
        exo_helptab[para + '_pl_qual'] = ['?' for j in range(len(exo_helptab))]
        for i in range(len(exo_helptab)):
            exo_helptab[para + '_pl_qual'][i] = assign_quality_elementwise(exo_helptab, para, i)
    return exo_helptab

def assign_quality(table, column = '',special_mode = ''):
    if special_mode == 'teff_st_spec':
        # assign quality based on reliability flags
        table = teff_st_spec_assign_quality(table)
    elif special_mode == 'exo':
        table = exo_assign_quality(table)
    elif special_mode == 'gaia_binary':
        table[column]=['B' if table['binary_flag'][j] == 'True' \
             else 'E' for j in range(len(table))]
    elif special_mode == 'wds_sep1':
        table[column] = ['C' if type(j) != np.ma.core.MaskedConstant
                         else 'E' for j in wds_mes_sep_ang1['sep_ang_obs_date']]
    elif special_mode == 'wds_sep2':
        # add a quality to sep1 which is better than sep2. because newer measurements should be better.
        table[column] = ['B' if type(j) != np.ma.core.MaskedConstant
                         else 'E' for j in wds_mes_sep_ang2['sep_ang_obs_date']]
    else:
        if column == 'coo_gal_qual':
            # not sure how to transform quality from ra and dec to gal
            qual = '?'
        elif special_mode in ['teff_st_phot','radius_st_flame','mass_st_flame']:
            # just assumption. to do -> think about it more
            qual = 'B'
        elif special_mode in ['model','wds_binary']:
            # quality assumption of model
            qual = 'C'
        elif special_mode == 'sim_binary':
            qual = 'D'
        table[column] == [qual for j in range(len(table))]

    return table

