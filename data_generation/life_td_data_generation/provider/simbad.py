""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np  #arrays
from astropy.table import setdiff, Table, join, vstack, unique, MaskedColumn
from datetime import datetime

#self created modules
from utils.io import save
from provider.utils import fetch_main_id, OidCreator, fill_sources_table, create_sources_table, query, nullvalues, \
    replace_value, create_provider_table
from provider.assign_quality_funcs import assign_quality
from sdata import empty_dict

#---------------define queries--------------------------------------
adql_queries_moduls = {
    'select_statement': """SELECT b.main_id,b.ra AS coo_ra,b.dec AS coo_dec,
        b.coo_err_angle, b.coo_err_maj, b.coo_err_min,b.oid,
        b.coo_bibcode AS coo_ref, b.coo_qual,b.sp_type AS sptype_string,
        b.sp_qual AS sptype_qual, b.sp_bibcode AS sptype_ref,
        b.plx_err, b.plx_value, b.plx_bibcode AS plx_ref,b.plx_qual,
        h_link.membership, h_link.parent AS parent_oid,
        h_link.link_bibcode AS h_link_ref, a.otypes,ids.ids,
        f.I as mag_i_value, f.J as mag_j_value, f.K as mag_k_value
        """,
    'tables_statement': """
    FROM basic AS b
        JOIN ids ON b.oid=ids.oidref
            JOIN alltypes AS a ON b.oid=a.oidref
                LEFT JOIN h_link ON b.oid=h_link.child
                    LEFT JOIN allfluxes AS f ON b.oid=f.oidref

    """}


def main_adql_queries(plx_cut):
    return adql_queries_moduls['select_statement'] + \
        adql_queries_moduls['tables_statement'] + \
        'WHERE b.plx_value >=' + str(plx_cut)


adql_upload_queries = {
    'sy_without_plx_but_child_with_upload': adql_queries_moduls['select_statement'] + \
                                            adql_queries_moduls['tables_statement'] + \
                                            """JOIN TAP_UPLOAD.t1 ON b.oid=t1.parent_oid
                                                 WHERE (b.plx_value IS NULL) AND (otype='**..')""",
    'pl_without_plx_but_host_with_upload': adql_queries_moduls['select_statement'] + \
                                           adql_queries_moduls['tables_statement'] + \
                                           """JOIN TAP_UPLOAD.t1 ON b.oid=t1.oid
                                             WHERE (b.plx_value IS NULL) AND (otype='Pl..')""",
    'ids_from_upload': """SELECT id, t1.*
                          FROM ident
                                   JOIN TAP_UPLOAD.t1 ON oidref = t1.oid"""
}

def create_simbad_helpertable(distance_cut_in_pc, test_objects):
    """
    Creates helper table.
    
    :param float distance_cut_in_pc: Distance up to which stars are included.
    :param test_objects: Objects to be tested where they drop out of the criteria or make it till the end.
    :type test_objects: list(str) 
    :returns: Helper table and dictionary of database table names and tables.
    :rtype: astropy.table.table.Table, dict(str,astropy.table.table.Table)
    """
    plx_in_mas_cut = 1000. / distance_cut_in_pc
    #making cut a bit bigger for correct treatment of objects on boundary
    plx_cut = plx_in_mas_cut - plx_in_mas_cut / 10.

    sim = empty_dict.copy()
    sim['provider'] = create_provider_table('SIMBAD',
                                            "http://simbad.u-strasbg.fr:80/simbad/sim-tap",
                                            '2000A&AS..143....9W')


    #------------------querrying----------------------------------------
    #perform query for objects with in distance given
    sim_helptab = query(sim['provider']['provider_url'][0], main_adql_queries(plx_cut))
    #querries parent and children objects with no parallax value
    parents_without_plx = query(sim['provider']['provider_url'][0],
                                adql_upload_queries['sy_without_plx_but_child_with_upload'],
                                [sim_helptab])
    children_without_plx = query(sim['provider']['provider_url'][0],
                                 adql_upload_queries['pl_without_plx_but_host_with_upload'],
                                 [sim_helptab])

    test_objects = np.array(test_objects)
    if len(test_objects) > 0:
        print('in sim through plx query',
              test_objects[np.where(np.isin(test_objects,
                                            sim_helptab['main_id']))])
        print('in sim through child plx query',
              test_objects[np.where(np.isin(test_objects,
                                            parents_without_plx['main_id']))])
        print('in sim through parent plx query',
              test_objects[np.where(np.isin(test_objects,
                                            children_without_plx['main_id']))])

    #adding of no_parallax objects to rest of simbad query objects

    sim_helptab = vstack([sim_helptab, parents_without_plx])
    sim_helptab = vstack([sim_helptab, children_without_plx])

    print(' sorting object types...')

    #sorting from object type into star, system and planet type
    sim_helptab['type'] = ['None' for i in range(len(sim_helptab))]
    sim_helptab['binary_flag'] = ['False' for i in range(len(sim_helptab))]
    to_remove_list = []
    removed_otypes = []
    for i in range(len(sim_helptab)):
        #planets
        if "Pl" in sim_helptab['otypes'][i]:
            sim_helptab['type'][i] = 'pl'
        #stars
        elif "*" in sim_helptab['otypes'][i]:
            #system containing multiple stars
            if "**" in sim_helptab['otypes'][i]:
                sim_helptab['type'][i] = 'sy'
                sim_helptab['binary_flag'][i] = 'True'
            #individual stars
            else:
                sim_helptab['type'][i] = 'st'
        else:
            removed_otypes.append(sim_helptab['otypes'][i])
            #most likely single brown dwarfs
            #storing information for later removal from table called simbad
            to_remove_list.append(i)
    #removing any objects that are neither planet, star nor system in type
    if to_remove_list != []:
        print('removing', len(removed_otypes), ' objects that had object types:',
              list(set(removed_otypes)))
        print('example object of them:', sim_helptab['main_id'][to_remove_list[0]])
        sim_helptab.remove_rows(to_remove_list)

    if len(test_objects) > 0:
        print('in sim through otype criteria',
              test_objects[np.where(np.isin(test_objects,
                                            sim_helptab['main_id']))])

    return sim_helptab, sim


def stars_in_multiple_system(cat, sim_h_link, all_objects):
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
    parents = sim_h_link['parent_main_id', 'main_id', 'h_link_ref'][:]
    parents.rename_column('main_id', 'child_main_id')
    parents.rename_column('parent_main_id', 'main_id')
    sy_wo_child = setdiff(cat['main_id', 'type', 'sptype_string'][:],
                          parents[:], keys=['main_id'])
    #that don t have children: sy_wo_child['main_id','type']
    #list of those with children
    sy_w_child = join(parents[:],
                      cat['main_id', 'type', 'sptype_string'][:],
                      keys=['main_id'])
    #list of those with children joined with type of child
    all_objects.rename_columns(['type', 'main_id'],
                               ['child_type', 'child_main_id'])
    sy_w_child = join(sy_w_child[:],
                      all_objects['child_type', 'child_main_id'][:],
                      keys=['child_main_id'], join_type='left')
    #remove all where type child is not pl
    sy_w_child_pl = sy_w_child[np.where(sy_w_child['child_type'] == 'pl')]
    if len(sy_w_child_pl) == 0:
        #no systems with child of type planet
        sy_wo_child_st = sy_wo_child
    else:
        #join with list of sy that dont habe children
        sy_wo_child_st = vstack([sy_wo_child[:], sy_w_child_pl[:]])
        sy_wo_child_st.remove_column('child_type')
    #systems that don t have children except planets: sy_wo_child_st
    #no + in sptype_string because that is another indication of binarity
    temp = [len(i.split('+')) == 1 for i in sy_wo_child_st['sptype_string']]
    #have it as an array of bools 
    temp = np.array(temp)
    #have it as lisit of indices 
    temp = list(np.where(temp == True)[0])
    single_sptype = sy_wo_child_st[:][temp]
    #and no + in spectral type: single_sptype['main_id','type']      
    cat['type'][np.where(np.isin(cat['main_id'],
                                 single_sptype['main_id']))] = \
        ['st' for j in range(len(cat[np.where(np.isin(cat['main_id'],
                                                      single_sptype['main_id']))]))]
    return cat


def creating_helpertable_stars(sim_helptab, sim):
    """
    Creates another helper table.
    
    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict(str,astropy.table.table.Table)
    :returns: Helper table.
    :rtype: astropy.table.table.Table
    """
    temp_stars = sim_helptab[np.where(sim_helptab['type'] != 'pl')]
    #removing double objects (in there due to multiple parents)
    stars = Table(unique(temp_stars, keys='main_id'), copy=True)
    return stars


def expanding_helpertable_stars(sim_helptab, sim, stars):
    """
    Adds data to helper table stars.
    
    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict(str,astropy.table.table.Table)
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.table.Table
    :returns: Helper table stars.
    :rtype: astropy.table.table.Table
    """
    #--------------------creating helper table sim_stars----------------
    #updating multiplicity object type
    #no children and sptype does not contain + -> type needs to be st

    #all objects in stars table: stars['main_id','type']
    stars[np.where(stars['type'] == 'sy')] = stars_in_multiple_system(
        stars[np.where(stars['type'] == 'sy')], sim['h_link'][:],
        sim_helptab['main_id', 'type'][:])

    # binary_flag 'True' for all stars with parents
    # meaning stars[main_id] in sim_h_link[child_main_id] 
    #-> stars[binary_flag]=='True'    
    stars['binary_flag'][np.where(np.isin(stars['main_id'],
                                          sim['h_link']['main_id']))] = \
        ['True' for j in range(len(stars[np.where(
            np.isin(stars['main_id'], sim['h_link']['main_id']))]))]

    #change null value of plx_qual
    stars['plx_qual'] = stars['plx_qual'].astype(object)
    stars = replace_value(stars, 'plx_qual', '', stars['plx_qual'].fill_value)

    for band in ['i', 'j', 'k']:
        #initiate some of the ref columns
        stars[f'mag_{band}_ref'] = MaskedColumn(dtype=object,
                                                length=len(stars),
                                                mask=[True for j in range(len(stars))])
        #add simbad reference where no other is given
        stars[f'mag_{band}_ref'][np.where(
            stars[f'mag_{band}_value'].mask == False)] = [
            sim['provider']['provider_bibcode'][0] for j in range(len(
                stars[f'mag_{band}_ref'][np.where(
                    stars[f'mag_{band}_value'].mask == False)]))]

    stars = replace_value(stars, 'plx_ref', '', sim['provider']['provider_bibcode'][0])
    stars = replace_value(stars, 'sptype_ref', '',
                          sim['provider']['provider_bibcode'][0])
    stars = replace_value(stars, 'coo_ref', '', sim['provider']['provider_bibcode'][0])

    stars['binary_ref'] = [sim['provider']['provider_bibcode'][0] for j in range(
        len(stars))]
    stars = assign_quality(stars,'binary_qual',special_mode='sim_binary')
    return stars


def create_ident_table(sim_helptab, sim):
    """
    Creates identifier table.
    
    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict(str,astropy.table.table.Table)
    :returns: Identifier table.
    :rtype: astropy.table.table.Table
    """

    sim_ident = query(sim['provider']['provider_url'][0],
                      adql_upload_queries['ids_from_upload'],
                      [sim_helptab['oid', 'main_id'][:].copy()])  #adds column id
    sim_ident['id_ref'] = [sim['provider']['provider_bibcode'][0] \
                           for j in range(len(sim_ident))]
    sim_ident.remove_column('oid')
    return sim_ident


def create_h_link_table(sim_helptab, sim, stars):
    """
    Creates hierarchical link table.
    
    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.table.Table
    :param sim: Dictionary of database table names and tables.
    :type sim: dict(str,astropy.table.table.Table)
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.table.Table
    :returns: Hierarchical link table.
    :rtype: astropy.table.table.Table
    """
    sim_h_link = sim_helptab['main_id', 'parent_oid', 'h_link_ref', 'membership']
    #sim_h_link=nullvalues(sim_h_link,'parent_oid',0,verbose=False)
    ###sim_h_link=nullvalues(sim_h_link,'membership',-1,verbose=False)

    # removing entries in h_link where parent objects are clusters or 
    # associations as we are 
    #only interested in hierarchical multiples. 
    sim_h_link = sim_h_link[np.where(np.isin(sim_h_link['parent_oid'],
                                             stars['oid']))]

    sim_h_link = fetch_main_id(sim_h_link, OidCreator(name='parent_main_id', colname='parent_oid'))
    sim_h_link.remove_column('parent_oid')
    #typeconversion needed as smallint fill value != int null value
    sim_h_link['membership'] = sim_h_link['membership'].astype(int)
    sim_h_link = nullvalues(sim_h_link, 'membership', 999999)
    sim_h_link = replace_value(sim_h_link, 'h_link_ref', '',
                               sim['provider']['provider_bibcode'][0])
    sim_h_link = unique(sim_h_link)
    return sim_h_link


def create_objects_table(sim_helptab, stars):
    """
    Creates objects table.
    
    :param sim_helptab: Main SIMBAD helper table.
    :type sim_helptab: astropy.table.table.Table
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.table.Table
    :returns: Objects table.
    :rtype: astropy.table.table.Table
    """
    #-----------------creating output table sim_planets-----------------
    temp_sim_planets = sim_helptab['main_id', 'ids',
    'type'][np.where(sim_helptab['type'] == 'pl')]
    sim_planets = Table(unique(
        temp_sim_planets, keys='main_id'), copy=True)
    #-----------------creating output table sim_objects-----------------
    sim_objects = vstack([sim_planets['main_id', 'ids', 'type'],
                          stars['main_id', 'ids', 'type']])
    sim_objects['ids'] = sim_objects['ids'].astype(object)
    #tbd: add identifier simbad main_id without leading * and whitespaces
    return sim_objects


def create_sim_sources_table(stars, sim):
    """
    Creates sources table.
    
    :param sim: Dictionary of database table names and tables.
    :type sim: dict(str,astropy.table.table.Table)
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.table.Table
    :returns: Sources table.
    :rtype: astropy.table.table.Table
    """
    #--------------creating output table sim_sources -------------------
    tables = [sim['provider'], stars, sim['h_link'], sim['ident']]
    #define header name of columns containing references data
    ref_columns = [['provider_bibcode'], ['coo_ref', 'plx_ref', 'mag_i_ref',
                                          'mag_j_ref', 'mag_k_ref', 'binary_ref', 'sptype_ref'],
                   ['h_link_ref'], ['id_ref']]
    sim_sources = create_sources_table(tables, ref_columns,
                                       sim['provider']['provider_name'][0])
    return sim_sources


def create_star_basic_table(stars):
    """
    Creates basic stellar data table.
    
    :param stars: Secondary SIMBAD helper table.
    :type stars: astropy.table.table.Table
    :returns: Basic stellar data table.
    :rtype: astropy.table.table.Table
    """
    sim_star_basic = stars['main_id', 'coo_ra', 'coo_dec', 'coo_err_angle',
    'coo_err_maj', 'coo_err_min', 'coo_qual', 'coo_ref',
    'mag_i_value', 'mag_i_ref', 'mag_j_value', 'mag_j_ref',
    'mag_k_value', 'mag_k_ref',
    'sptype_string', 'sptype_qual', 'sptype_ref',
    'plx_value', 'plx_err', 'plx_qual', 'plx_ref']
    # changing type from object to string for later join functions
    sim_star_basic['sptype_string'] = sim_star_basic['sptype_string'].astype(str)
    sim_star_basic['sptype_qual'] = sim_star_basic['sptype_qual'].astype(str)
    sim_star_basic['sptype_ref'] = sim_star_basic['sptype_ref'].astype(str)
    return sim_star_basic


def provider_simbad(distance_cut_in_pc, test_objects=[]):
    """
    Optains and arranges SIMBAD data.
    
    :param float distance_cut_in_pc: Distance up to which stars are included.
    :param test_objects: Objects to be tested where they drop out of the criteria or make it till the end.
    :type test_objects: list(str)
    :returns: Dictionary with names and astropy tables containing
        reference data, provider data, object data, identifier data, object to 
        object relation data, basic stellar data and binarity data.
    :rtype: dict(str,astropy.table.table.Table)
    """

    sim_helptab, sim = create_simbad_helpertable(distance_cut_in_pc, test_objects)
    stars = creating_helpertable_stars(sim_helptab, sim)
    sim['ident'] = create_ident_table(sim_helptab, sim)
    sim['h_link'] = create_h_link_table(sim_helptab, sim, stars)
    stars = expanding_helpertable_stars(sim_helptab, sim, stars)
    sim['objects'] = create_objects_table(sim_helptab, stars)
    sim['sources'] = create_sim_sources_table(stars, sim)
    sim['star_basic'] = create_star_basic_table(stars)
    sim['mes_binary'] = stars['main_id', 'binary_flag', 'binary_qual', 'binary_ref']

    save(list(sim.values()), ['sim_' + element for element in list(sim.keys())])
    return sim
