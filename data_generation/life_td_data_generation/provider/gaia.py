""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from pyvo.dal import TAPService
from astropy.table import Table, vstack, setdiff, join
from datetime import datetime

#self created modules
from utils.io import save
from provider.utils import fetch_main_id, IdentifierCreator, sources_table, ids_from_ident, replace_value
import sdata as sdc


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
    gaia_provider=Table()
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
        service = TAPService(gaia_provider['provider_url'][0])
        result=service.run_sync(adql_query.format(**locals()), maxrec=160000)
        gaia=result.to_table()

    else:
        gaia=query(gaia_provider['provider_url'][0],adql_query) 
        
    gaia.rename_columns(['mass_flame','radius_flame'],
                        ['mass_st_value','radius_st_value'])
    gaia['gaia_id']=['Gaia DR3 '+str(gaia['source_id'][j]) for j in range(len(gaia))]
    gaia['ref']=['2022arXiv220800211G' for j in range(len(gaia))]#dr3 paper
    
    #---------------gaia_ident-----------------------
    gaia_sim_idmatch=fetch_main_id(gaia['gaia_id','ref'],
                           IdentifierCreator(name='main_id',colname='gaia_id')) 
    #should be gaia_id, main_id, ref minus 40 objects that have only gaia_id
    gaia_ident=gaia_sim_idmatch.copy()
    gaia_ident.rename_columns(['gaia_id','ref'],['id','id_ref'])
    gaia_ident['id_ref']=gaia_ident['id_ref'].astype(str)
    #creating simbad main_id ident rows
    sim_main_id_ident=Table()
    sim_main_id_ident['main_id']=gaia_ident['main_id']
    sim_main_id_ident['id']=gaia_ident['main_id']
    sim_main_id_ident['id_ref']=['2000A&AS..143....9W' for j in range(len(gaia_ident))]
    gaia_ident=vstack([gaia_ident,sim_main_id_ident])
    #now need to add the 40 objects that have only gaia_identifiers
    #for setdiff need both columns to be same type
    for col in gaia_sim_idmatch.colnames:
        gaia_sim_idmatch[col]=gaia_sim_idmatch[col].astype(str)
    gaia_only_id=setdiff(gaia['gaia_id','ref'],gaia_sim_idmatch['gaia_id','ref'])
    gaia_only_id['main_id']=gaia_only_id['gaia_id']
    gaia_only_id.rename_columns(['gaia_id','ref'],['id','id_ref'])
    #for vstack need both columns to be same type
    for col in gaia_ident.colnames:
        gaia_ident[col]=gaia_ident[col].astype(str)
    gaia_ident=vstack([gaia_ident,gaia_only_id])
    #add main_id to gaia table
    gaia=join(gaia_ident['main_id','id'],gaia,
                       keys_left='id', keys_right='gaia_id')
    gaia.remove_column('id')

    #-----------------gaia_objects------------------
    gaia_objects=Table(names=['main_id','ids'],dtype=[object,object])
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
    gaia_objects=join(gaia_objects,gaia['main_id','nss_solution_type'],join_type='left')
    gaia_objects['type'][np.where(gaia_objects['nss_solution_type']!='')]=['sy' for j in range(len(
            gaia_objects['type'][np.where(gaia_objects['nss_solution_type']!='')]))]
    gaia_objects.remove_column('nss_solution_type')

    #gaia_mes_binary
    gaia_mes_binary=gaia_objects['main_id','type']
    # tbd add binary flag True to children of system objects once I get h_link 
    # info from gaia')
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
    
    gaia_mes_teff_st=vstack([gaia_mes_teff_st,temp])
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
    gaia_sources=Table()
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
        save([gaia_list_of_tables[i][:]],['gaia_'+table_names[i]])
    return gaia_list_of_tables
