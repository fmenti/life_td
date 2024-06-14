"""
Generates the data for the LIFE Target Database.

"""


import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables

#self created modules
import sdata as sdc
from utils.utils import stringtoobject, load
from provider.exo import provider_exo
from provider.gaia import provider_gaia
from provider.life import provider_life
from provider.sdb import provider_gk
from provider.simbad import provider_simbad
from provider.wds import provider_wds
from building import building

data_path='../../data/'

empty=sdc.provider('empty')
table_names=empty.table_names
#tbd change provider to provider_info to minimize confusion in bulilding function with providers list


def create_life_td(distance_cut_in_pc):
    """
    Generates the life_td data.
    
    Calls the modules provider, helperfunctions and building.
    
    :param distance_cut_in_pc: Distance cut of the objects in parsec.
    :type distance_cut_in_pc: float
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """
    
    print(f'Generating life_td data with distance cut of {distance_cut_in_pc} pc')
    #------------------------obtain data from external sources---------------------
    empty_provider=[ap.table.Table() for i in range(len(table_names))]
    
    sim=provider_simbad(empty_provider[:],distance_cut_in_pc)
    gk=provider_gk(table_names,empty_provider[:],distance_cut_in_pc)
    wds=provider_wds(table_names,empty_provider[:],False)
    exo=provider_exo(table_names,empty_provider[:],temp=False)
    life=provider_life(table_names,empty_provider[:])
    gaia=provider_gaia(table_names,empty_provider[:],distance_cut_in_pc)
    
    for i in range(len(sim)):
        sim[i]=stringtoobject(sim[i])
        gk[i]=stringtoobject(gk[i])
        wds[i]=stringtoobject(wds[i])
        exo[i]=stringtoobject(exo[i])
        life[i]=stringtoobject(life[i])
        gaia[i]=stringtoobject(gaia[i])

    #------------------------combine data from external sources---------
    database_tables=building([sim[:],gk[:],exo[:],life[:],gaia[:],wds[:]],table_names,empty_provider[:])
    return sim, gk, wds, exo, life, gaia, database_tables


def load_life_td():
    """
    Loads previously created life_td data.
    
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """    
    
    print(f'Loading life_td generated data')
    
    sim=load(['sim_' + direction for direction in table_names])
    gk=load(['gk_' + direction for direction in table_names])
    wds=load(['wds_' + direction for direction in table_names])
    exo=load(['exo_' + direction for direction in table_names])
    life=load(['life_' + direction for direction in table_names])
    gaia=load(['gaia_' + direction for direction in table_names])
    database_tables=load(table_names,location=data_path)
    
    for i in range(len(sim)):
        sim[i]=stringtoobject(sim[i])
        gk[i]=stringtoobject(gk[i])
        wds[i]=stringtoobject(wds[i])
        exo[i]=stringtoobject(exo[i])
        life[i]=stringtoobject(life[i])
        gaia[i]=stringtoobject(gaia[i])
        database_tables[i]=stringtoobject(database_tables[i])
    return sim, gk, wds, exo, life, gaia, database_tables

def partial_create(distance_cut_in_pc,create=['sim', 'gk', 'wds', 'exo', 'life', 'gaia']):
    """
    Partially generates, partially loads life_td data.
    
    Generates the in the create list specified life_td data, loads the 
    rest and builds everything together.
    
    :param distance_cut_in_pc: Distance cut of the objects in parsec.
    :type distance_cut_in_pc: float
    :param create: If one or more of 'sim', 'gk', 'wds', 'exo', 'life',
        and 'gaia'  are present, generates those tables, the missing 
        ones are loaded.
    :type create: list(str)
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """
    
    print(f'Building life_td data with distance cut of {distance_cut_in_pc} pc')
    
    empty_provider=[ap.table.Table() for i in range(len(table_names))]
    
    if 'sim' in create:
        sim=provider_simbad(empty_provider[:],distance_cut_in_pc)
    else:
        sim=load(['sim_' + direction for direction in table_names])
    for i in range(len(sim)):
        sim[i]=stringtoobject(sim[i])
    
    if 'gk' in create:
        gk=provider_gk(table_names,empty_provider[:],distance_cut_in_pc)
    else:
        gk=load(['gk_' + direction for direction in table_names])
    for i in range(len(gk)):
        gk[i]=stringtoobject(gk[i])
    
    if 'wds' in create:
        wds=provider_wds(table_names,empty_provider[:],False)
    else:
        wds=load(['wds_' + direction for direction in table_names])
    for i in range(len(wds)):
        wds[i]=stringtoobject(wds[i])
    
    if 'exo' in create:
        exo=provider_exo(table_names,empty_provider[:],temp=False)
    else:
        exo=load(['exo_' + direction for direction in table_names])
    for i in range(len(exo)):
        exo[i]=stringtoobject(exo[i])
    
    if 'life' in create:
        life=provider_life(table_names,empty_provider[:])
    else:
        life=load(['life_' + direction for direction in table_names])
    for i in range(len(life)):
        life[i]=stringtoobject(life[i])
    
    if 'gaia' in create:
        gaia=provider_gaia(table_names,empty_provider[:],distance_cut_in_pc)
    else:
        gaia=load(['gaia_' + direction for direction in table_names])
    for i in range(len(gaia)):
        gaia[i]=stringtoobject(gaia[i])
            
    #------------------------combine data from external sources---------
    database_tables=building([sim[:],gk[:],exo[:],life[:],gaia[:],wds[:]],table_names,empty_provider[:])
    return sim, gk, wds, exo, life, gaia, database_tables
    
if (__name__ == '__main__'):
    print('Executing as standalone script')
    create_life_td(5.)

