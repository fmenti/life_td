"""
Generates the data for the LIFE Target Database.

"""


import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
import importlib #reloading external functions after modification

#self created modules
import sdc
import helperfunctions as hf
import provider as p
import building as b
importlib.reload(sdc)#reload module after changing it
importlib.reload(hf)#reload module after changing it
importlib.reload(p)#reload module after changing it
importlib.reload(b)#reload module after changing it


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
    
    sim=p.provider_simbad(empty_provider[:],distance_cut_in_pc)
    gk=p.provider_gk(table_names,empty_provider[:],distance_cut_in_pc)
    wds=p.provider_wds(table_names,empty_provider[:],False)
    exo=p.provider_exo(table_names,empty_provider[:],temp=False)
    life=p.provider_life(table_names,empty_provider[:])
    gaia=p.provider_gaia(table_names,empty_provider[:],distance_cut_in_pc)
    
    for i in range(len(sim)):
        sim[i]=hf.stringtoobject(sim[i])
        gk[i]=hf.stringtoobject(gk[i])
        wds[i]=hf.stringtoobject(wds[i])
        exo[i]=hf.stringtoobject(exo[i])
        life[i]=hf.stringtoobject(life[i])
        gaia[i]=hf.stringtoobject(gaia[i])

    #------------------------combine data from external sources---------
    database_tables=b.building([sim[:],gk[:],exo[:],life[:],gaia[:],wds[:]],table_names,empty_provider[:])
    return sim, gk, wds, exo, life, gaia, database_tables


def load_life_td():
    """
    Loads previously created life_td data.
    
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """    
    
    print(f'Loading life_td generated data with distance cut of {distance_cut_in_pc} pc')
    
    sim=hf.load(['sim_' + direction for direction in table_names])
    gk=hf.load(['gk_' + direction for direction in table_names])
    wds=hf.load(['wds_' + direction for direction in table_names])
    exo=hf.load(['exo_' + direction for direction in table_names])
    life=hf.load(['life_' + direction for direction in table_names])
    gaia=hf.load(['gaia_' + direction for direction in table_names])
    database_tables=hf.load(table_names,location='../../data/')
    
    for i in range(len(sim)):
        sim[i]=hf.stringtoobject(sim[i])
        gk[i]=hf.stringtoobject(gk[i])
        wds[i]=hf.stringtoobject(wds[i])
        exo[i]=hf.stringtoobject(exo[i])
        life[i]=hf.stringtoobject(life[i])
        gaia[i]=hf.stringtoobject(gaia[i])
        database_tables[i]=hf.stringtoobject(database_tables[i])
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
        sim=p.provider_simbad(empty_provider[:],distance_cut_in_pc)
    else:
        sim=hf.load(['sim_' + direction for direction in table_names])
    for i in range(len(sim)):
        sim[i]=hf.stringtoobject(sim[i])
    
    if 'gk' in create:
        gk=p.provider_gk(table_names,empty_provider[:],distance_cut_in_pc)
    else:
        gk=hf.load(['gk_' + direction for direction in table_names])
    for i in range(len(gk)):
        gk[i]=hf.stringtoobject(gk[i])
    
    if 'wds' in create:
        wds=p.provider_wds(table_names,empty_provider[:],False)
    else:
        wds=hf.load(['wds_' + direction for direction in table_names])
    for i in range(len(wds)):
        wds[i]=hf.stringtoobject(wds[i])
    
    if 'exo' in create:
        exo=p.provider_exo(table_names,empty_provider[:],temp=False)
    else:
        exo=hf.load(['exo_' + direction for direction in table_names])
    for i in range(len(exo)):
        exo[i]=hf.stringtoobject(exo[i])
    
    if 'life' in create:
        life=p.provider_life(table_names,empty_provider[:])
    else:
        life=hf.load(['life_' + direction for direction in table_names])
    for i in range(len(life)):
        life[i]=hf.stringtoobject(life[i])
    
    if 'gaia' in create:
        gaia=p.provider_gaia(table_names,empty_provider[:],distance_cut_in_pc)
    else:
        gaia=hf.load(['gaia_' + direction for direction in table_names])
    for i in range(len(gaia)):
        gaia[i]=hf.stringtoobject(gaia[i])
            
    #------------------------combine data from external sources---------
    database_tables=b.building([sim[:],gk[:],exo[:],life[:],gaia[:],wds[:]],table_names,empty_provider[:])
    return sim, gk, wds, exo, life, gaia, database_tables
    
if (__name__ == '__main__'):
    print('Executing as standalone script')
    create_life_td(5.)

