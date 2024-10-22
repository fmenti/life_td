"""
Generates the data for the LIFE Target Database.

"""

from astropy.table import Table

#self created modules
from sdata import empty_cat, empty_provider_tables_dict
from utils.io import stringtoobject, load, Path, string_to_object_whole_dict
from provider.exo import provider_exo
from provider.gaia import provider_gaia
from provider.life import provider_life
from provider.sdb import provider_sdb
from provider.simbad import provider_simbad
from provider.wds import provider_wds
from building import building


#tbd change provider to provider_info to minimize confusion in bulilding function with providers list

def load_life_td():
    """
    Loads previously created life_td data.
    
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """    
    
    print(f'Loading life_td generated data')
    
    sim=load(['sim_' + direction for direction in table_names])
    sdb=load(['sdb_' + direction for direction in table_names])
    wds=load(['wds_' + direction for direction in table_names])
    exo=load(['exo_' + direction for direction in table_names])
    life=load(['life_' + direction for direction in table_names])
    gaia=load(['gaia_' + direction for direction in table_names])
    database_tables=load(table_names,location=Path().data)
    
    for i in range(len(sim)):
        sim[i]=stringtoobject(sim[i])
        sdb[i]=stringtoobject(sdb[i])
        wds[i]=stringtoobject(wds[i])
        exo[i]=stringtoobject(exo[i])
        life[i]=stringtoobject(life[i])
        gaia[i]=stringtoobject(gaia[i])
        database_tables[i]=stringtoobject(database_tables[i])
    return sim, sdb, wds, exo, life, gaia, database_tables

def load_cat(provider_name):
    cat=empty_cat.copy()
    prov=load([provider_name+'_' + direction for direction in list(cat.keys())])  
    for i,table in enumerate(list(cat.keys())):
        cat[table] = prov[i] 
    return cat

def partial_create(distance_cut_in_pc,create=[]):
    """
    Partially generates, partially loads life_td data.
    
    Generates the in the create list specified life_td data, loads the 
    rest and builds everything together.
    
    :param distance_cut_in_pc: Distance cut of the objects in parsec.
    :type distance_cut_in_pc: float
    :param create: If one or more of 'sim', 'sdb', 'wds', 'exo', 'life',
        and 'gaia'  are present, generates those tables, the missing 
        ones are loaded.
    :type create: list(str)
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """
    
    print(f'Building life_td data with distance cut of {distance_cut_in_pc} pc')
    
    provider_tables_dict=empty_provider_tables_dict.copy()
      
    functions = [provider_simbad,provider_sdb,provider_wds,provider_exo,provider_life,provider_gaia]
    arguments = [(distance_cut_in_pc),(distance_cut_in_pc),(False),(True),(),(distance_cut_in_pc)]
    
    for i,prov in enumerate(list(provider_tables_dict.keys())):
        #provider_tables_dict[prov]=string_to_object_whole_dict(provider_tables_dict[prov])
        if prov in create:
            if arguments[i]!=():
                cat=functions[i](arguments[i])
            else:
                cat=functions[i]()
        else:
            cat = load_cat(prov)
        provider_tables_dict[prov]=string_to_object_whole_dict(cat)
    
            
    #------------------------combine data from external sources---------
    database_tables=building(provider_tables_dict)
    return provider_tables_dict, database_tables

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
    provider_tables_dict, database_tables= partial_create(distance_cut_in_pc,
                                   create=['sim', 'sdb', 'wds', 'exo', 'life', 'gaia'])
    return provider_tables_dict, database_tables
    
if (__name__ == '__main__'):
    print('Executing as standalone script')
    create_life_td(5.)

