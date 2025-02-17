"""
Generates the data for the LIFE Target Database.

"""

from astropy.table import Table
from astropy import io

#self created modules
from sdata import empty_dict, empty_provider_tables_dict, table_names
from utils.io import stringtoobject, load, Path, string_to_object_whole_dict
from provider.exo import provider_exo
from provider.gaia import provider_gaia
from provider.life import provider_life
from provider.sdb import provider_sdb
from provider.simbad import provider_simbad
from provider.wds import provider_wds
from building import building


#tbd change provider to provider_info to minimize confusion in bulilding function with providers list
def load_cat(provider_name=''):
    cat=empty_dict.copy()
    if provider_name=='':
        prov=load([direction for direction in list(cat.keys())],location=Path().data)  
    else:
        prov=load([provider_name+'_' + direction for direction in list(cat.keys())])  
    for i,table in enumerate(list(cat.keys())):
        cat[table] = prov[i] 
    return cat


def load_life_td():
    """
    Loads previously created life_td data.
    
    :return: life_td data in different tables
    :rtype: list(astropy.table.table.Table)
    """    
    
    print(f'Loading life_td generated data')

    provider_tables_dict=empty_provider_tables_dict.copy()

    for i,prov in enumerate(list(provider_tables_dict.keys())):
        print(f'Loading {prov} data')
        cat = load_cat(prov)

    db = load_cat('')
    database_tables=string_to_object_whole_dict(db)
        
    return provider_tables_dict, database_tables

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

    # to do: solve error: 
    # ValueError: Unit 'MJup' not supported by the CDS standard. Did you mean Mjup, YMjup, ZMjup, yMjup or zMjup? same for 'ercent'
    with open(Path().data+'distance_cut.txt', 'w') as f:  
	# Write the floating point number to the file  
    	f.write(str(distance_cut_in_pc))

    
    print(f'Building life_td data with distance cut of {distance_cut_in_pc} pc')
    
    provider_tables_dict=empty_provider_tables_dict.copy()
    data=io.votable.parse_single_table(
        Path().additional_data+"sdb_30pc_09_02_2024.xml").to_table()
    functions = [provider_simbad,provider_sdb,provider_wds,provider_exo,provider_life,provider_gaia]
    arguments = [(distance_cut_in_pc),(distance_cut_in_pc,data),(False),(),(),(distance_cut_in_pc)]
    
    for i,prov in enumerate(list(provider_tables_dict.keys())):
        #provider_tables_dict[prov]=string_to_object_whole_dict(provider_tables_dict[prov])
        if prov in create:
            if arguments[i]!=():
                try:
                    cat=functions[i](arguments[i])
                except:
                    cat=functions[i](*arguments[i])
            else:
                cat=functions[i]()
        else:
            print(f'Loading {prov} data')
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
    
    #------------------------obtain data from external sources---------------------
    provider_tables_dict, database_tables= partial_create(distance_cut_in_pc,
                                   create=['sim', 'sdb', 'wds', 'exo', 'life', 'gaia'])
    return provider_tables_dict, database_tables
    
if (__name__ == '__main__'):
    print('Executing as standalone script')
    create_life_td(5.)

