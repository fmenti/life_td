import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
from datetime import datetime
import importlib #reloading external functions after modification

#self created modules
import helperfunctions as hf
import provider as p
import building as b
importlib.reload(hf)#reload module after changing it
importlib.reload(p)#reload module after changing it
importlib.reload(b)#reload module after changing it



table_names=['sources','objects','provider','ident','h_link','star_basic',
              'planet_basic','disk_basic','mes_mass_pl',
              'mes_teff_st','mes_radius_st','mes_mass_st','mes_binary',
             'mes_sep_ang','best_h_link']
#tbd change provider to provider_info to minimize confusion in bulilding function with providers list


def create_life_td(distance_cut_in_pc):
    
    #------------------------obtain data from external sources---------------------
    empty_provider=[ap.table.Table() for i in range(len(table_names))]
    
    sim=p.provider_simbad(table_names,empty_provider[:],distance_cut_in_pc)
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

    #------------------------combine data from external sources-------------------
    database_tables=b.building([sim[:],gk[:],exo[:],life[:],gaia[:],wds[:]],table_names,empty_provider[:])
    return sim, gk, wds, exo, life, gaia, database_tables


def load_life_td():
    #loading instead of running data:
    
    
    sim=hf.load(['sim_' + direction for direction in table_names])
    gk=hf.load(['gk_' + direction for direction in table_names])
    wds=hf.load(['wds_' + direction for direction in table_names])
    exo=hf.load(['exo_' + direction for direction in table_names])
    life=hf.load(['life_' + direction for direction in table_names])
    gaia=hf.load(['gaia_' + direction for direction in table_names])
    database_tables=hf.load(table_names)
    
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
    empty_provider=[ap.table.Table() for i in range(len(table_names))]
    
    if 'sim' in create:
        sim=p.provider_simbad(table_names,empty_provider[:],distance_cut_in_pc)
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
            
    #------------------------combine data from external sources-------------------
    database_tables=b.building([sim[:],gk[:],exo[:],life[:],gaia[:],wds[:]],table_names,empty_provider[:])
    return sim, gk, wds, exo, life, gaia, database_tables