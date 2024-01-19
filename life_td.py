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
              'mes_teff_st','mes_radius_st','mes_mass_st','mes_binary','mes_sep_ang']


def create_life_td(distance_cut_in_pc):
    
    #------------------------obtain data from external sources---------------------
    empty_provider=[ap.table.Table() for i in range(len(table_names))]
    
    sim=p.provider_simbad(table_names,empty_provider[:],distance_cut_in_pc)
    wds=p.provider_wds(table_names,empty_provider[:],distance_cut_in_pc,False)
    gk=p.provider_gk(table_names,empty_provider[:],distance_cut_in_pc)
    exo=p.provider_exo(table_names,empty_provider[:],distance_cut_in_pc,temp=False)
    life=p.provider_life(table_names,empty_provider[:],distance_cut_in_pc)
    gaia=p.provider_gaia(table_names,empty_provider[:],distance_cut_in_pc)

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