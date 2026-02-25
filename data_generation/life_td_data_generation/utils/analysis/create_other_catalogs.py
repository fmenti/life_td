from typing import Any

import numpy as np #arrays
import astropy as ap #votables
import importlib #reloading external functions after modification
import sys

sys.path.append('../../implementation/life/data_generation/life_td_data_generation')

#self created modules
from utils.io import save, load, stringtoobject, Path
import utils.analysis.analysis as la
importlib.reload(la)#reload module after changing it
import provider as p
importlib.reload(p)#reload module after changing it

para3=['sim_ra', 'sim_dec', 'sim_plx', 'distance', 'gal_coord_l', 'gal_coord_b',
       'mod_Teff', 'mod_R', 'mod_M', 'sep_phys',  'sim_i', 'sim_j']
para4=['coo_ra', 'coo_dec', 'plx_value', 'dist_st_value', 'coo_gal_l', 'coo_gal_b',
       'teff_st_value', 'radius_st_value','mass_st_value','sep_phys_value',  'mag_i_value', 'mag_j_value']

para_hwo=['ra', 'dec', 'sy_plx', 'sy_dist', 'st_teff', 'st_rad', 'st_mass']
para_starcat4=['coo_ra', 'coo_dec', 'plx_value', 'dist_st_value', 'teff_st_value',
               'radius_st_value','mass_st_value']


def create_updated_ltc3():
    """
    This function updates the StarCat3 with newer main_id's.

    :returns: LIFE-StarCat3 with updated main_id from SIMBAD.
    :rtype: astropy.table.table.Table
    """

    #load LTC3
    ltc3=ap.io.ascii.read("data/LIFE-StarCat3.csv")
    ltc3=stringtoobject(ltc3,3000) #len 1732

    ltc3_update=p.fetch_main_id(ltc3,'sim_name',name='main_id',oid=False)
    #print(ltc3_update)#1728 -> lost 4 objects
    lost=ltc3[np.where(np.invert(np.in1d(ltc3['sim_name'],ltc3_update['sim_name'])))]
    #L  326-61, L  380-78,  L  755-19, LP  399-68 -> maybe were truncated?

    print('To compare the version 3 to version 4 we had to updated the simbad main identifiers of v3.')
    print('Of the originally ',len(ltc3), 'objects ',len(lost),'could not be matched that way.')
    print('After matching on position those could be retrieved again.')
    newltc3=ltc3
    newltc3['sim_name'][np.where(newltc3['sim_name']=='L  326-61')]='L 326-21'
    newltc3['sim_name'][np.where(newltc3['sim_name']=='L  380-78')]='L 380-8'
    newltc3['sim_name'][np.where(newltc3['sim_name']=='L  755-19')]='LP 755-19'
    newltc3['sim_name'][np.where(newltc3['sim_name']=='LP  399-68')]='L 399-68'
    newltc3=p.fetch_main_id(newltc3,'sim_name',name='main_id',oid=False)
    print('newltc3',newltc3)
    save([newltc3],['newltc3'],location='data/')
    [StarCat4]=load(['StarCat4'],location='data/')
    in_3_but_not_in_4=newltc3[np.where(np.invert(np.in1d(newltc3['main_id'],StarCat4['main_id'])))]
    print('in_3_but_not_in_4',in_3_but_not_in_4)
    save([in_3_but_not_in_4],['in_3_but_not_in_4'],location='data/')
    return newltc3


def create_hwo():
    """
    Prepares hwo catalog for comparison to LIFE-StarCat4.

    Reads the HWO SPORES catalog file, adds SIMBAD main identifiers and transforms
    StarCat4 comparable columns from type string to float.

    :returns: Table containing HWO catalog data.
    :rtype: astropy.table.table.Table
    """

    spores_raw=ap.io.ascii.read("data/spores_catalog_v1.2.0.csv")#,header_start=0,)
    spores_raw=stringtoobject(spores_raw,3000)
    colnames=list(spores_raw[0])
    colnames=colnames[1:]
    data=spores_raw[1:].copy()

    data.remove_columns('col1')
    spores=ap.table.Table(data,names=colnames)
    print(spores)
    #add main_id column
    print(spores.colnames)
    temp=spores['tic_id','hip_name'].copy()
    hwo=p.fetch_main_id(temp,'tic_id',name='main_id',oid=False)
    hwo.rename_column('hip_name','temp')
    hwo=ap.table.join(spores,hwo,keys='tic_id')
    hwo.remove_column('temp')
    print(hwo)
    for para in para_hwo:
            hwo[para]=hwo[para].astype(float)
    save([hwo],['hwo'],location='data/')
    return hwo

def create_hpic():
    HPIC= ap.io.ascii.read('data/full_HPIC.csv')
    print(HPIC)
    print(HPIC.colnames)
    HPIC=HPIC[np.where(HPIC['sy_dist']!='null')]#143 null values in distance
    HPIC=HPIC[np.where(HPIC['st_spectype']!='null')]#of those about 2600 have no spectral type given
    HPIC['sy_dist']=HPIC['sy_dist'].astype(float)
    save([HPIC],['hpic'],location='data/')
    return HPIC

