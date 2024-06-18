""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np #arrays
from astropy import units, io, coordinates
from astropy.table import Table, unique, join, MaskedColumn, Column
from datetime import datetime

#self created modules
from utils.utils import save, load
from provider.utils import sources_table, replace_value
import sdata as sdc

additional_data_path='../../additional_data/'


def lum_class(nr,sptype):
    """
    Extracts luminocity class.
    
    :param int nr: index number
    :param str sptype: spectral type
    :returns: luminocity class
    :rtype: str
    """
    
    lum_class=sptype[nr]
    if len(sptype)>nr+1 and sptype[nr+1] in ['I','V']:
        lum_class=sptype[nr:nr+2]
        if len(sptype)>nr+2 and sptype[nr+2] in ['I','V']:
            lum_class=sptype[nr:nr+3]
    return lum_class

def sptype_string_to_class(temp,ref):
    """
    Extracts stellar parameters from spectral type string one.

    This function extracts the temperature class, temperature class number
    and luminocity class information from the spectral type string (e.g. 
    M5V to M, 5 and V). It stores that information in the for this purpose
    generated new columns. Only objects of temperature class O, B, A, F,
    G, K, and M are processed. Only objects of luminocity class IV, V and VI
    are processed.

    :param temp: Table containing spectral type information in
        the column sptype_string.
    :type temp: astropy.table.table.Table
    :param str ref: Designates origin of data.
    :returns: Table like temp with additional columns class_temp,
        class_temp_nr, class_lum and class_ref.
    :rtype: astropy.table.table.Table
    """

    temp['class_temp']=MaskedColumn(dtype=object,length=len(temp))
    temp['class_temp_nr']=MaskedColumn(dtype=object,length=len(temp))
    temp['class_lum']=MaskedColumn(dtype=object,length=len(temp))
    temp['class_ref']=MaskedColumn(dtype=object,length=len(temp))
    #tbd: rewrite code using recoursive function

    for i in range(len(temp)):
        #sorting out objects like M5V+K7V
        #strip d for spectral types starting with small d because it is an old annotation for dwarf star
        if len(temp['sptype_string'][i])>0:
            if temp['sptype_string'][i][0]=='d':
                temp['class_lum'][i]='V'
        sptype=temp['sptype_string'][i].strip('d')
        
        if (len(sptype.split('+'))==1 and
                #sorting out entries like ''
                len(sptype)>0 and 
                # sorting out brown dwarfs i.e. T1V
                sptype[0] in ['O','B','A','F','G','K','M']):
            #assigning temperature class and reference
            temp['class_temp'][i]=sptype[0]
            temp['class_ref'][i]=ref
            if len(sptype)>1 and sptype[1] in ['0','1','2','3','4','5','6','7','8','9']:
                temp['class_temp_nr'][i]=sptype[1]
                #distinguishing between objects like K5V and K5.5V
                if len(sptype)>2 and sptype[2]=='.':
                    temp['class_temp_nr'][i]=sptype[1:4]
                    if len(sptype)>4 and sptype[4] in ['I','V']:
                        temp['class_lum'][i]=lum_class(4,sptype)
                    #make sure sptypes starting with d don't get class_lum overwritten
                    else:
                        temp['class_lum'][i]='V'
                elif len(sptype)>2 and sptype[2] in ['I','V']:
                    temp['class_lum'][i]=lum_class(2,sptype)
                #tbd add assumption of V if nothing given. valid because V is longest 
                #evolution stage so most stars will be in V
                else:
                    temp['class_lum'][i]='V'
            else:
                temp['class_lum'][i]='V'
        else:
            temp['class_temp'][i]='?'
            temp['class_temp_nr'][i]='?'
            temp['class_lum'][i]='?'
            temp['class_ref'][i]='?'
    return temp

def realspectype(cat):
    """
    Removes rows not containing main sequence stars.

    Removes rows of cat where elements in column named 'sim_sptype' are
    either '', 'nan' or start with an other letter than the main sequence
    spectral type classification.

    :param cat: Table containing 'sim_sptype' column
    :type cat: astropy.table.table.Table
    :returns: Table, param cat with undesired rows removed
    :rtype: astropy.table.table.Table
    """
    ms_tempclass=np.array(['O','B','A','F','G','K','M'])
    ms_temp=cat[np.where(np.in1d(cat['class_temp'],ms_tempclass))]
    
    ms_lumclass=np.array(['V'])
    ms=ms_temp[np.where(np.in1d(ms_temp['class_lum'],ms_lumclass))]
    
    return ms

def model_param():
    """
    Loads and cleans up model file.

    Loads the table of Eric E. Mamajek containing stellar parameters 
    modeled from spectral types. Cleans up the columns for spectral 
    type, effective temperature radius and mass.

    :returns: Table of the 4 parameters as columns
    :rtype: astropy.table.table.Table
    """

    EEM_table=io.ascii.read(additional_data_path+"Mamajek2022-04-16.csv")['SpT','Teff','R_Rsun','Msun']
    EEM_table.rename_columns(['R_Rsun','Msun'],['Radius','Mass'])
    EEM_table=replace_value(EEM_table,'Radius',' ...','nan')
    EEM_table=replace_value(EEM_table,'Mass',' ...','nan')
    EEM_table=replace_value(EEM_table,'Mass',' ....','nan')
    EEM_table['Teff'].unit=units.K
    EEM_table['Radius'].unit=units.R_sun
    EEM_table['Mass'].unit=units.M_sun       
    io.votable.writeto(io.votable.from_table(EEM_table), \
                          f'{additional_data_path}model_param.xml')#saving votable
    return EEM_table

def match_sptype(cat,model_param,sptypestring='sim_sptype',teffstring='mod_Teff',\
                 rstring='mod_R',mstring='mod_M'):
    """
    Assigns modeled parameter values.

    Matches the spectral types with the ones in Mamajek's table and 
    includes the modeled effective Temperature,
    stellar radius and stellar mass into the catalog.

    :param cat: astropy table containing spectral type information
    :type cat: astropy.table.table.Table
    :param str sptypestring: Column name where the spectral 
        type information is located
    :param str teffstring: Column name for the stellar effective 
        temperature column
    :param str rstring: Column name for the stellar radius column
    :param str mstring: Column name for the stellar mass column
    :returns: Table cat with added new columns for 
        effective temperature, radius and mass filled with model values
    :rtype: astropy.table.table.Table
    """

    #initiating columns with right units
    
    arr=np.zeros(len(cat))
    cat[teffstring]=arr*np.nan*units.K
    cat[teffstring]=MaskedColumn(mask=np.full(len(cat),True), \
                                         length=len(cat),unit=units.K)
    cat[rstring]=arr*np.nan*units.R_sun
    cat[mstring]=arr*np.nan*units.M_sun
    #go through all spectral types in cat
    for j in range(len(cat[sptypestring])): 
        # for all the entries that are not empty
        if cat[sptypestring][j]!='':
            #remove first d coming from old notation for dwarf meaning main sequence star
            cat[sptypestring][j]=cat[sptypestring][j].strip('d')
            #go through the model spectral types of Mamajek 
            for i in range(len(model_param['SpT'])): 
                #match first two letters
                if model_param['SpT'][i][:2]==cat[sptypestring][j][:2]: 
                        cat[teffstring][j]=model_param['Teff'][i]
                        cat[rstring][j]=model_param['Radius'][i]
                        cat[mstring][j]=model_param['Mass'][i]
            #as the model does not cover all spectral types on .5 accuracy, check those separately
            if cat[sptypestring][j][2:4]=='.5':
                for i in range(len(model_param['SpT'])):
                    # match first four letters
                    if model_param['SpT'][i][:4]==cat[sptypestring][j][:4]:
                        cat[teffstring][j]=model_param['Teff'][i]
                        cat[rstring][j]=model_param['Radius'][i]
                        cat[mstring][j]=model_param['Mass'][i] 
        else:
            cat[sptypestring][j]='None' 
    return cat

def spec(cat):
    """
    Runs the spectral type related functions realspectype and match_sptype. 

    It also removes all empty columns of the effective temperature, removes 
    rows that are not main sequence, removes rows with binary subtype and 
    non unique simbad name.

    :param cat: astropy table containing columns named 
        'sim_sptype','sim_name' and 'sim_otypes'
    :type cat: astropy.table.table.Table
    :returns: Catalog of mainsequence stars with unique 
        simbad names, no binary subtypes and modeled parameters.
    :rtype: astropy.table.table.Table
    """   

    #Do I even need realspectype function? I can just take cat where class_temp not empty
    cat=realspectype(cat)
    #model_param=io.votable.parse_single_table(\
        #f"catalogs/model_param.xml").to_table()
    mp=model_param()#create model table as votable
    cat=match_sptype(cat,mp,sptypestring='sptype_string')
    cat.remove_rows([np.where(cat['mod_Teff'].mask==True)])
    cat.remove_rows([np.where(np.isnan(cat['mod_Teff']))])
    cat=unique(cat, keys='main_id')
    return cat

def provider_life(table_names,life_list_of_tables):
    """
    Loads SIMBAD data and postprocesses it. 
    
    Postprocessing enables to provide more useful information. It uses a model
    from Eric E. Mamajek to predict temperature, mass and radius from the simbad 
    spectral type data.
    
    :param table_names: Contains the names for the output tables.
    :type table_names: list(str)
    :param life_list_of_tables: List of same length as table_names containing
        empty astropy tables.
    :type life_list_of_tables: list(astropy.table.table.Table)
    :returns: List of astropy table containing
        reference data, provider data, basic stellar data, stellar effective
        temperature data, stellar radius data and stellar mass data.
    :rtype: list(astropy.table.table.Table)
    """
    
    #---------------define provider-------------------------------------
    life_provider=Table()
    life_provider['provider_name']=['LIFE']
    life_provider['provider_url']=['www.life-space-mission.com']
    life_provider['provider_bibcode']=['2022A&A...664A..21Q']
    life_provider['provider_access']=datetime.now().strftime('%Y-%m-%d')
    
    print('Creating ',life_provider['provider_name'][0],' tables ...')
    #---------------------star_basic----------------
    #galactic coordinates:  transformed from simbad ircs coordinates using astropy
    [life_star_basic]=load(['sim_star_basic'])
    ircs_coord=coordinates.SkyCoord(\
            ra=life_star_basic['coo_ra'],dec=life_star_basic['coo_dec'],frame='icrs')
    gal_coord=ircs_coord.galactic
    life_star_basic['coo_gal_l']=gal_coord.l.deg*units.degree
    life_star_basic['coo_gal_b']=gal_coord.b.deg*units.degree
    life_star_basic['dist_st_value']=1000./life_star_basic['plx_value'] 
    life_star_basic['dist_st_value']=np.round(life_star_basic['dist_st_value'],2)
    #null value treatment: plx_value has masked entries therefore distance_values too
    #ref:
    life_star_basic['dist_st_ref']=MaskedColumn(dtype=object,length=len(life_star_basic),
                                    mask=[True for j in range(len(life_star_basic))])
    life_star_basic['dist_st_ref'][np.where(life_star_basic['dist_st_value'].mask==False)]= \
            [life_provider['provider_name'][0] for j in range(len(
                life_star_basic['dist_st_ref'][np.where(life_star_basic['dist_st_value'].mask==False)]))]
    # can I do the same transformation with the errors? -> try on some examples and compare to simbad ones
    life_star_basic['coo_gal_err_angle']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_err_maj']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_err_min']=[-1
                        for j in range(len(life_star_basic))]
    life_star_basic['coo_gal_qual']=['?'
                        for j in range(len(life_star_basic))]
    life_star_basic['main_id']=life_star_basic['main_id'].astype(str)
    # source
    # transformed from simbad ircs coordinates using astropy
    life_star_basic['coo_gal_ref']=Column(dtype=object,length=len(life_star_basic))
    life_star_basic['coo_gal_ref']=life_provider['provider_name'][0] 
    #for all entries since coo_gal column not masked column
             
    life_star_basic['coo_gal_ref']=life_star_basic['coo_gal_ref'].astype(str)
    life_star_basic=life_star_basic['main_id','coo_gal_l','coo_gal_b','coo_gal_err_angle',
                                   'coo_gal_err_maj','coo_gal_err_min','coo_gal_qual',
                                   'coo_gal_ref','dist_st_value','dist_st_ref','sptype_string']
    

    life_star_basic=sptype_string_to_class(life_star_basic,life_provider['provider_name'][0])
    
    #-----------measurement tables -------------------------------------
    #applying model from E. E. Mamajek on SIMBAD spectral type

    [sim_objects]=load(['sim_objects'],stringtoobjects=False)
    
    stars=sim_objects[np.where(sim_objects['type']=='st')]
    cat=join(stars,life_star_basic)
    cat=spec(cat['main_id','sptype_string','class_lum','class_temp'])
    #if I take only st objects from sim_star_basic I don't loose objects during realspectype
    life_mes_teff_st=cat['main_id','mod_Teff']
    life_mes_teff_st.rename_column('mod_Teff','teff_st_value')
    life_mes_teff_st['teff_st_qual']=['D' for i in range(len(life_mes_teff_st))]
    life_mes_teff_st['teff_st_ref']=['2013ApJS..208....9P' for i in range(len(life_mes_teff_st))]
    
    life_mes_radius_st=cat['main_id','mod_R']
    life_mes_radius_st.rename_column('mod_R','radius_st_value')
    life_mes_radius_st['radius_st_qual']=['D' for i in range(len(life_mes_radius_st))]
    life_mes_radius_st['radius_st_ref']=['2013ApJS..208....9P' for i in range(len(life_mes_radius_st))]
    
    life_mes_mass_st=cat['main_id','mod_M']
    life_mes_mass_st.rename_column('mod_M','mass_st_value')
    life_mes_mass_st['mass_st_qual']=['D' for i in range(len(life_mes_mass_st))]
    life_mes_mass_st['mass_st_ref']=['2013ApJS..208....9P' for i in range(len(life_mes_mass_st))]
    
    #specifying stars cocerning multiplicity
    #main sequence simbad object type: MS*, MS? -> luminocity class
    #Interacting binaries and close CPM systems: **, **?
    
    #-----------------sources table-------------------------------------
    life_sources=Table()
    tables=[life_provider,life_star_basic,life_mes_teff_st,life_mes_radius_st,life_mes_mass_st]
    ref_columns=[['provider_bibcode'],['coo_gal_ref'],['teff_st_ref'],['radius_st_ref'],['mass_st_ref']]
    for cat,ref in zip(tables,ref_columns):
        life_sources=sources_table(cat,ref,life_provider['provider_name'][0],life_sources)
    
    #removing this column because I had to adapt it where there was a leadin d entry but change not useful for db just for 
    #life parameter creation
    life_star_basic.remove_column('sptype_string')
    
    for i in range(len(table_names)):
        if table_names[i]=='sources': life_list_of_tables[i]=life_sources
        if table_names[i]=='provider': life_list_of_tables[i]=life_provider
        if table_names[i]=='star_basic': life_list_of_tables[i]=life_star_basic
        if table_names[i]=='mes_teff_st': life_list_of_tables[i]=life_mes_teff_st
        if table_names[i]=='mes_radius_st': life_list_of_tables[i]=life_mes_radius_st
        if table_names[i]=='mes_mass_st': life_list_of_tables[i]=life_mes_mass_st
        save([life_list_of_tables[i][:]],['life_'+table_names[i]])
    return life_list_of_tables
