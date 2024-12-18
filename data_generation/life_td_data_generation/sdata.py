"""
Dictionaries
"""
from astropy.table import Table

paras_dict = {
    'sources' : [],
    'objects': [],
    'provider': [],
    'ident': ['id'],
    'best_h_link': ['h_link'],
    'star_basic': ['coo','plx','dist_st','coo_gal','mag_i',
                                      'mag_j','mag_k','class','sptype'],
    'planet_basic': ['mass_pl'],
    'disk_basic': ['rad'],
    'mes_mass_pl': ['mass_pl'],
    'mes_teff_st': ['teff_st'],
    'mes_radius_st': ['radius_st'],
    'mes_mass_st': ['mass_st'],
    'mes_binary': ['binary'],
    'mes_sep_ang': ['sep_ang'],
    'h_link': ['h_link']    
}

empty_cat = {
    'sources' : Table(),
    'objects': Table(),
    'provider': Table(),
    'ident': Table(),
    'best_h_link': Table(),
    'star_basic': Table(),
    'planet_basic': Table(),
    'disk_basic': Table(),
    'mes_mass_pl': Table(),
    'mes_teff_st': Table(),
    'mes_radius_st': Table(),
    'mes_mass_st': Table(),
    'mes_binary': Table(),
    'mes_sep_ang': Table(),
    'h_link': Table()
}

empty_cat_wit_columns = {
    'sources' : Table(names=['ref','provider_name','source_id'],
                        dtype=[object,object,int]),
    'objects': Table(names=['type','ids','main_id'],
                            dtype=[object,object,object]),
    'provider': Table(names=['provider_name','provider_url',
                             'provider_bibcode','provider_access'],
                            dtype=[object,object,object,object]),
    'ident': Table(names=['object_idref','id','id_source_idref'],
                            dtype=[int,object,int]),
    'best_h_link': Table(names=['child_object_idref','parent_object_idref',
                             'h_link_source_idref','h_link_ref','membership'],
                            dtype=[int,int,int,object,int]),
    'star_basic': Table(names=['object_idref',
                               'coo_ra','coo_dec','coo_err_angle','coo_err_maj',
                             'coo_err_min','coo_qual','coo_source_idref',
                             'coo_ref',
                         'coo_gal_l','coo_gal_b','coo_gal_err_angle',
                             'coo_gal_err_maj','coo_gal_err_min','coo_gal_qual',
                             'coo_gal_source_idref','coo_gal_ref',
                        'mag_i_value','mag_i_err','mag_i_qual',
                             'mag_i_source_idref','mag_i_ref',
                        'mag_j_value','mag_j_err','mag_j_qual',
                             'mag_j_source_idref','mag_j_ref',
                        'mag_k_value','mag_k_err','mag_k_qual',
                             'mag_k_source_idref','mag_k_ref',
                              'plx_value','plx_err','plx_qual',
                             'plx_source_idref','plx_ref',
                              'dist_st_value','dist_st_err','dist_st_qual',
                             'dist_st_source_idref','dist_st_ref',
                              'sptype_string','sptype_err','sptype_qual',
                            'sptype_source_idref','sptype_ref',
                              'class_temp','class_temp_nr','class_lum',
                            'class_source_idref','class_ref',
                              'teff_st_value','teff_st_err','teff_st_qual',
                            'teff_st_source_idref','teff_st_ref',
                              'radius_st_value','radius_st_err','radius_st_qual',
                            'radius_st_source_idref','radius_st_ref',
                              'mass_st_value','mass_st_err','mass_st_qual',
                            'mass_st_source_idref','mass_st_ref',
                              'binary_flag','binary_qual','binary_source_idref',
                            'binary_ref',
                              'sep_ang_value','sep_ang_err','sep_ang_obs_date',
                            'sep_ang_qual','sep_ang_source_idref',
                            'sep_ang_ref'],
                            dtype=[int,
                                   float,float,float,float,float,object,int,object,
                            float,float,float,float,float,object,int,object,
                            float,float,object,int,object,
                            float,float,object,int,object,
                            float,float,object,int,object,
                                  float,float,object,int,object,
                                  float,float,object,int,object,
                                  object,float,object,int,object,
                                  object,object,object,int,object,
                                  float,float,object,int,object,
                                  float,float,object,int,object,
                                  float,float,object,int,object,
                                  object,object,int,object,
                                  float,float,int,object,int,object]),
    'planet_basic': Table(names=['object_idref',
                                 'mass_pl_value','mass_pl_err','mass_pl_rel',
                            'mass_pl_qual','mass_pl_source_idref',
                            'mass_pl_ref'],
                          dtype=[int,
                                 float,float,object,object,int,object]),
    'disk_basic': Table(names=['object_idref',
                               'rad_value','rad_err','rad_rel','rad_qual',
                            'rad_source_idref','rad_ref'],
                        dtype=[int,
                              float,float,object,object,int,object]
                       ),
    'mes_mass_pl': Table(names=['object_idref',
                                'mass_pl_value','mass_pl_err','mass_pl_rel',
                            'mass_pl_qual','mass_pl_source_idref',
                            'mass_pl_ref'],
                        dtype=[int,
                              float,float,object,object,int,object]),
    'mes_teff_st': Table(names=['object_idref',
                                'teff_st_value','teff_st_err','teff_st_qual',
                            'teff_st_source_idref','teff_st_ref'],
                        dtype=[int,
                               float,float,object,int,object]),
    'mes_radius_st': Table(names=['object_idref',
                                 'radius_st_value','radius_st_err','radius_st_qual',
                            'radius_st_source_idref','radius_st_ref'],
                        dtype=[int,
                               float,float,object,int,object]),
    'mes_mass_st': Table(names=['object_idref',
                               'mass_st_value','mass_st_err','mass_st_qual',
                            'mass_st_source_idref','mass_st_ref'],
                        dtype=[int,
                              float,float,object,int,object]),
    'mes_binary': Table(names=['object_idref',
                              'binary_flag','binary_qual','binary_source_idref',
                            'binary_ref'],
                        dtype=[int,
                              object,object,int,object]),
    'mes_sep_ang': Table(names=['object_idref',
                               'sep_ang_value','sep_ang_err','sep_ang_obs_date',
                            'sep_ang_qual','sep_ang_source_idref',
                            'sep_ang_ref'],
                        dtype=[int,
                              float,float,int,object,int,object]),
    'h_link': Table(names=['child_object_idref','parent_object_idref',
                             'h_link_source_idref','h_link_ref','membership'],
                        dtype=[int,int,int,object,int])
}

empty_provider_tables_dict = {
    'sim': empty_cat.copy(),
    'sdb': empty_cat.copy(),
    'exo': empty_cat.copy(),
    'life': empty_cat.copy(),
    'gaia': empty_cat.copy(),
    'wds': empty_cat.copy()
}                 
        
        
