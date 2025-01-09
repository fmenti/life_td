from life_td import *

#def integration_test_create_life_td():
#    provider_tables_dict, database_tables = create_life_td(5.)
    
    #curves... database_tables['mes_mass_st']
    #assert that difference data to curve below value
    
#def integration_test_load...

def integration_test_fake_provider_data():
    #sdb=...
    fake_sdb = {
    'sources' : Table(data=['priv. comm.','Grant Kennedy Disks'],names=['ref','provider_name'],
                        dtype=[object,object]),
    'objects': Table(data=names=['type','ids','main_id'],
                            dtype=[object,object,object]),
    'provider': Table(names=['provider_name','provider_url',
                             'provider_bibcode','provider_access'],
                            dtype=[object,object,object,object]),
    'ident': Table(names=['object_idref','id','id_source_idref'],
                            dtype=[int,object,int]),
    'best_h_link': Table(),
    'star_basic': Table(),
    'planet_basic': Table(),
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
    #sim=...
    #provider_tables_dict, database_tables = create_life_td(5.)
    #assert
    
    