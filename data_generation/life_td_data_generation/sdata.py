"""
Structured data class
"""

import astropy as ap

import utils as hf

#now create a class that returns an astropy table
class structured_data:
    """
    A class that structures data for ingestion into life_td DaCHS.
    
    :param str name: Name of class object
    :param table_names: Names of the tables
    :type table_names: list(str)
    :param list_of_tables: Tables
    :type list_of_tables: list(astropy.table.table.Table)
    """
    
    def __init__(self,name):
        """Constructor method"""
        
        self.name = name
        # tbd: instead of table_names and list_of_tables lists have them 
        # both in one dictionary
        # self.dict= dict(table_name = , table = )
        # current issue, how to initialize an emty dictionary
        self.table_names = []
        self.list_of_tables = []
    
    def add_table(self,table_name):
        """
        Adds a table and its table_name to
        
        :param str table_name: Name of the to be added table
        """
        
        self.table_name = table_name
        self.table_names.append(table_name)
        self.list_of_tables.append(ap.table.Table())
        
    def table(self,table_name):
        """
        Takes table_name, returns corresponding table
        
        :param str table_name: Name of table
        :returns: table corresponding to table_name
        """
        
        self.table_name = table_name
        return self.list_of_tables[self.table_names.index(self.table_name)]
    
    def add_parameter(self,table_name,parameter,para_type):
        """
        Adds a parameter as a column to the specified table.
        
        Careful, only adding columns, not updating of them.
        
        :param str table_name: Name of table
        :param parameter: Name of parameter to be added.
        :type parameter: str
        :param para_type: Data types for parameter.
        :type para_type: dtype
        """
        
        col=ap.table.Column(name=parameter,dtype=para_type)
        self.list_of_tables[self.table_names.index(table_name)].add_column(col)
    
    def add_parameters(self,table_name,parameters,para_types):
        """
        Adds the specified parameters to the specified table.
        
        :param str table_name: Name of table
        :param parameters: Names of parameters to be added.
        :type parameters: list(str)
        :param para_types: Data types for parameter.
        :type para_types: list(dtype)
        """
        
        for para,ptype in zip(parameters,para_types):
            self.add_parameter(table_name,para,ptype)
    
    def save(self,path='../../data/'):
        """
        Saves the individual tables.
        
        :param str path: Directory for saving location, defaults to
            ../../data/
        """
        
        hf.save(self.list_of_tables,self.table_names,
                location=path+self.name+'_')
        
    def load(self,path='../../data/'):
        """
        Loads the individual tables.
        
        :param str path: Directory for loading, defaults to ../../data/
        """
        
        print('tbd')
        
        
class provider(structured_data):
    """
    Template to later store in data from different providers
    
    :param str name: Name of class object
    """
    
    def __init__(self,name):
        """Constructor method"""
        
        super().__init__(name)
        
        # table for provenance information of data
        self.add_table('sources')
        self.add_parameters('sources',
                        ['ref','provider_name','source_id'],
                        [object,object,int])
                        
        # table for object identification
        self.add_table('objects')
        # database internal object identifier
        self.add_parameter('objects','object_id',int)
        self.add_parameters('objects',
                            ['type','ids','main_id'],
                            [object,object,object])

        # table for collection of metadata on provider
        self.add_table('provider')
        self.add_parameters('provider',
                            ['provider_name','provider_url',
                             'provider_bibcode','provider_access'],
                            [object,object,object,object])

        # table for collection of individual identifiers
        self.add_table('ident')
        self.add_parameters('ident',
                            ['object_idref','id','id_source_idref'],
                            [int,object,int])

        # table for collection of hierarchical link meaning relation 
        # between objects. Only best relation for same object pair permitted.
        self.add_table('best_h_link')
        self.add_parameters('best_h_link',
                            ['child_object_idref','parent_object_idref',
                             'h_link_source_idref','h_link_ref','membership'],
                            [int,int,int,object,int])

        # table for collection of basic stellar data
        self.add_table('star_basic')
        # database internal object identifier
        self.add_parameter('star_basic','object_idref',int)
        
        # coordinate data
        self.add_parameters('star_basic',
                            ['coo_ra','coo_dec','coo_err_angle','coo_err_maj',
                             'coo_err_min','coo_qual','coo_source_idref',
                             'coo_ref'],
                            [float,float,float,float,float,object,int,object])
        # galactical coordinate data
        self.add_parameters('star_basic',
                            ['coo_gal_l','coo_gal_b','coo_gal_err_angle',
                             'coo_gal_err_maj','coo_gal_err_min','coo_gal_qual',
                             'coo_gal_source_idref','coo_gal_ref'],
                            [float,float,float,float,float,object,int,object])

        # magnitude data
        self.add_parameters('star_basic',
                            ['mag_i_value','mag_i_err','mag_i_qual',
                             'mag_i_source_idref','mag_i_ref'],
                            [float,float,object,int,object])
        self.add_parameters('star_basic',
                            ['mag_j_value','mag_j_err','mag_j_qual',
                             'mag_j_source_idref','mag_j_ref'],
                            [float,float,object,int,object])
        self.add_parameters('star_basic',
                            ['mag_k_value','mag_k_err','mag_k_qual',
                             'mag_k_source_idref','mag_k_ref'],
                            [float,float,object,int,object])

        # parallax data
        self.add_parameters('star_basic',
                            ['plx_value','plx_err','plx_qual',
                             'plx_source_idref','plx_ref'],
                            [float,float,object,int,object])
        # distance data
        self.add_parameters('star_basic',
                            ['dist_st_value','dist_st_err','dist_st_qual',
                             'dist_st_source_idref','dist_st_ref'],
                            [float,float,object,int,object])
        
        # spectral type data, as concatenaded string
        self.add_parameters('star_basic',
                            ['sptype_string','sptype_err','sptype_qual',
                            'sptype_source_idref','sptype_ref'],
                            [object,float,object,int,object])
        # spectral type data, as temperature and luminosity class
        self.add_parameters('star_basic',
                            ['class_temp','class_temp_nr','class_lum',
                            'class_source_idref','class_ref'],
                            [object,object,object,int,object])
                
        # effective temperature data
        self.add_parameters('star_basic',
                            ['teff_st_value','teff_st_err','teff_st_qual',
                            'teff_st_source_idref','teff_st_ref'],
                            [float,float,object,int,object])
                            
        # radius data
        self.add_parameters('star_basic',
                            ['radius_st_value','radius_st_err','radius_st_qual',
                            'radius_st_source_idref','radius_st_ref'],
                            [float,float,object,int,object])
        
        # mass data
        self.add_parameters('star_basic',
                            ['mass_st_value','mass_st_err','mass_st_qual',
                            'mass_st_source_idref','mass_st_ref'],
                            [float,float,object,int,object])                  
        
        # binarity data
        self.add_parameters('star_basic',
                            ['binary_flag','binary_qual','binary_source_idref',
                            'binary_ref'],
                            [object,object,int,object])
        # binary angular separation data                    
        self.add_parameters('star_basic',
                            ['sep_ang_value','sep_ang_err','sep_ang_obs_date',
                            'sep_ang_qual','sep_ang_source_idref',
                            'sep_ang_ref'],
                            [float,float,int,object,int,object])                   
        
        # table for collection of basic planetary data
        self.add_table('planet_basic')
        
        # database internal object identifier
        self.add_parameter('planet_basic','object_idref',int)
                
        # mass data
        self.add_parameters('planet_basic',
                            ['mass_pl_value','mass_pl_err','mass_pl_rel',
                            'mass_pl_qual','mass_pl_source_idref',
                            'mass_pl_ref'],
                            [float,float,object,object,int,object])
                                    
        # table for collection of basic disk data
        self.add_table('disk_basic')
        # database internal object identifier
        self.add_parameter('disk_basic','object_idref',int)
        
        # black body radius data
        self.add_parameters('disk_basic',
                            ['rad_value','rad_err','rad_rel','rad_qual',
                            'rad_source_idref','rad_ref'],
                            [float,float,object,object,int,object])
        
        # tables for multiple measurements per object
        # table for collecting planetary mass data
        self.add_table('mes_mass_pl')
        # database internal object identifier
        self.add_parameter('mes_mass_pl','object_idref',int)
        # mass data
        self.add_parameters('mes_mass_pl',
                            ['mass_pl_value','mass_pl_err','mass_pl_rel',
                            'mass_pl_qual','mass_pl_source_idref',
                            'mass_pl_ref'],
                            [float,float,object,object,int,object])
        
        # table for collecting stellar effective temperature data
        self.add_table('mes_teff_st')
        # database internal object identifier
        self.add_parameter('mes_teff_st','object_idref',int)
        # effective temperature data
        self.add_parameters('mes_teff_st',
                            ['teff_st_value','teff_st_err','teff_st_qual',
                            'teff_st_source_idref','teff_st_ref'],
                            [float,float,object,int,object])
                            
        # table for collecting stellar radius data
        self.add_table('mes_radius_st')
        # database internal object identifier
        self.add_parameter('mes_radius_st','object_idref',int)
        # radius data
        self.add_parameters('mes_radius_st',
                            ['radius_st_value','radius_st_err','radius_st_qual',
                            'radius_st_source_idref','radius_st_ref'],
                            [float,float,object,int,object])
                            
        # table for collecting stellar mass data
        self.add_table('mes_mass_st')
        # database internal object identifier
        self.add_parameter('mes_mass_st','object_idref',int)
        # mass data
        self.add_parameters('mes_mass_st',
                            ['mass_st_value','mass_st_err','mass_st_qual',
                            'mass_st_source_idref','mass_st_ref'],
                            [float,float,object,int,object]) 
                            
        # table for collecting binarity data
        self.add_table('mes_binary')
        # database internal object identifier
        self.add_parameter('mes_binary','object_idref',int)
        # binarity data
        self.add_parameters('mes_binary',
                            ['binary_flag','binary_qual','binary_source_idref',
                            'binary_ref'],
                            [object,object,int,object])
                            
        # table for collecting binary angular separation data
        self.add_table('mes_sep_ang')
        # database internal object identifier
        self.add_parameter('mes_sep_ang','object_idref',int)
        # binary angular separation data                    
        self.add_parameters('mes_sep_ang',
                            ['sep_ang_value','sep_ang_err','sep_ang_obs_date',
                            'sep_ang_qual','sep_ang_source_idref',
                            'sep_ang_ref'],
                            [float,float,int,object,int,object]) 
                              
        # table for collection of hierarchical link meaning relation 
        # between objects data. multiple relations for same object pair ok.
        self.add_table('h_link')
        # binary angular separation data                    
        self.add_parameters('h_link',
                            ['child_object_idref','parent_object_idref',
                             'h_link_source_idref','h_link_ref','membership'],
                            [int,int,int,object,int])                   
        
        
