"""
Helper functions for the creation and analysis of the LIFE Target 
    Database.
"""

from numpy import dtype
from astropy.io import votable

class Path:
    def __init__(self):
        """Constructor method"""
        self.data='../../data/'
        self.additional_data='../../additional_data/'
        self.plot='../../plots/'
    

def objecttostring(cat):
    for i in cat.colnames:
        if cat[i].dtype == object: #object =adaptable length string
            #transform the type into string
            cat[i] = cat[i].astype(str)
    return cat
    
def save(cats,names,location=Path().additional_data):
    """
    This functions saves the tables given as list in the cats parameter.
    
    :param cats: Python list of astropy table to be saved.
    :type cats: list(astropy.table.table.Table)
    :param names: Contains names for saving location
        of tables in cats.
    :type names: list(str)
    :param str location: Defaults to ../../data/additional_data/
    """
    
    #go through all the elements in both lists
    for cat,path in zip(cats,names):
        temp=cat.copy()
        #for each column header
        temp = objecttostring(temp)
        #save the table
        votable.writeto(
        	    votable.from_table(temp), f'{location}{path}.xml')
    return

def stringtoobject(cat,number=100):
    """
    This function changes string type columns to object type.
    
    The later has the advantace of allowing strings of varying length.
    Without it truncation of entries is a risk.
    
    :param cat: Table with at least two columns
    :type cat: astropy.table.table.Table
    :param int number: Length of longest string type element in the table.
        Default is 100.
    :returns: Table with all string type columns transformed into object type ones.
    :rtype: astropy.table.table.Table
    """
    
    # defining string types as calling them string does not work and 
    # instead the type name <U3 is needed for a string of length 3
    stringtypes=[dtype(f'<U{j}') for j in range(1,number)]
    #for each column header
    for i in cat.colnames:
        #if the type of the column is string
        if cat[i].dtype in stringtypes:
            #transform the type into object
            cat[i] = cat[i].astype(object)
    return cat

def string_to_object_whole_dict(dictionary,number=100):
    for table_name in list(dictionary.keys()):
        dictionary[table_name] = stringtoobject(dictionary[table_name],number)
    return dictionary

def load(paths,stringtoobjects=True,location=Path().additional_data):
    """    
    This function loads xml tables. 
    
    :param paths: Filenames.
    :type paths: list(str)
    :param stringtoobjects: Wheter stringtoobject function should be 
        called.
    :type stringtoobjects: bool
    :param location: Folder to save the file in, default is ../../data/additional_data/
    :type location: str
    :returns: Loaded tables.
    :rtype: list(astropy.table.table.Table)
    """
    
    #initialize return parameter as list
    cats=[]
    #go through all the elements in the paths list
    for path in paths:
        #read the saved data into the cats lists as astropy votable element
        to_append=votable.parse_single_table(f'{location}{path}.xml')
        cats.append(to_append.to_table())
    #go through all the tables in the cats list
    if stringtoobjects:
        for cat in cats:
            cat=stringtoobject(cat,3000)
    return cats