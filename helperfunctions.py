"""
Helper functions for the creation and analysis of the LIFE Target database.
Author: Franziska Menti 22.12.2023
"""

import numpy as np #arrays
import astropy as ap #votables

#-------------------global helper functions------------------------------------
def save(cats,names):
    """
    This functions saves the tables given as python list in the cats parameter.
    The saving location is 'data/{name}.xml' where path is given in the paths
    parameter.
    :param cats: Python list of astropy table to be saved.
    :param names: Python list of strings containing names for saving location
        of tables in cats.
    """
    #go through all the elements in both lists
    for cat,path in zip(cats,names):
        #for each column header
        for i in cat.colnames:
            #if the type of the column is object (=adaptable length string)
            if cat[i].dtype == object:
                #transform the type into string
                cat[i] = cat[i].astype(str)
        #save the table
        ap.io.votable.writeto(
        	    ap.io.votable.from_table(cat), f'data/{path}.xml')
    return

def stringtoobject(cat,number=100):
    """
    This function changes from string to object format.
    The later has the advantace of allowing strings of varying length.
    Without it truncation of entries is a risk.
    :param cat: Astropy table.
    :param number: Length of longest string type element in the table.
        Default is 100.
    :return cat: Astropy table with all string columns transformed into
        object type ones.
    """
    #defining string types as calling them string does not work and instead
    #the type name <U3 is needed for a string of length 3
    stringtypes=[np.dtype(f'<U{j}') for j in range(1,number)]
    #for each column header
    for i in cat.colnames:
        #if the type of the column is string
        if cat[i].dtype in stringtypes:
            #transform the type into object
            cat[i] = cat[i].astype(object)
    return cat

def load(paths,stringtoobjects=True):
    """
    This function loads the tables saved in XML format at saving locations
    specified in paths. If stringtoobject is True the function 
    stringtoobjects is invoked.
    :param paths: Python list of saving locations.
    :return cats: Python list of loaded astropy tables.
    """
    #initialize return parameter as list
    cats=[]
    #go through all the elements in the paths list
    for path in paths:
        #read the saved data into the cats lists as astropy votable element
        to_append=ap.io.votable.parse_single_table(f'data/{path}.xml')
        cats.append(to_append.to_table())
    #go through all the tables in the cats list
    if stringtoobjects:
        for cat in cats:
            cat=stringtoobject(cat,3000)
    return cats