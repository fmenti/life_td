"""
Helper functions for the creation and analysis of the LIFE Target database.
Author: Franziska Menti 22.12.2023
"""

import numpy as np #arrays
import astropy as ap #votables
import matplotlib.pyplot as plt

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
        	    ap.io.votable.from_table(cat), f'../data/{path}.xml')
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
        to_append=ap.io.votable.parse_single_table(f'../data/{path}.xml')
        cats.append(to_append.to_table())
    #go through all the tables in the cats list
    if stringtoobjects:
        for cat in cats:
            cat=stringtoobject(cat,3000)
    return cats

def object_contained(stars,cat,details=False):
    """
    This function checks which of the star identifiers in stars are not 
    present in cat. If details is True then the amount and individual objects
    which are in stars but not cat are printed. Returns only star identifiers 
    in stars that are also present in cat.
    :param stars: numpy string array containing star identifiers
    :param cat: numpy string array containing star identifiers
    :param details: Bool, default is False.
    :return stars: numpy string array containing star identifiers 
        which are also present in cat
    """
    lost=stars[np.where(np.invert(np.in1d(stars,cat)))]
    if len(lost)>0:
        print('This criterion was not met by:',len(lost))
        if details:
            print(lost)
    stars=stars[np.where(np.invert(np.in1d(stars,lost)))]
    return stars

def compare_catalogs(cat1,cat2,cat1_idname,cat2_idname,cat1_paranames,cat2_paranames):
    """
    This function
    :param cat1: astropy table to be compared to cat2
    :param cat2: astropy table to be compared to cat1
    :param cat1_idname: column name of cat 1 column containing object identifiers
    :param cat2_idname: column name of cat 2 column containing object identifiers
    :param cat1_paranames: column names of cat 1 columns containing measurements to be compared to cat 2
    :param cat2_paranames: column names of cat 2 columns containing measurements 
        to be compared to cat 1. Order of column names needs to be same as in cat1_parameters
    """
    print('lenght cat 1:',len(cat1))
    print('lenght cat 2:',len(cat2))
    common=cat1[np.where(np.in1d(cat1[cat1_idname],cat2[cat2_idname]))]
    print('number of common objects:',len(common))
    print('number of objects in cat 1 but not cat 2:',len(cat1)-len(common))
    print('number of objects in cat 2 but not cat 1:',len(cat2)-len(common))
    common2=cat2[np.where(np.in1d(cat2[cat2_idname],cat1[cat1_idname]))]
    difference=ap.table.join(common, common2,keys_left=cat1_idname,keys_right=cat2_idname)
        # still need to do some detail work on colnames. I can call all paranames of 1 with 1 in end
    for i in range(len(cat1_paranames)):
        print(f'Distribution of parameter {cat1_paranames[i]}')
        if cat1_paranames[i]!=cat2_paranames[i]:
            print(f' respectively {cat2_paranames[i]}')
            
        # removing objects with masked or nan values:
        print('tbd remove objects with masked or nan values')
        #histogram with three subplots (cat1,cat2,common)
        fig, ((ax1,ax2,ax3)) = plt.subplots(1,3,\
                       figsize = [16, 5],\
                       gridspec_kw={'hspace': 0, 'wspace': 0},sharey=True) 
        
        ax1.set_xlabel(cat1_paranames[i])
        ax1.set_ylabel('Number of objects')
        ax1.hist(cat1[cat1_paranames[i]],histtype='bar',log=True)
        #ax1.set_xscale('log')
    
        ax2.set_xlabel(cat2_paranames[i])
        ax2.hist(cat2[cat2_paranames[i]],histtype='bar',log=True)
                 
        ax3.set_xlabel(cat1_paranames[i])
        ax3.hist(common[cat1_paranames[i]],histtype='bar',log=True)
        plt.show()
        
        plt.figure
        plt.xlabel(cat1_paranames[i])
        plt.ylabel('Number of objects')
        plt.hist(cat1[cat1_paranames[i]],histtype='bar',log=True,color='C0', alpha=0.8, label='cat1')
        plt.hist(cat2[cat2_paranames[i]],histtype='bar',log=True,color='C1', alpha=0.8, label='cat2')
        plt.hist(common[cat1_paranames[i]],histtype='bar',log=True,color='C2', alpha=0.8, label='common')
        plt.legend()
        plt.show()
        print('cat1, cat2 and common objects with values of cat1')
        
        #mean
        print('cat 1 mean and std:',np.mean(cat1[cat1_paranames[i]]),np.std(cat1[cat1_paranames[i]]))
        #std (The standard deviation is the square root of the average of 
        #   the squared deviations from the mean, i.e., std = sqrt(mean(x)), where x = abs(a - a.mean())**2)
        print('cat 2 mean and std:',np.mean(cat2[cat2_paranames[i]]),np.std(cat2[cat2_paranames[i]]))
        print('common objects cat 1 values mean and std:',\
              np.mean(common[cat1_paranames[i]]),np.std(common[cat1_paranames[i]]))
        print('common objects cat 2 values mean and std:',\
              np.mean(common2[cat2_paranames[i]]),np.std(common2[cat2_paranames[i]]))
        print(f'Distribution of difference of parameter {cat1_paranames[i]}')
        if cat1_paranames[i]!=cat2_paranames[i]:
            print(f' respectively {cat2_paranames[i]}')
            difference[f'diff_{cat1_paranames[i]}']=difference[cat1_paranames[i]]-difference[cat2_paranames[i]]
        else:
            difference[f'diff_{cat1_paranames[i]}']=difference[cat1_paranames[i]+'_1']-difference[cat2_paranames[i]+'_2']
        print('differences of common objects values mean and std:',\
              np.mean(difference[f'diff_{cat1_paranames[i]}']),np.std(difference[f'diff_{cat1_paranames[i]}']))
        plt.figure
        plt.xlabel(f'diff_{cat1_paranames[i]}')
        plt.ylabel('Number of objects')
        plt.hist(difference[f'diff_{cat1_paranames[i]}'],histtype='bar',log=True)
        plt.show()
    return
