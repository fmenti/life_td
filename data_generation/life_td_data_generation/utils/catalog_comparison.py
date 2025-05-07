"""
Helper functions for the creation and analysis of the LIFE Target 
    Database.
"""

import numpy as np #arrays
import astropy as ap #votables
import matplotlib.pyplot as plt


def detail_criteria(database_tables, l):
    """
    This code scans for reasons, why the given objects would not be included in StarCat4

    """
    print(
        'tbd: file starcat4 analysis where this code is used to give specific reason why a given object was not included into the cat4')
    # to do:
    # include get main id search
    # this code could use some refactoring. and also printing the corresponding objects why it was not included.
    # write tests
    cat_i = database_tables['ident']
    cat_o = database_tables['objects']
    cat_h = database_tables['best_h_link']  # -> use best_h_link
    children = []
    not_found = []
    number_of_stellar_siblings = []
    multiple_parents = []  # multiple components in wds
    single_child = []
    binary = []
    higher_order_multiple = []
    star_without_parent = []
    system_without_child = []
    for name in l:
        if len(cat_i['main_id'][np.where(cat_i['id'] == name)]) > 0:
            main_id = cat_i['main_id'][np.where(cat_i['id'] == name)][0]
            if cat_o['type'][np.where(cat_o['main_id'] == main_id)][0] == 'sy':
                if len(cat_h[np.where(cat_h['parent_main_id'] == main_id)]) > 0:
                    # print('system object but with found child object',main_id)
                    # print(cat_h[np.where(cat_h['parent_main_id']==main_id)])
                    for child in cat_h['child_main_id'][np.where(cat_h['parent_main_id'] == main_id)]:
                        children.append(child)
                else:
                    system_without_child.append(name)
            elif cat_o['type'][np.where(cat_o['main_id'] == main_id)][0] == 'st':
                # if it has parents:
                if len(cat_h[np.where(cat_h['child_main_id'] == main_id)]) > 0:
                    # print('system object but with found child object',main_id)
                    # print(cat_h[np.where(cat_h['parent_main_id']==main_id)])
                    if len(cat_h[np.where(cat_h['child_main_id'] == main_id)]) > 1:
                        multiple_parents.append(name)
                    elif len(cat_h[np.where(cat_h['child_main_id'] == main_id)]) == 1:
                        st_sib = 0
                        nestled = False
                        parent = cat_h['parent_main_id'][np.where(cat_h['child_main_id'] == main_id)]
                        # wrong code here, parent can't be type st
                        for sibling in cat_h['child_main_id'][np.where(cat_h['parent_main_id'] == parent)]:
                            #print(cat_h[np.where(cat_h['child_main_id'] == sibling)])
                            if cat_o['type'][np.where(cat_o['main_id'] == sibling)][0] == 'sy':
                                nestled = True
                            elif cat_o['type'][np.where(cat_o['main_id'] == sibling)][0] == 'st':
                                st_sib += 1
                            number_of_stellar_siblings.append(st_sib)
                        if st_sib == 1 and nestled == False:
                            single_child.append(name)
                        if st_sib > 2 or nestled == True:
                            higher_order_multiple.append(name)
                        if st_sib == 2 and nestled == False:
                            binary.append(name)
                else:
                    star_without_parent.append(name)
            else:
                print('no sys or st \n', name, cat_o['main_id', 'type'][np.where(cat_o['main_id'] == main_id)][0])
        else:
            # print('object not found',name)
            not_found.append(name)
    print('Of the', len(l), 'objects given:')
    print('Some are not in cat 4 because they are either:')
    print('system_without_child', system_without_child, len(system_without_child))
    print('star_without_parent', star_without_parent, len(star_without_parent))
    print('not found', not_found, '\n')
    print('Some were not included because they are systems and instead one should look at their child objects:')
    print('children', children, '\n')
    print('some are non trivial binaries:')
    print('multiple parents', multiple_parents, len(multiple_parents))
    print('higher_order_multiple', higher_order_multiple, len(higher_order_multiple), '\n')
    # print('# stellar siblings',number_of_stellar_siblings)
    if len(single_child) > 0:
        print('single_child', single_child, len(single_child))
    print('And the reminder have conpanions that don t fit the spectral type requirements')
    print('trivial binary', binary, len(binary))
    return


def object_contained(stars,cat,details=False):
    """
    Checks which of the entries in stars is also in cat.
    
    This function checks which of the star identifiers in stars are not 
    present in cat. If details is True then the amount and individual objects
    which are in stars but not cat are printed. Returns only star identifiers 
    in stars that are also present in cat.
    
    :param stars: numpy string array containing star identifiers
    :type stars:
    :param cat: numpy string array containing star identifiers
    :type cat:
    :param bool details: Defaults to False.
    :returns: numpy string array containing star identifiers 
        which are also present in cat
    :rtype: 
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
    This function compares two catalogs.
    
    :param cat1:  to be compared to cat2
    :type cat1: astropy.table.table.Table
    :param cat2:  to be compared to cat1
    :type cat2: astropy.table.table.Table
    :param str cat1_idname: column name of cat 1 column containing object 
        identifiers
    :param str cat2_idname: column name of cat 2 column containing object 
        identifiers
    :param cat1_paranames: column names of cat 1 columns containing 
        measurements to be compared to cat 2
    :type cat1_paranames: list(str)
    :param cat2_paranames: column names of cat 2 columns containing measurements 
        to be compared to cat 1. Order of column names needs to be same as in cat1_parameters
    :type cat2_paranames: list(str)
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

def create_common(cat1,cat2):
    """
    Takes two catalogs, removes objects that are not in the other one.
    
    :param cat1: First catalog.
    :type cat1: astropy.table.table.Table
    :param cat2: Second catalog.
    :type cat2: astropy.table.table.Table
    :returns: Touple of cat1 containing only objects that are also in 
        cat2 and cat2 containing only objects that are also in cat1.
    :rtype: list(astropy.table.table.Table)
    """
    
    common_stars=np.intersect1d(list(cat1['main_id']),list(cat2['main_id']))
    common_cat1=cat1[np.where(np.in1d(cat1['main_id'],common_stars))]
    common_cat2=cat2[np.where(np.in1d(cat2['main_id'],common_stars))]   
    common_cat2=stringtoobject(common_cat2,number=1000)
    common_cat1=stringtoobject(common_cat1,number=1000)
    return common_cat1, common_cat2
