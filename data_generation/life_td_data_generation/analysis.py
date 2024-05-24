import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables
from datetime import datetime
import matplotlib.pyplot as plt

#self created modules
import utils as hf

###############################################################################
#-------------------------Sanity tests--------------------------------------------
###############################################################################
    
def show(provider,
         table='objects',
         columns=[],
         wherecol='main_id',
         whereobj='* zet02 Ret'):
    """
    This function prints the specified columns of a table.
    
    :param provider: database_tables, sim, wds, gaia, gk, life
    :param table:
    :param columns: if empty prints all
    :param wehercol:
    :param whereobj:
    """
    
    table_names=['sources','objects','provider','ident','h_link','star_basic',
              'planet_basic','disk_basic','mes_mass_pl',
              'mes_teff_st','mes_radius_st','mes_mass_st','mes_binary','mes_sep_ang','best_h_link']
    
    cat=provider[table_names.index(table)]

    if columns==[]:
        print(cat[np.where(cat[wherecol]==whereobj)])
    else:
        print(cat[columns][np.where(cat[wherecol]==whereobj)])
    return
    
    
def sanitytest(cat,colname):
    """
    Prints histogram of parameter for quick sanity test.
    
    :param cat: Table containing column colname.
    :type cat: astropy.table.table.Table
    :param str colname: Column name
    """
    
    arr=cat[colname]
    if len(arr)==0 or len(arr[np.where(arr!=1e20)])==0:
        print(colname,'is empty')
        return
    elif max(arr)==1e20:
        arr=arr[np.where(arr!=1e20)]
    plt.figure
    plt.xlabel(f'{colname}')
    plt.ylabel('Number of objects')
    plt.hist(arr,histtype='bar',log=True)
    plt.show()
    return



def disthist(arr,names,path,max_dist,xaxis):
    '''
    Makes a histogram of the stellar distance distribution.
    
    :param arr: ndarray containing the distance for 
        each star for multiple sets of stars
    :param names: list containing the labels for the sets
    :param path: location where the plot is saved
    '''
    
    n_bins=np.arange(0,max_dist+1,2.5)
    
    plt.figure
    plt.xlabel(xaxis)
    plt.ylabel('Number of objects')
    plt.title('Distance distribution')
    plt.hist(arr,n_bins,histtype='bar',log=True)
    plt.legend(names)
    #plt.savefig(path)
    plt.show()
    return


def spechist(spectypes,mute=False):
    '''
    Makes a histogram of the spectral type distribution.
    
    :param spectypes: array containing a spectral type for each star
    :type spectypes: 
    :param bool mute: if true prints the number of total stars, total 
        spectral types and spec and specdist
    :returns: array containig the spectral type abreviations 
        O B A F G K M. array containing the number of stars for
        each spectral type in spec.
    :rtype: touple(np.array)
    '''
    
    spec=np.array(['O','B','A','F','G','K','M'])
    s=len(spec)
    specdist =np.zeros(s)
    for i in range(len(spectypes)):
        if spectypes[i] not in ['','nan']:
            for j in range(0,s):
                if spec[j] == spectypes[i][0]: 
                    specdist[j] += 1
    if mute==False:
        print('total stars:',np.shape(spectypes)[0],\
              'total spectral types:',int(np.sum(specdist)+1))
        print('spectral type distribution:\n',spec,specdist)
    specdist=specdist.astype(int)
    return spec, specdist

def final_plot(stars,labels,distance_cut_in_pc,path='../plots/final_plot.png', \
                color=['tab:blue','tab:orange','tab:green']):
    """
    Plots spectral distribution in two subplots.
    
    Makes plot with two subfigures, first a histogram of spectral 
    type and then one of spectral type for each distance sample of 
    0-5, 5-10, 10-15 and 15-20pc.
    
    :param stars: list of astropy.table.table.Table of shape (,2) 
        with columns first spectral type and then distance in pc.
    :param labels: list containing the labels for the plot
    :param path: location to save the plot
    """
    
    steps=int(distance_cut_in_pc/5)

    n_legend=len(stars)
    spec=np.array(['O','B','A','F','G','K','M'])
    n_temp_class=len(spec)
    specdist=np.zeros((n_legend,n_temp_class))
    
    #prepares subfigure display
    fig, ((ax1,ax2)) = plt.subplots(1,2,sharey='row',\
                       figsize = [16, 5],\
                       gridspec_kw={'hspace': 0, 'wspace': 0})   
    
    x=np.arange(n_temp_class)
    #['r','b','yellow']
    width = 0.15
    
    #first subfigure
    for i in range(n_legend):
        specdist[i]=spechist(stars[i][0][:],mute=True)[1]
        #print(specdist[i])
        ax1.bar(x[2:]-n_legend*width/2+i*width,specdist[i][2:],
                width, align='center',label=labels[i],color=color[i])#,color=color[i])

    ax1.set_yscale('log') 
    #problem beim zweiten plot: min() arg is an empty sequence
    plt.sca(ax1)
    plt.xticks(x[2:], spec[2:])
    ax1.set_title(f"Spectral type distribution")
    ax1.set_ylabel('Number of stars')
    ax1.set_xlabel('Spectral types')
    ax1.legend(loc='upper left')

    #second subfigure 
    index = np.arange(steps) #wie viele bars pro label es haben wird
    #n = np.arange(7)#wieviele labels
    sub_specdist=np.zeros((n_legend,steps,n_temp_class))
           
    for i in range(n_legend):
        for j in range(steps):
            #here have an error because lifestarcat only goes up to 20pc
            if j*5.> max(stars[i][1][:]):
                sub_specdist[i][j]=[0.]*7
            else: 
                sub_specdist[i][j]=spechist(
                        stars[i][np.where((stars[i][1][:] >j*5.)*\
                        (stars[i][1][:] < (j+1)*5.))][0][:],mute=True)[1]
        for l in np.arange(n_temp_class)[2:]:
            #print(sub_specdist[i][:,l],
            #5*index-((n+(n+2)*(s-2))*width)+((n+2)*l+i)*width) 
            #problem if in one bin all 0
            #make if clause for it
            ax2.bar((steps+1)*index-((n_legend+(n_legend+2)*(n_temp_class-2))*width)+((n_legend+2)*l+i)*width,\
                    tuple(sub_specdist[i][:,l]),width,color=color[i])
         
    ax2.set_xlabel('Distance [pc]')
    plt.sca(ax2)
    if distance_cut_in_pc==20.:
        xticks_name=['0-5','5-10','10-15','15-20'][:steps]
    elif distance_cut_in_pc==25.:
        xticks_name=['0-5','5-10','10-15','15-20','20-25'][:steps]
    elif distance_cut_in_pc==30.:
        xticks_name=['0-5','5-10','10-15','15-20','20-25','25-30'][:steps]
    plt.xticks((steps+1)*index, (xticks_name))
    ax2.set_title(f"Spectral type and distance distribution")  
    plt.savefig(path, dpi=300)
    return

def spechistplot(stars,name,path=''):
    """
    Makes a histogram of the spectral distribution of the stars sample.
    
    :param stars: nd array containing at least one astropy table 
        with column spectral type
    :param name: list containing the labels for the plot
    :param path: location to save the plot
    """
    
    n=len(stars)
    spec=np.array(['O','B','A','F','G','K','M'])
    s=len(spec)
    specdist=np.zeros((n,s))
    
    fig, ax = plt.subplots()
    x=np.arange(s)
    width = 0.15

    for i in range(n):
        specdist[i]=spechist(stars[i],mute=True)[1]
        ax.bar(x[2:]-n*width/2+i*width,specdist[i][2:],
               width,label=name[i])
    
    ax.set_yscale('log')
    ax.set_ylabel('Number of objects')
    ax.set_title('Spectral type distribution')
    ax.set_xticks(x[2:])
    ax.set_xticklabels(spec[2:])
    ax.legend()
    fig.tight_layout()
    plt.savefig('../plots/'+path, dpi=300)
    plt.show()
    return

def objecthistplot(cat,name,path=''):
    """
    Makes a histogram of the object type distribution of the cat sample.
    
    :param cat: nd array containing at least one astropy table 
        with column object type
    :param name: list containing the labels for the plot
    :param path: location to save the plot
    """
    
    spec=np.array(['System','Star','Exoplanet','Disk'])
    
    plt.figure
    plt.title('Object type distribution')
    plt.xlabel('Number of objects')
    plt.hist(cat,histtype='bar',log=True,orientation='horizontal')
    plt.yticks(np.arange(4),spec)
    plt.savefig('../plots/'+path, dpi=300)
    plt.show()
    return

def sanity_tests(table_names,database_tables, distance_cut_in_pc):
    """
    TBD
    """
    
    data=database_tables
    #data=exo
    print('looking at table data and metadata \n')
    for i in range(len(table_names)):
        print(table_names[i],i,'\n')
        print(data[i].info)
        print(data[i][0:5],'\n','\n')

    print('looking at plots')

    ###problematic if I change order of database_tables tables. need to make it independent of that using table_names list
    cats=[database_tables[table_names.index('star_basic')],
      database_tables[table_names.index('disk_basic')],
      database_tables[table_names.index('mes_mass_pl')]]
    colnames=[['coo_ra','coo_dec','coo_err_angle',
           'coo_err_maj','coo_err_min','coo_qual',
           'mag_i_value','mag_j_value',
           'mag_k_value',
           'plx_value','plx_err','plx_qual'],
              ['rad_value','rad_err'],
               ['mass_pl_value','mass_pl_err']]
    for i in [0,1,2]:
        cat=cats[i]
        for colname in colnames[i]:
            sanitytest(cat,colname)

    disthist(database_tables[5]['dist_st_value'],['star_basic'],'test',distance_cut_in_pc,'dist_st_value')

    ltc3=ap.io.ascii.read("../data/additional_data/LIFE-StarCat3.csv")
    ltc3=hf.stringtoobject(ltc3,3000)
    print(ltc3['distance'])
    ltc3['class_temp']=ap.table.MaskedColumn(dtype=object,length=len(ltc3))
    for i in range(len(ltc3)):
        #sorting out entries like '', DA2.9, T1V
        if len(ltc3['sim_sptype'][i])>0 and ltc3['sim_sptype'][i][0] in ['O','B','A','F','G','K','M']:
            ltc3['class_temp'][i]=ltc3['sim_sptype'][i][0]
        else:
            ltc3['class_temp'][i]='?'
    table=database_tables[5]['class_temp','dist_st_value','binary_flag'][np.where(database_tables[5]['class_temp']!='?')]
    table=table['class_temp','dist_st_value'][np.where(table['binary_flag']=='False')]
    final_plot([table,
            ltc3['class_temp','distance'][np.where(ltc3['class_temp']!='?')]],['LTC4_singles','LTC3'],distance_cut_in_pc)



    table=database_tables[5]['class_temp','dist_st_value','binary_flag'][np.where(database_tables[5]['class_temp']!='?')]
    singles=table['class_temp'][np.where(table['binary_flag']=='False')]
    multiples=table['class_temp'][np.where(table['binary_flag']=='True')]
    total=table['class_temp']
    spechistplot([singles, multiples,total],['singles','multiples','total'])


    table=database_tables[1]['type']
    objecthistplot([table],['objects'],'objecthist')
    return
