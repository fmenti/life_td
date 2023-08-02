import matplotlib.pyplot as plt
import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables

#code snippet not running outside of environment of data_acquisition.py
###############################################################################
#-------------------------Sanity tests--------------------------------------------
###############################################################################

data=database_tables
#data=exo
print('looking at table data and metadata \n')
for i in range(len(table_names)):
    print(table_names[i],i,'\n')
    print(data[i].info)
    print(data[i][0:5],'\n','\n')
    
print('looking at plots')

def sanitytest(cat,colname):
    arr=cat[colname]
    plt.figure
    plt.xlabel(f'{colname}')
    plt.ylabel('Number of objects')
    plt.hist(arr,histtype='bar',log=True)
    plt.show()
    return

cats=[database_tables[5],database_tables[7],database_tables[8],database_tables[9]]
colnames=[['coo_ra','coo_dec','coo_err_angle','coo_err_maj','coo_err_min',
                'plx_value','plx_err'],
              ['rad_value','rad_err'],
                ['dist_st_value','dist_st_err'],
               ['mass_pl_value','mass_pl_err']]
for i in [0,1,2,3]:
    cat=cats[i]
    for colname in colnames[i]:
        sanitytest(cat,colname)

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

disthist(database_tables[5]['dist_st_value'],['star_basic'],'test',distance_cut_in_pc,'dist_st_value')

def spechist(spectypes,mute=False):
    '''
    Makes a histogram of the spectral type distribution.
    :param spectypes: array containing a spectral type for each star
    :param mute: if true prints the number of total stars, total 
        spectral types and spec and specdist
    :return spec: array containig the spectral type abreviations 
        O B A F G K M
    :return specdist: array containing the number of stars for
        each spectral type in spec
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
    #plt.savefig(path)
    plt.show()
    return

def final_plot(stars,labels,distance_cut_in_pc,path='results/final_plot.png',color=['tab:blue','tab:orange','tab:green']):
    """
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
    xticks_name=['0-5','5-10','10-15','15-20','20-25'][:steps]
    plt.xticks((steps+1)*index, (xticks_name))
    ax2.set_title(f"Spectral type and distance distribution")  
    #plt.savefig(path)
    return
ltc3=ap.io.ascii.read("data/LIFE-StarCat3.csv")
ltc3=stringtoobject(ltc3,3000)
print(ltc3['distance'])
ltc3['class_temp']=ap.table.MaskedColumn(dtype=object,length=len(ltc3))
for i in range(len(ltc3)):
    #sorting out entries like '', DA2.9, T1V
    if len(ltc3['sim_sptype'][i])>0 and ltc3['sim_sptype'][i][0] in ['O','B','A','F','G','K','M']:
        ltc3['class_temp'][i]=ltc3['sim_sptype'][i][0]
    else:
        ltc3['class_temp'][i]='?'

final_plot([database_tables[5]['class_temp','dist_st_value'][np.where(database_tables[5]['class_temp']!='?')],
            ltc3['class_temp','distance'][np.where(ltc3['class_temp']!='?')]],['star_basic','ltc3'],distance_cut_in_pc)
