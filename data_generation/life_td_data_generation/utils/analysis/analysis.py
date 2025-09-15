from dataclasses import dataclass
import numpy as np #arrays
from astropy import io
from astropy.table import MaskedColumn
import matplotlib.pyplot as plt
from enum import Enum
from scipy.optimize import curve_fit
from scipy.stats import norm
from astropy.io.ascii import read

#self created modules
from utils.io import stringtoobject, Path
from sdata import empty_dict


def linfit(x_data, y_data,x):
    poly = np.polyfit(x_data,y_data,1)
    return poly[0]*x+poly[1]

def fitfunction(name,x,y,arr):
    if name=='model':
        def model_f(x,a,b,c):
            #return a*(x-b)**(2)+c
          return a * np.exp(-b * x) + c
          #return a*x**(-b)+c
        popt, pcov = curve_fit(model_f, x, y, p0=[ 2.,   1, 1.])
        a_opt, b_opt, c_opt = popt
        x_model = np.linspace(min(x), max(x), 100)
        y_model = model_f(x_model, a_opt, b_opt, c_opt) 
        print(popt)
        #assert a_opt <3 and a_opt >1
        #assert b_opt <1 and b_opt >0
        #assert c_opt <2 and a_opt >1
        return x_model,y_model
    elif name == 'poly':
        poly=np.polyfit(x,y,2)
        y_fit=np.polyval(poly, x)
        return x,y_fit
    elif name == 'norm':
        mu, std = norm.fit(arr)
        y_fit = norm.pdf(x, mu, std)
        print(mu,std)
        #assert mu <10 and mu >9
        #assert std <4 and std >3
        return x,y_fit
    
empty=empty_dict.copy()
table_names=list(empty.keys())

class SpectralType(Enum):
    O = 'O'
    B = 'B'
    A = 'A'
    F = 'F'
    G = 'G'
    K = 'K'
    M = 'M'

def prepare_table(cat,columns):
    if columns==[]:
        table=cat
    else:
        table=cat[columns]
    return table

def match_table(table,wherecol,whereobj):
    if wherecol!='' or whereobj!='':
        table=table[np.where(table[wherecol]==whereobj)]        
    return table

def show(provider,
         table='objects',
         columns=[],
         wherecol='',
         whereobj='',
         notobjs=[]):
    """
    This function prints the specified columns of a table.
    
    :param provider: database_tables, sim, wds, gaia, gk, life
    :param table:
    :param columns: if empty prints all
    :param wehercol:
    :param whereobj:
    """
    
    cat=provider[table_names.index(table)]
    table=prepare_table(cat,columns)
    
    table=match_table(table,wherecol,whereobj)
    #at one point also excluding some  
    #if notobjs!=[]:
    #    for element in notobjs:
    #        table=match_table(table,wherecol,whereobj)
    print(table)    
    return

def is_starnames(arr):
    return type(arr[0])==np.str_
def array_only_fill_values(arr):
    return len(arr[np.where(arr!=1e20)])==0
def remove_fill_values(arr):
    return arr[np.where(arr!=1e20)]

def different_data(arr):
    if is_starnames(arr):
        pass   
    else:
        if len(arr)==0 or array_only_fill_values(arr):
            pass
        if max(arr)==1e20:
            arr=remove_fill_values(arr)
    return arr
    
#TBD: durch funktionalen test ersetzen    
def sanitytest(cat,colname):
    """
    Prints histogram of parameter for quick sanity test.
    
    :param cat: Table containing column colname.
    :type cat: astropy.table.table.Table
    :param str colname: Column name
    """
    
    arr=cat[colname]
    # because python and np disagree on arr!=1e20 result
    
    different_data(arr)

    plt.figure()
    plt.xlabel(f'{colname}')
    plt.ylabel('Number of objects')
    plt.hist(arr,histtype='bar',log=True)
    plt.show()
    return



def stellar_distance_histogram(arr,names,path,max_dist,xaxis):
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

def spectral_type_histogram(spectypes):
    '''
    Makes a histogram of the spectral type distribution.
    
    :param spectypes: array containing a spectral type for each star
    :type spectypes: 
    :returns: array containig the number of stars for
        each spectral type in spec.
    :rtype: np.array
    '''
    
    specdist =np.zeros(len(SpectralType))
    for spectype in spectypes:
        if spectype not in ['','nan']:
            for j,spectraltype in enumerate(SpectralType):
                if spectype[0] == spectraltype.name:
                    specdist[j] += 1
    specdist=specdist.astype(int)
    return specdist

def myround(x, base=5):
    return base * round(x/base)

def get_distance_cut(distance_samples):
    """
    Gets the maximum distance of the sample of catalogs.
    """
    distance_cut=0
    for i in range(len(distance_samples)):
        if max(distance_samples[i])> distance_cut:
            distance_cut=max(distance_samples[i])
    return distance_cut

def define_xticks(distance_cut: float):
    """
    Prepares the xticks of the second subfigure in the final_plot.
    
    :param float distance_cut: Max distance of sample in parsec.
    :returns:
    :rtype: 
    """
    
    if distance_cut<=35.:
        xticks_total=['0-5','5-10','10-15','15-20','20-25','25-30']
        stepsize=5.
    elif distance_cut<=50.:
        xticks_total=['0-10','10-20','20-30','30-40','40-50']
        stepsize=10.
    else:
        print('Error: distance bigger than 50pc not programmed.')
    
    xticks=xticks_total[:int(distance_cut/stepsize)]   
    return xticks, stepsize

@dataclass
class Plotparas:
    width = 0.15
    color = ['tab:blue','tab:cyan','tab:green','tab:olive']

def x_position(x,n_samples,sample_index):
    """
    Calculates x position for plots so that not all samples plotted over each other.
    """
    width_of_samples=n_samples*Plotparas.width
    sample_location=sample_index*Plotparas.width
    # if I change sign in front of sample_location result is that catalogs get shown in order backwards
    return x-width_of_samples/2+sample_location

def tight_plot(x,spec):
    return x[1:],spec[1:]

def full_subplot(x,spec,i,labels):
    x1,y1=tight_plot(x,spec)
    plt.bar(x1,y1,Plotparas.width, align='center',
                color=Plotparas.color[i],log=True, edgecolor='black',
                label=labels[i])
    return

def distancecut_subplot(x,spec,i,hatch):
    x1,y1=tight_plot(x,spec)
    plt.bar(x1,y1,Plotparas.width, align='center',
                color=Plotparas.color[i],log=True, edgecolor='black',
                hatch=hatch)
    return

def spectral_type_histogram_catalog_comparison(stellar_catalogs,labels,path=Path().plot+'sthcc.png'):
    """
    Creates the figure for the RNAAS article.
    
    Spectral type distribution of catalogs with shading for amount below 20pc.
    """
    spectral_type_samples=[stellar_cat['spec'] for stellar_cat in stellar_catalogs]
        
    #ititializes the data containers
    spec=[spectraltype.name for spectraltype in SpectralType]
    #initializes array to contain spectral distribution of total samples
    specdist_total=np.zeros((len(spectral_type_samples),len(spec)))
    #initializes array to contain spectral distribution of sub samples
    specdist_sub=np.empty_like(specdist_total)
    
    plt.figure(figsize=(8, 4))
    x=np.arange(len(spec))
    
    
    for i in range(len(spectral_type_samples)):
        sample_total=stellar_catalogs[i]['spec']
        specdist_total[i]=spectral_type_histogram(sample_total)
                
        x_pos=x_position(x,len(spectral_type_samples),i)
        
        full_subplot(x_pos,specdist_total[i],i,labels)

        #now make same for only within 20pc sample
        sample_sub=sample_total[np.where(stellar_catalogs[i]['dist']<20.)]
        specdist_sub[i]=spectral_type_histogram(sample_sub)
        
        distancecut_subplot(x_pos,specdist_sub[i],i,'//')
        
        sample_sub=sample_total[np.where(stellar_catalogs[i]['dist']<10.)]
        specdist_sub[i]=spectral_type_histogram(sample_sub)
        
        distancecut_subplot(x_pos,specdist_sub[i],i,'////')

        
    #creating a single label for the tree hatch barplots  
    plt.bar([1],[0],log=True,edgecolor='black',hatch='//',label=['d < 20pc'],facecolor="none")  
    plt.bar([1],[0],log=True,edgecolor='black',hatch='////',label=['d < 10pc'],facecolor="none") 
    
    xt,yt=tight_plot(x,spec)
    plt.xticks(xt,yt)
    plt.title(f"Spectral type distribution")
    plt.ylabel('Number of stars')
    plt.xlabel('Spectral types')
    plt.legend(loc='upper left')
    plt.savefig(path, dpi=300)
    return

def final_plot(stars,labels,path=Path().plot+'final_plot.png'):
    """
    Plots spectral distribution in two subplots.
    
    This plot is a bit of a mess, in the future use the function
    spectral_type_histogram_catalog_comparison instead.
    Makes plot with two subfigures, first a histogram of spectral 
    type and then one of spectral type for each distance sample of 
    0-5, 5-10, 10-15 and 15-20pc.
    
    :param stars: list of astropy.table.table.Table of shape (,2) 
        with columns first spectral type and then distance in pc.
    :param labels: list containing the labels for the plot
    :param path: location to save the plot
    """
    #tbd: replace stars with two variables to prevent having [0] hardcoded stuff
    # list of arrays of spectral type info: spectral_type_samples
    # list of arrays of spectral type info: distance_samples
    
    spectral_type_samples=[stellar_cat[stellar_cat.colnames[0]] for stellar_cat in stars]
    distance_samples=[stellar_cat[stellar_cat.colnames[1]] for stellar_cat in stars]
    
    
    distance_cut_in_pc=myround(get_distance_cut(distance_samples))
    xticks, stepsize=define_xticks(distance_cut_in_pc)
    
    #variables for both subfigures
    width = 0.15 #with of the bars
    color=['tab:blue','tab:orange','tab:green']
    
    #ititializes the data containers
    spec=[spectraltype.name for spectraltype in SpectralType]
    specdist=np.zeros((len(spectral_type_samples),len(spec)))
    
    #prepares subfigure display
    fig, ((ax1,ax2)) = plt.subplots(1,2,sharey='row',\
                       figsize = [16, 5],\
                       gridspec_kw={'hspace': 0, 'wspace': 0})  
    
    #variables for first subfigure
    x=np.arange(len(spec))
    
    
    #first subfigure
    #remove spectraltypes that are not present (O and B) for tighter plot
    # => x -> x[2:], specdist[i] -> specdist[i][2:]
    
    #first_sub_figure(ax1)
    for i in range(len(spectral_type_samples)):
        specdist[i]=spectral_type_histogram(spectral_type_samples[i])
        #shift x so that not all samples plotted over each other
        ax1.bar(x_position(x[2:],len(spectral_type_samples),width,i),
                specdist[i][2:],
                width, align='center',label=labels[i],color=color[i])

    ax1.set_yscale('log') 
    #problem beim zweiten plot: min() arg is an empty sequence
    plt.sca(ax1)
    plt.xticks(x[2:], spec[2:])
    ax1.set_title(f"Spectral type distribution")
    ax1.set_ylabel('Number of stars')
    ax1.set_xlabel('Spectral types')
    ax1.legend(loc='upper left')

    #second subfigure 
    index = np.arange(len(xticks)) #wie viele bars pro label es haben wird
    #n = np.arange(7)#wieviele labels
    sub_specdist=np.zeros((len(spectral_type_samples),len(xticks),len(spec)))
    #def distance_binned_spectral_distributions(spectral_type_samples,distance_samples):
        
        #return sub_specdist
           
    for i, spectral_type_sample in enumerate(spectral_type_samples):
        #what am I doing here? I create the distance binned spectral distributions
        for j in range(len(xticks)):
            #here have an error because lifestarcat only goes up to 20pc
            if j*stepsize> max(stars[i][1][:]):
                sub_specdist[i][j]=[0.]*7
            else: 
                sub_specdist[i][j]=spectral_type_histogram(
                        stars[i][np.where((stars[i][1][:] >j*stepsize)*\
                        (stars[i][1][:] < (j+1)*stepsize))][0][:])
        #here I put in the axis bars
        for l in np.arange(len(spec))[2:]:
            #print(sub_specdist[i][:,l],
            #5*index-((n+(n+2)*(s-2))*width)+((n+2)*l+i)*width) 
            #problem if in one bin all 0
            #make if clause for it
            
            ax2.bar(x_position((len(xticks)+1)*index,len(stars)+ (len(stars)+2) * (len(spec)-2),width, (len(stars)+2)*l+i),\
                    tuple(sub_specdist[i][:,l]),width,color=color[i])
         
    ax2.set_xlabel('Distance [pc]')
    plt.sca(ax2)
    plt.xticks((len(xticks)+1)*index, (xticks))#(xticks_name))
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
    spec=[spectraltype.name for spectraltype in SpectralType]
    s=len(spec)
    specdist=np.zeros((n,s))
    
    fig, ax = plt.subplots()
    x=np.arange(s)
    width = 0.15

    for i in range(n):
        specdist[i]=spectral_type_histogram(stars[i])
        ax.bar(x[2:]-n*width/2+i*width,specdist[i][2:],
               width,label=name[i])
    
    ax.set_yscale('log')
    ax.set_ylabel('Number of objects')
    ax.set_title('Spectral type distribution')
    ax.set_xticks(x[2:])
    ax.set_xticklabels(spec[2:])
    ax.legend()
    fig.tight_layout()
    plt.savefig(Path().plot+path, dpi=300)
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
    
    plt.figure()
    plt.title('Object type distribution')
    plt.xlabel('Number of objects')
    plt.hist(cat,histtype='bar',log=True,orientation='horizontal')
    plt.yticks(np.arange(4),spec)
    plt.savefig(Path().plot+path, dpi=300)
    plt.show()
    return

def sanity_tests(database_tables, distance_cut_in_pc=30.,StarCat3=False):
    """
    Performs some sanity tests.
    
    TBD: split into four functions that can be run separately.
    
    :param database_tables:
    :type database_tables:
    :param float distance_cut_in_pc:
    :param bool StarCat3: Dafaults to False. If true plots StarCat3.
    """
    
    data=database_tables
    #data=exo
    print('looking at table data and metadata \n')
    for i in range(len(table_names)):
        print(table_names[i],i,'\n')
        print(data[i].info(['attributes', 'stats']))
        print(data[i][0:5],'\n','\n')

    print('looking at plots')

    ###problematic if I change order of database_tables tables. need to 
    #make it independent of that using table_names list
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

    stellar_distance_histogram(database_tables[5]['dist_st_value'],['star_basic'],'test',
            distance_cut_in_pc,'dist_st_value')
    
    table=database_tables[5]['class_temp','dist_st_value',
                                'binary_flag'][np.where(
                                    database_tables[5]['class_temp']!='?')]
    table=table['class_temp','dist_st_value'][np.where(
            table['binary_flag']=='False')]
            
    if StarCat3:
        ltc3=io.ascii.read(Path().additional_data+"LIFE-StarCat3.csv")
        ltc3=stringtoobject(ltc3,3000)
        print(ltc3['distance'])
        ltc3['class_temp']=MaskedColumn(dtype=object,length=len(ltc3))
        for i in range(len(ltc3)):
            #sorting out entries like '', DA2.9, T1V
            if len(ltc3['sim_sptype'][i])>0 and ltc3['sim_sptype'][i][0] in \
                     ['O','B','A','F','G','K','M']:
                ltc3['class_temp'][i]=ltc3['sim_sptype'][i][0]
            else:
                ltc3['class_temp'][i]='?'
                
        final_plot([table,ltc3['class_temp','distance'][np.where(
                ltc3['class_temp']!='?')]],
                ['LTC4_singles','LTC3'],distance_cut_in_pc)
    else:
        final_plot([table],['single_stars'],distance_cut_in_pc)


    table=database_tables[5]['class_temp','dist_st_value',
            'binary_flag'][np.where(
                database_tables[5]['class_temp']!='?')]
    singles=table['class_temp'][np.where(table['binary_flag']=='False')]
    multiples=table['class_temp'][np.where(table['binary_flag']=='True')]
    total=table['class_temp']
    spechistplot([singles, multiples,total],['singles','multiples','total'])


    table=database_tables[1]['type']
    objecthistplot([table],['objects'],'objecthist')
    return
