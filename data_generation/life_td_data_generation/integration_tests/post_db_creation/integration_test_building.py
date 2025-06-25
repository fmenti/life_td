from utils.io import load, Path
from utils.analysis.analysis import different_data
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.table import Table

#not showing plots
#show_kw=True #<--- added
show_kw=False
    
if show_kw:   #<--- added
    plt.ion() #<--- added

with open(Path().data+'distance_cut.txt', 'r') as h:  
    distance_cut=float(h.readlines()[0])

def get_data(table_name,colname):
    # loading the correct table
    [table]=load([table_name],location=Path().data)
    # exctracting the correct columns
    arr_with_potential_fill_values=table[colname]
    # removing fill values
    data=different_data(arr_with_potential_fill_values)
    return data


def model_exp_decay(x,a,b,c):
    return a * np.exp(-b * x) + c


def plot_data_and_fit(title,data,p0):
    #trying to not fit the scatter but the histogram directly
    #issue, scatter points offset from histogram bars
    #current issue, not meaningful zoom because y_model too high
    
    plt.figure()
    plt.title(title)
    bins=10
    
    bin_heights, bin_borders,patches=plt.hist(data, bins, density=True,label='pdf of data')
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    #careful density=True means probability density (area under histogram integrates to 1
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, bins)

    popt, pcov = curve_fit(model_exp_decay, bin_centers, bin_heights, p0=p0)#1,8,0])
    #print(popt)
    a_opt, b_opt, c_opt = popt
    #use here bin_centers[0] instead of min(x) to prevent first y_model value being so high
    # that one does nearly not see any data
    x_model = np.linspace(bin_centers[0], max(x), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt) 
    
    plt.plot(x_model, y_model, color='r', label='fit')
    plt.ylabel('amount of objects')
    plt.legend(loc='upper right')
    plt.show()
    return

def ravsdec(x_label,y_label,x,y):
    plt.figure()
    ra=np.linspace(0,360)
    plt.scatter(x,y,s=2)
    #plt.scatter(within45deg['coo_ra'],within45deg['coo_dec'],s=2)
    #ecliptic plane in equatorial coordinates
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title('coordinate distribution')
    #plt.savefig('plots/quasi_aitoff', dpi=300)
    plt.show()
    return

def norm_fit(data,title):
    bins=10
    plt.hist(data, bins, density=True,label='pdf of data')#, alpha=0.6, color='g')
    #do I need to keep density=True here or can I use non normalized display?

    mu, std = norm.fit(data)
    y_fit = norm.pdf(np.sort(data), mu, std)
    
    plt.title('mag_i_value')
    plt.plot(np.sort(data), y_fit, color='r',label='fit')
    plt.legend()
    plt.show()
    return mu

def test_data_makes_sense_main_id():
    # experimental numbers
    dist_dist = Table(data=[ [    5,    10,    30],
                         [  420,  1728, 22107],
                         [  205,   931, 15869],
                         [  101,   413,  4625],
                         [  113,   366,  1464],
                         [    6,    18,   149]],
                 names=['dist','total','st','sy','pl','di'], 
                 dtype=[float, float,float,float,float,float])
    
    [table]=load(['objects'],location=Path().data)
    data=table['main_id','type']

    spec=np.array(['System','Star','Exoplanet','Disk'])
    
    plt.figure()
    plt.title(f'Object type distribution up to {distance_cut} pc')
    plt.xlabel('Number of objects')
    plt.hist(data['type'],histtype='bar',log=True,orientation='horizontal')
    plt.yticks(np.arange(4),spec)
    #plt.savefig(Path().plot+path, dpi=300)
    plt.show()

    total=len(data)**(1/3)/distance_cut
    st=len(data[np.where(data['type']=='st')])**(1/3)/distance_cut
    sy=len(data[np.where(data['type']=='sy')])**(1/3)/distance_cut
    pl=len(data[np.where(data['type']=='pl')])**(1/3)/distance_cut
    di=len(data[np.where(data['type']=='di')])**(1/3)/distance_cut
    
    assert total <2 and total > 0.5
    assert st <1.5 and st > 0.5
    assert sy <1.5 and sy > 0.4
    assert pl <1.5 and pl > 0.3
    assert di <0.4 and di > 0.1
    #tbd make this via analytical and not just experimental numbers

def test_data_makes_sense_mass_st():
    #data
    data=get_data('mes_mass_st','mass_st_value')
    
    plot_data_and_fit('mes_mass_st',data,[1,8,0])

    #assert
    assert max(data) < 60 # O3V
    assert min(data) > 0.074 # brown dwarf

def test_data_makes_sense_temp_st():
    #want a scatter plot x axis distance and y axis temperature
    #not sure what to assert
    #not sure how to get data -> what to do about fill values?
    #data
    # loading the correct table
    [table]=load(['star_basic'],location=Path().data)
    # exctracting the correct columns
    arr=table['dist_st_value','teff_st_value']
    arr2=arr[np.where(arr['dist_st_value']!=1e20)]
    data=arr2[np.where(arr2['teff_st_value']!=1e20)]
    
    #plt.figure()
    fig, ax = plt.subplots(figsize = (9, 6))
    
    ax.scatter(data['dist_st_value'],data['teff_st_value'],s=2)
    ax.set_yscale("log")
    
    ax.set_xlabel('Distance [pc]')
    ax.set_ylabel('Temperature [K]')

    #assert
    assert False


def test_data_makes_sense_mass_pl():
    #data
    data=get_data('mes_mass_pl','mass_pl_value')

    plot_data_and_fit('mes_mass_pl',data,[1,1,0]) 

    #assert
    assert max(data) < 75 # m star
    assert min(data) > 0 
    #not working for distance cut 20



def test_data_makes_sense_coo():
    #data
    data_ra=get_data('star_basic','coo_ra')
    data_dec=get_data('star_basic','coo_dec')

    #make sky plot
    ravsdec('coo_ra','coo_dec',data_ra,data_dec)

    assert min(data_ra)>0
    assert max(data_ra)<360
    assert min(data_dec)>-90
    assert max(data_dec)<90

def test_data_makes_sense_coo_gal():
    #data
    data_l=get_data('star_basic','coo_gal_l')
    data_b=get_data('star_basic','coo_gal_b')

    #make sky plot
    ravsdec('coo_gal_l','coo_gal_b',data_l,data_b)

    assert min(data_l)>0
    assert max(data_l)<360
    assert min(data_b)>-90
    assert max(data_b)<90



def test_data_makes_sense_mag_i():
    #data
    data=get_data('star_basic','mag_i_value')
     
    mu = norm_fit(data,'mag_i_value')
    
    assert mu <11 and mu > 7

def test_data_makes_sense_mag_j():
    #data
    data=get_data('star_basic','mag_j_value')
    
    mu = norm_fit(data,'mag_j_value')
    
    assert mu <11 and mu > 7

def test_data_makes_sense_mag_k():
    #data
    data=get_data('star_basic','mag_k_value')
    
    mu = norm_fit(data,'mag_k_value')

    assert mu <11 and mu > 6#failed at 5pc

def test_data_makes_sense_plx():
    #data
    data=get_data('star_basic','plx_value')

    plot_data_and_fit('plx_value',data,[0.03,0.01,0.001])

    #assert
    assert min(data) > 28 # 35 pc equivalent marcsec
    assert min(data) < 1000 # 1 pc equivalent marcsec

def test_data_makes_sense_dist_st():
    #data
    data=get_data('star_basic','dist_st_value')

    plot_data_and_fit('plx_value',data,[1,-1,0])

    #assert
    assert min(data) <  35 
    assert min(data) > 1 



    
 

