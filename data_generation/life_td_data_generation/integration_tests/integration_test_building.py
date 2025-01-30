from utils.io import load, Path
from utils.analysis import different_data, linfit
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
    [table]=load([table_name],location=Path().data)
    arr_with_potential_fill_values=table[colname]
    data=different_data(arr_with_potential_fill_values)
    return data

def get_x_and_y(data):
    bins=10
    
    plt.figure()
    y, bins2,patches=plt.hist(data, bins, density=True)
    #careful density=True means probability density (area under histogram integrates to 1
    xmin, xmax = plt.xlim()
    plt.close()
    x = np.linspace(xmin, xmax, bins)
    return x,y

def model_exp_decay(x,a,b,c):
    return a * np.exp(-b * x) + c

def plot_data_and_fit(title,x,y,x_model,y_model):
    plt.figure()
    plt.title(title)
    plt.scatter(x,y,label='data')
    plt.plot(x_model, y_model, color='r', label='fit')
    plt.legend(loc='upper right')
    plt.show()
    return

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
    
    x,y=get_x_and_y(data)

    popt, pcov = curve_fit(model_exp_decay, x, y, p0=[1,8,0])
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt) 

    plot_data_and_fit('mes_mass_st',x,y,x_model,y_model)

    #assert
    assert max(data) < 60 # O3V
    assert min(data) > 0.074 # brown dwarf

def test_data_makes_sense_mass_pl():
    #data
    data=get_data('mes_mass_pl','mass_pl_value')

    x,y=get_x_and_y(data)

    popt, pcov = curve_fit(model_exp_decay, x, y, p0=[1,1,0])
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt) 

    plot_data_and_fit('mes_mass_pl',x,y,x_model,y_model)   

    #assert
    assert max(data) < 75 # m star
    assert min(data) > 0 

def ravsdec(x_label,y_label,x,y):
    plt.figure()
    ra=np.linspace(0,360)
    plt.scatter(x,y,s=2)
    #plt.scatter(within45deg['coo_ra'],within45deg['coo_dec'],s=2)
    #ecliptic plane in equatorial coordinates
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #plt.savefig('plots/quasi_aitoff', dpi=300)
    plt.show()
    return

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

def norm_fit(data,title):
    bins=10
    y, bins2,patches=plt.hist(data, bins, density=True)#, alpha=0.6, color='g')
    #do I need to keep density=True here or can I use non normalized display?

    #xmin, xmax = plt.xlim()
    #x = np.linspace(xmin, xmax, bins)  

    mu, std = norm.fit(data)
    #y_fit = norm.pdf(x, mu, std)
    y_fit = norm.pdf(np.sort(data), mu, std)
    
    #plt.plot(x, y_fit, color='r')
    plt.title('mag_i_value')
    plt.plot(np.sort(data), y_fit, color='r')
    plt.show()
    return mu

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

    assert mu <11 and mu > 7

def test_data_makes_sense_plx():
    #data
    data=get_data('star_basic','plx_value')
    
    x,y=get_x_and_y(data)
    popt, pcov = curve_fit(model_exp_decay, x, y, p0=[0.03,0.01,0.001])
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt) 

    plot_data_and_fit('plx_value',x,y,x_model,y_model)
    #assert
    assert min(data) > 28 # 35 pc equivalent marcsec
    assert min(data) < 1000 # 1 pc equivalent marcsec



    
 

