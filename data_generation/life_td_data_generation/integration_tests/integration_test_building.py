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

def model_exp_decay(x,a,b,c):
    return a * np.exp(-b * x) + c

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
    
    bins=10

    fig=plt.figure()
    y, bins2,patches=plt.hist(data, bins, density=True)
    #careful density=True means probability density (area under histogram integrates to 1
    plt.close(fig)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, bins)

    popt, pcov = curve_fit(model_exp_decay, x, y, p0=[1,8,0])
    print(popt)
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt) 
 
    plt.figure()
    plt.title('mes_mass_st')
    plt.scatter(x,y)
    plt.plot(x_model, y_model, color='r')
    plt.show()

    # values for dist were tested before
    a = linfit([30,5],[1.2,3.0],distance_cut)
    b = linfit([30,5],[7.9,17.4],distance_cut)
    c = linfit([30,5],[-0.02,0.12],distance_cut)

    #assert
    for val,val_opt in zip([a,b,c],[a_opt,b_opt,c_opt]):
        assert val-abs(val/10.) < val_opt
        assert val_opt <val+abs(val/10.) 

def test_data_makes_sense_mass_pl():
    #data
    data=get_data('mes_mass_pl','mass_pl_value')

    bins=10
    fig=plt.figure()
    y, bins2,patches=plt.hist(data, bins)#, density=True) #had to do this to solve overflow encountered in exp
    #careful density=True means probability density (area under histogram integrates to 1
    plt.close(fig)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, bins)

    popt, pcov = curve_fit(model_exp_decay, x, y, p0=[1,1,0])
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt) 
 
    plt.figure()
    plt.title('mes_mass_pl')
    plt.scatter(x,y)
    plt.plot(x_model, y_model, color='r')
    plt.show()    

    # values for dist were tested before
    a = linfit([30,5],[960,74],distance_cut)
    b = linfit([30,5],[19,36],distance_cut)
    c = linfit([30,5],[-0.24,1.6],distance_cut)

    #assert
    for val,val_opt in zip([a,b,c],[a_opt,b_opt,c_opt]):
        assert val-abs(val/10.) < val_opt
        assert val_opt <val+abs(val/10.)

def test_data_makes_sense_coo():
    #data
    [table]=load(['star_basic'],location=Path().data)
    data_ra=table['coo_ra']
    data_dec=table['coo_dec']

    arr_ra=different_data(data_ra)
    arr_dec=different_data(data_dec)

    assert min(arr_ra)>0
    assert max(arr_ra)<360
    assert min(arr_dec)>-90
    assert max(arr_dec)<90

def test_data_makes_sense_coo_gal():
    #data
    [table]=load(['star_basic'],location=Path().data)
    data_l=table['coo_gal_l']
    data_b=table['coo_gal_b']

    arr_l=different_data(data_l)
    arr_b=different_data(data_b)

    assert min(arr_l)>0
    assert max(arr_l)<360
    assert min(arr_b)>-90
    assert max(arr_b)<90

def test_data_makes_sense_mag_i():
    #data
    [table]=load(['star_basic'],location=Path().data)
    data=table['mag_i_value']
    
    arr=different_data(data)

    bins=10
    y, bins2,patches=plt.hist(arr, bins, density=True, alpha=0.6, color='g')

    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, bins)  

    mu, std = norm.fit(arr)
    y_fit = norm.pdf(x, mu, std)

    plt.plot(x, y_fit, color='r')
    plt.show()

    m = linfit([30,5],[10,8],distance_cut)
    s = linfit([30,5],[3.4,3.8],distance_cut)
    
    assert mu <m+m/10. and mu > m-m/10.
    assert std <s+s/10. and std >s-s/10.


    
 

