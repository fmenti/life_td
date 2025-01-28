from utils.io import load, Path
from utils.analysis import different_data
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.table import Table

#not showing plots
show_kw=True #<--- added
#show_kw=False
    
if show_kw:   #<--- added
    plt.ion() #<--- added

def approxvalue(dist, val,x):
    poly = np.polyfit(dist,val,1)
    return poly[0]*x+poly[1]

with open(Path().data+'distance_cut.txt', 'r') as h:  
    distance_cut=float(h.readlines()[0])

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
    [table]=load(['mes_mass_st'],location=Path().data)
    arr_with_potential_fill_values=table['mass_st_value']
    data=different_data(arr_with_potential_fill_values)
    bins=10

    fig=plt.figure()
    y, bins2,patches=plt.hist(data, bins, density=True)
    #careful density=True means probability density (area under histogram integrates to 1
    plt.close(fig)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, bins)

    def model_f(x,a,b,c):
        return a * np.exp(-b * x) + c

    popt, pcov = curve_fit(model_f, x, y, p0=[1,2,0])
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_f(x_model, a_opt, b_opt, c_opt) 
 
    plt.figure()
    plt.title('mes_mass_st')
    plt.scatter(x,y)
    plt.plot(x_model, y_model, color='r')
    plt.show()

    # values for dist were tested before
    a = approxvalue([30,5],[0.5,3.0],distance_cut)
    b = approxvalue([30,5],[1,18],distance_cut)
    c = approxvalue([30,5],[-0.1,0.2],distance_cut)

    #assert
    assert a_opt <a+a/10. and a_opt >a-a/10.
    assert b_opt <b+b/10. and b_opt >b-b/10.
    assert c_opt <c+c/10. and a_opt >c-c/10.

def test_data_makes_sense_mass_pl():
    #data
    [table]=load(['mes_mass_pl'],location=Path().data)
    data=table['mass_pl_value']

    bins=10
    counts, bin_edges = np.histogram(data,bins)

    y=counts

    x = np.linspace(bin_edges[0],bin_edges[-1],bins)

    def model_f(x,a,b,c):
        return a*x**(-b)+c

    popt, pcov = curve_fit(model_f, x, y, p0=[1.5,0.5,1.5])
    a_opt, b_opt, c_opt = popt
    x_model = np.linspace(min(x), max(x), 100)
    y_model = model_f(x_model, a_opt, b_opt, c_opt) 
 
    plt.figure()
    plt.title('mes_mass_pl')
    plt.scatter(x,y)
    plt.plot(x_model, y_model, color='r')
    plt.show()    

    #assert
    assert a_opt <3 and a_opt >1
    assert b_opt <1 and b_opt >0
    assert c_opt <2 and a_opt >1

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

    m = approxvalue([30,5],[10,8],distance_cut)
    s = approxvalue([30,5],[4,3.8],distance_cut)
    
    assert mu <m+m/10. and mu > m-m/10.
    assert std <s+s/10. and std >s-s/10.


    
 

