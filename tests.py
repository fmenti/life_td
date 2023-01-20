import matplotlib.pyplot as plt
import numpy as np #arrays
import pyvo as vo #catalog query
import astropy as ap #votables

#code snippet not running outside of environment of data_acquisition.py
def sanitytest(cat,colname):  
    arr=cat[colname]
    plt.figure
    plt.xlabel(f'{colname}')
    plt.ylabel('Number of objects')
    plt.hist(arr,histtype='bar',log=True)
    plt.show()
    return

cats=[star_basic,disk_basic,mesDist,mesMass]
colnames=[['coo_ra','coo_dec','coo_err_angle','coo_err_maj','coo_err_min',
                'plx_value','plx_err'],
              ['rad_value','rad_err'],
                ['dist_value','dist_err'],
               ['mass_val','mass_err']]
for i in [0,1,2,3]:
    cat=cats[i]
    for colname in colnames[i]:
        sanitytest(cat,colname)