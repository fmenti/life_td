import numpy as np  # Used for arrays
import astropy as ap  # Used for votables
from utils.starcat4 import starcat_creation

def test_star_cat_looks_fine():

    starcat4 = starcat_creation(30)
    # why do I have difficulties running
    # this function and why do I not get the print I want? because was called integration test instead of test...
    #print(starcat4)
    K_in_30=starcat4[np.where(starcat4['class_temp']=='K')]
    K_in_20=K_in_30[np.where(K_in_30['dist_st_value']<20.)]
    print(K_in_20) #126, felix reports 154 which is concerning because in ltc3 we had 245
    #find out which ones are missing and why, maybe start with G stars as there are fewer
    G_in_30 = starcat4[np.where(starcat4['class_temp'] == 'G')]
    G_in_20 = G_in_30[np.where(G_in_30['dist_st_value'] < 20.)]
    print(len(G_in_20)) #34 felix reports 35 and in ltcv 65

    assert type(starcat4) == type(ap.table.Table())

def test_star_cat_contains_specific_objects():

    starcat4 = starcat_creation(30)

    assert len(starcat4[np.where(starcat4['main_id'] == 'TRAPPIST-1')]) == 1

