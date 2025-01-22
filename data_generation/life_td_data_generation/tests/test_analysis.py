from utils.analysis import *
from astropy.table import Table



# TODO wirte a test for each cat data
def test_different_data_standard():
    cat = Table(data=[[39, 40, 2.5],
                      [1e20, 90, 81]],
                names=['coo_dec','coo_err_angle'], 
                dtype=[float, float])

    colname = 'coo_dec'
    arr=cat[colname]
    
    arr=different_data(arr)

    assert len(arr) == 3
    assert arr[0] == 39
    assert max(arr) == 40
    
def test_different_data_fill_values():
    cat = Table(data=[[39, 40, 2.5],
                      [1e20, 90, 81]],
                names=['coo_dec','coo_err_angle'], 
                dtype=[float, float])

    colname = 'coo_err_angle'
    arr=cat[colname]
    
    arr=different_data(arr)

    assert len(arr) == 2
    assert arr[0] == 90
    assert max(arr) == 90
    


    #next: do same like this function again for the different data columns so that all diff_data cases get tested.

#def test_2():
#    assert False
