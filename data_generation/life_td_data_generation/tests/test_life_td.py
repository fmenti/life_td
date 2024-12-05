from life_td import *
from sdata import empty_cat_wit_columns
from utils.io import save

def test_load_cat():
    cat=empty_cat_wit_columns.copy()
    save(list(cat.values()),['test_cat_'+ element for element in list(cat.keys())])
    provider_name='test_cat'
    loaded=load_cat(provider_name)
    assert loaded==cat

