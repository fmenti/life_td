from utils.analysis.catalog_comparison import testobject_dropout
from astropy.table import Table

def test_testobject_dropout():
    test_objects = Table(
        data=[['name 2', 'name 3'], [1e20, 90]],
        names=["main_id", "coo_err_angle"],
        dtype=[object, float],)
    parent_sample = Table(
        data=[['name 1', 'name 2', 'name 4'], [1e20, 90, 81]],
        names=["main_id", "coo_err_angle"],
        dtype=[object, float],)

    dropout, test_without_dropout = testobject_dropout(test_objects['main_id'],parent_sample['main_id'])

    assert dropout == ['name 3']
    assert test_without_dropout == ['name 2']
