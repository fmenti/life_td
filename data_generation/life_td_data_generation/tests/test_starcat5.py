from catalog.starcat5 import (
    choose_service,
    unresolved_binaries,
)
from astropy.table import Table
import numpy as np

def test_choose_service():
    assert choose_service("")=="http://localhost:8080/tap"
    assert choose_service("heid") == "http://dc.zah.uni-heidelberg.de/tap"
    assert choose_service("gvo") == "http://dc.g-vo.org/tap"


def test_unresolved_binaries():
    # data
    #still need to add child_main_id mask -> actually no, that comes from join

    systems = Table((np.array([1,2]),
                     np.array(['system1','system2']),
                     np.array(['True','True'])),
                names=('object_id','main_id', 'binary_flag'),
                dtype=[int, object, object])
    children = Table((np.array([1,1]),
                     np.array(['star1','system2']),
                     np.array(['st','sy'])),
                names=('parent_object_idref','child_main_id','child_type'),
                dtype=[int, object, object])
    stars = Table((np.array(['star1','star2']),
                     np.array(['system1','']),
                     np.array(['True','False'])),
                names=('main_id', 'parent_main_id', 'binary_flag'),
                dtype=[object, object, object])
    ub_columns=stars.colnames+['unresolved_binaries']

    # execute
    ub = unresolved_binaries(systems, children, stars)

    assert ub_columns == ub.colnames
    assert 'system2' in ub['main_id']
    assert ub['unresolved_binaries'][np.where(ub['main_id']=='star1')]=='False'
    assert ub['unresolved_binaries'][np.where(ub['main_id']=='system2')]=='True'




