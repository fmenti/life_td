from catalog.starcat5 import (
    choose_service,
    add_unresolved_binaries,
    flag_non_main_sequence_stars,
)
from astropy.table import Table, MaskedColumn
import numpy as np
from utils.analysis.analysis import match_table

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
    ub = add_unresolved_binaries(systems, children, stars)

    assert ub_columns == ub.colnames
    assert 'system2' in ub['main_id']
    assert ub['unresolved_binaries'][np.where(ub['main_id']=='star1')]=='False'
    assert ub['unresolved_binaries'][np.where(ub['main_id']=='system2')]=='True'

def test_flag_non_main_sequence_stars():
    data = Table((np.array(['ms_star',
                            'no_temp_class','no_lumclass','no_mass',
                            'wrong_temp_class','wrong_lumclass']),
                  MaskedColumn(np.array(['M','','K','G','Y','F']),
                               mask=[False,True,False,False,False,False]),
                  MaskedColumn(np.array(['V','V','','V','V','VI']),
                               mask=[False,False,True,False,False,False]),
                  MaskedColumn(np.array([0.5,1.,2.,999,1.5,3.]),
                               mask=[False,False,False,True,False,False])),
                names=('main_id', 'class_temp', 'class_lum','mass_st_value'),
                dtype=[object, object, object, float])

    catalog = flag_non_main_sequence_stars(data)
    ids = [['no_temp_class','wrong_temp_class'],
           ['no_lumclass','wrong_lumclass'],['no_mass']]
    for id,col in zip(ids,
                      ['ms_temp_class','ms_lum_class','mass_flag']):
        for identifier in id:
            if identifier == 'no_mass':
                false=False
            else:
                false ='False'
            assert match_table(catalog,'main_id',
                               identifier)[col]==false

    #now check that columns have only a certain amount of trues
    assert len(match_table(catalog,'ms_temp_class','True'))==4
    assert len(match_table(catalog,'ms_lum_class','True'))==4
    assert len(match_table(catalog,'mass_flag',True))==5

