from catalog.starcat5 import (
    choose_service,
    add_unresolved_binaries,
    flag_non_main_sequence_stars,
    sorting_number_of_id,
    flag_trivial_binaries,
    deal_with_separation,
    assign_critical_separation,
    crit_sep,
    apply_stability_constraint,
    flag_ecliptic,
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

def test_sorting_number_of_id():
    input_column = np.array(['test','test','test2','test3'])
    match_column = np.array(['test','test2','test3','test4','test'])

    flag_array = sorting_number_of_id(input_column,1,match_column)
    print(flag_array)
    assert match_column[flag_array][0]=='test2'
    assert match_column[flag_array][1] == 'test3'
    assert len(match_column[flag_array])==2

def test_flag_trivial_binaries():
    #data
    main_id = np.array(['star1','star2','single_star',
                               'starA','multiparent_star','multiparent_star',
                               'tripleA','tripleB','tripleC'])
    parent_main_id = np.array(['system1','system1','',
                               'systemAB','parent1','parent2',
                               'triple_parent','triple_parent','triple_parent'])
    catalog = Table((main_id,parent_main_id,
                     np.array(['True','True','False','True','True','True',
                               'True','True','True'])),
                names=('main_id', 'parent_main_id', 'binary_flag'),
                dtype=[object, object, object])
    children = Table((np.array(['star1','star2','systemAB','star',
                                'tripleA','tripleB','tripleC']),
                     np.array(['st','st','sy','st','st','st','st'])),
                names=('child_main_id','child_type'),
                dtype=[object, object])

    singles, multiples = flag_trivial_binaries(catalog,children)

    cols = ['higher_order_multiples','single_parent','trivial_binaries']
    assert multiples.colnames == catalog.colnames + cols
    assert singles.colnames == catalog.colnames
    assert len(singles)==1
    assert len(singles[np.where(singles['binary_flag']=='True')])==0
    assert len(multiples[np.where(multiples['binary_flag']=='False')])==0
    assert match_table(multiples,'higher_order_multiples',
                       True)['main_id']=='starA'
    assert match_table(multiples, 'single_parent',
                       False)['main_id'][0] == 'multiparent_star'
    assert len(match_table(multiples, 'single_parent',
                       False)) == 2
    assert len(match_table(multiples, 'trivial_binaries',
                       True)) == 2

def test_deal_with_separation():
    # data
    multiples = Table((MaskedColumn([1.,100.,np.ma.masked]),
                            np.array([True,True, False]),
                            np.array([5.,20.,10.])),
                names=('sep_ang_value','sep_flag','dist_st_value'),
                dtype=[float, bool, float])

    # execute
    result = deal_with_separation(multiples)

    # assert
    assert result.colnames == ['sep_ang_value','sep_flag','dist_st_value','sep_phys_value']
    assert result['sep_phys_value'][0] == 5.
    assert result['sep_phys_value'][1] == 2000.
    assert result['sep_phys_value'].mask[2] == True

def test_assign_critical_separation():
    # data
    multiples = Table((np.array(['B_parent', 'A_parent','A_parent',
                                 'B_parent','C_parent']),
                       np.array([True, True, True, True, False]),
                       np.array([5., 20., 10.,3., 1.]),
                       np.array([1., 2., 3.,4.,5.])),
                      names=('parent_main_id', 'suitable_companions',
                             'sep_phys_value', 'mass_st_value'),
                      dtype=[object, bool, float,float])

    # execute
    result, HZstability = assign_critical_separation(multiples)

    # assert
    crit_sep_B_mass1 = crit_sep(0,4/(1+4),5)[0]
    row_B_mass_1 = match_table(HZstability,'mass_st_value',1.)

    assert result.colnames == ['parent_main_id', 'suitable_companions',
                             'sep_phys_value', 'mass_st_value']
    assert HZstability.colnames == ['parent_main_id', 'suitable_companions',
                               'sep_phys_value', 'mass_st_value', 'a_crit_s']
    assert len(HZstability) == 4
    assert result['parent_main_id'][0] == 'A_parent'
    assert HZstability['parent_main_id'][0] == 'A_parent'
    assert row_B_mass_1['a_crit_s'][0] == crit_sep_B_mass1

def test_apply_stability_constraint():
    # data
    HZstability = Table((np.array(['A_parent', 'A_parent','B_parent',
                                 'B_parent']),
                       np.array(['star1','star2','star3','star4']),
                       np.array([25., 20., 19.,3.])),
                      names=('parent_main_id', 'main_id','a_crit_s'),
                      dtype=[object, object, float])

    # execute
    result = apply_stability_constraint(HZstability,a_max=10.0)

    # verify
    assert len(result)==2
    assert result['parent_main_id'][0] == 'A_parent'

def test_flag_ecliptic():
    StarCat5 = Table((np.array(['star1','star2','star3','star4']),
                      np.array([0., 0., 90., 90.]),
                       np.array([0., 50., 0.,-25.])),
                      names=('main_id', 'coo_ra','coo_dec'),
                      dtype=[object, float, float])
    angle = 45
    StarCat5[f'ecliptic_pm{angle}deg'] = flag_ecliptic(
        angle, StarCat5["coo_ra"], StarCat5["coo_dec"]
    )
    assert 'ecliptic_pm45deg' in StarCat5.colnames
    assert StarCat5['ecliptic_pm45deg'][0]=='True'
    assert StarCat5['ecliptic_pm45deg'][1]=='False'
    assert StarCat5['ecliptic_pm45deg'][2]=='True'
    assert StarCat5['ecliptic_pm45deg'][3]=='False'
