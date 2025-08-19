from provider.simbad import *
from astropy.table import Table, setdiff, MaskedColumn
import provider.simbad as simbad_module
from sdata import empty_dict
import pytest
from utils.io import stringtoobject

@pytest.fixture
def query_returns(monkeypatch):
    # Shared setup for tests that need the helper table
    distance_cut_in_pc = 5
    test_objects = []

    membership = MaskedColumn([np.ma.masked,np.ma.masked,100], name='membership', mask=[True, True,False])
    plx_qual = MaskedColumn(data=['A','A','C'], mask=[False, False, False], fill_value='N')
    t0 = Table({
        'main_id': np.array(['star1','star2','star3'],dtype=object),
        'coo_ra': [160.,165.,100.],
        'coo_dec': [35.,45.,30.],
        'coo_err_angle': [90.,17.,np.ma.masked],
        'coo_err_maj': [0.02,180.,np.ma.masked],
        'coo_err_min': [0.01,160.,np.ma.masked],
        'oid': [1,2,3],
        'coo_ref': np.array(['ref_coo','ref_coo','ref_coo2'],dtype=object),
        'coo_qual': ['A','A','C'],
        'sptype_string': np.array(['M7.0V','','M5.0V+M9'],dtype=object),
        'sptype_qual': ['C','','B'],
        'sptype_ref': np.array(['ref_sptype','','ref_sptype2'],dtype=object),
        'plx_err': [0.05,0.9,3],
        'plx_value': [190.,230.,375.],
        'plx_ref': np.array(['ref_plx','ref_plx','ref_plx2'],dtype=object),
        'plx_qual': plx_qual,
        'membership': membership,
        'parent_oid': [4,np.ma.masked,4],
        'h_link_ref': np.array(['ref_hlink','','ref_hlink'],dtype=object),
        'otypes': np.array(['*|**','BD?','*|**'],dtype=object),
        'ids': np.array(['star1|id1','star2|id2','star3|id3'],dtype=object),
        'mag_i_value': [np.ma.masked,np.ma.masked,7.],
        'mag_j_value': [8.,np.ma.masked,6.],
        'mag_k_value': [9.,np.ma.masked,8.]
    })
    parents_without_plx = t0[:0].copy()
    parents_without_plx.add_row(['system1',np.ma.masked,np.ma.masked,np.ma.masked,np.ma.masked,np.ma.masked,4,
                                 '','','','','',np.ma.masked,np.ma.masked,'','',np.ma.masked,np.ma.masked,
                                 '','**','system1',np.ma.masked,np.ma.masked,np.ma.masked])
    children_without_plx = t0[:0].copy()
    children_without_plx.add_row(['planet1',np.ma.masked,np.ma.masked,np.ma.masked,np.ma.masked,np.ma.masked,5,
                                 '','','','','',np.ma.masked,np.ma.masked,'','',np.ma.masked,3,
                                 'ref_hlink2','Pl','planet1',np.ma.masked,np.ma.masked,np.ma.masked])

    ident = Table({
        'main_id': np.array(['star1','star1','star3','star3','system1','planet1'],dtype=object),
        'oid': [1,1,3,3,4,5],
        'id': np.array(['star1','id1','star3','id3','system1','planet1'],dtype=object),
    })

    tables_iter = iter([t0, parents_without_plx, children_without_plx,ident])

    def fake_query(url, adql, uploads=None):
        return next(tables_iter)

    monkeypatch.setattr(simbad_module, 'query', fake_query)

    result_sim_helptab, result_sim = create_simbad_helpertable(distance_cut_in_pc, test_objects)

    result_sim['ident'] = create_ident_table(result_sim_helptab, result_sim)

    return result_sim_helptab, result_sim, t0, result_sim['ident']


def test_create_simbad_helpertable(query_returns):
    result_sim_helptab, result_sim, t0, ident = query_returns

    # Verify
    # all helptab columns are present
    assert set(result_sim_helptab.colnames) == set(t0.colnames+['type','binary_flag'])
    assert set(result_sim_helptab['type']) == set(['sy','sy','sy','pl'])
    # objects got correctly added and removed from helptab
    assert len(result_sim_helptab) == 4# stars, systems and planets but removed star2
    # provider table is correctly filled in
    assert len(result_sim['provider']) == 1

def test_creating_helpertable_stars(query_returns):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns

    # Execute
    result_stars = creating_helpertable_stars(result_sim_helptab, result_sim)

    # Verify
    expected_stars = result_sim_helptab
    expected_stars.remove_row(-1)
    for null,colname in zip([999999,999999,999999,999999,999999,999999,999999,999999],
                            ['membership','coo_err_angle','coo_err_maj',
                             'coo_err_min','parent_oid','mag_i_value','mag_j_value','mag_k_value']):
        result_stars= nullvalues(result_stars, colname, null)
        expected_stars = nullvalues(expected_stars, colname, null)
    assert len(result_stars) == 3# stars and systems but removed planets
    assert len(setdiff(result_stars, expected_stars)) == 0

def test_create_ident_table(query_returns):
    # Execute
    result_sim_helptab, result_sim, t0, return_ident = query_returns

    # Verify
    ref = result_sim['provider']['provider_bibcode'][0]
    expected_ident = Table({
        'main_id': np.array(['star1', 'star1','star3', 'star3','system1', 'planet1'],dtype=object),
        'id_ref': np.array([ref, ref, ref, ref, ref, ref],dtype=object),
        'id': np.array(['id1', 'star1','id3', 'star3','system1', 'planet1'],dtype=object)
    })
    print(return_ident.info)
    print(expected_ident.info)
    assert len(setdiff(return_ident, expected_ident)) == 0

def test_create_h_link_table(query_returns,monkeypatch):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    result_stars = creating_helpertable_stars(result_sim_helptab, result_sim)
    main_id_table = Table({
        'main_id': np.array(['star1','star3', 'planet1'],dtype=object),
        'parent_main_id': np.array(['system1','system1', 'star3'],dtype=object),
        'h_link_ref': np.array(['ref_hlink','ref_hlink','ref_hlink2'],dtype=object),
        'membership': [np.ma.masked,100,np.ma.masked],
        'parent_oid': [1,4,3]
    })

    def fake_fetch_main_id(sim_h_link, id_creator):
        return main_id_table

    monkeypatch.setattr(simbad_module, 'fetch_main_id', fake_fetch_main_id)

    # Execute
    result_h_link = create_h_link_table(result_sim_helptab, result_sim, result_stars)
    #issue: wrong null value in membership column. find out what simbad returns, str, need to have them as object.
    # sim returns int and masked entries

    expected_h_link = Table({
        'main_id': np.array(['star1','star3', 'planet1'],dtype=object),
        'parent_main_id': np.array(['system1','system1', 'star3'],dtype=object),
        'h_link_ref': np.array(['ref_hlink','ref_hlink','ref_hlink2'],dtype=object),
        'membership': [999999, 100,999999]
    })

    # Verify
    assert len(setdiff(result_h_link, expected_h_link)) == 0

def test_stars_in_multiple_system():
    # Data
    cat = Table({
        'main_id': np.array(['system1','star1','star2'],dtype=object),
        'type': np.array(['sy','sy','sy'],dtype=object),
        'sptype_string': np.array(['M7.0V+M5','M7.0V','M5.0V+M9'],dtype=object)
    })
    sim_h_link = Table({
        'main_id': np.array(['star1','star2','planet1'],dtype=object),
        'parent_main_id': np.array(['system1','system1','star1'],dtype=object),
        'h_link_ref': np.array(['ref1','ref1','ref2'],dtype=object)
    })
    all_objects = Table({
        'type': np.array(['sy','sy','sy','pl'],dtype=object),
        'main_id': np.array(['system1','star1','star2','planet1'],dtype=object)
    })
    # Execute
    return_cat = stars_in_multiple_system(cat, sim_h_link, all_objects)

    # Verify
    expected_cat = Table({
        'main_id': np.array(['system1','star1','star2'],dtype=object),
        'type': np.array(['sy','st','sy'],dtype=object),
        'sptype_string': np.array(['M7.0V+M5','M7.0V','M5.0V+M9'],dtype=object)
    })
    assert len(setdiff(return_cat, expected_cat)) == 0

def expanded_helptab_stars(expected_stars,ref):
    expected_stars.remove_row(-1)  # removing planet object row
    # doing stars in multiple systems function
    expected_stars['binary_flag'] = ['True', 'True', 'True']
    # expected_stars['plx_qual'][np.where(expected_stars['main_id']=='system1')] = '?'
    expected_stars['plx_qual']=expected_stars['plx_qual'].astype(object)
    expected_stars['plx_ref'][np.where(expected_stars['main_id'] == 'system1')] = ref
    # for magnitude do ref

    expected_stars['mag_i_ref'] = MaskedColumn(
        data=np.array(['', ref, ''],
                      dtype=object),
        mask=[True, False, True], fill_value='N')
    expected_stars['mag_j_ref'] = MaskedColumn(
        data=np.array([ref, ref, ''],
                      dtype=object),
        mask=[False, False, True], fill_value='N')
    expected_stars['mag_k_ref'] = MaskedColumn(
        data=np.array([ref, ref, ''],
                      dtype=object),
        mask=[False, False, True], fill_value='N')
    # replacing null values in plx , coo and sptype with sim prvoder
    expected_stars['sptype_ref'][np.where(expected_stars['main_id'] == 'system1')] = ref
    expected_stars['coo_ref'][np.where(expected_stars['main_id'] == 'system1')] = ref
    expected_stars['binary_ref'] = [ref, ref, ref]
    expected_stars['binary_qual'] = ['D', 'D', 'D']
    expected_stars['type'][np.where(expected_stars['main_id'] == 'star1')] = 'st'  # stars in multiple systems function
    return expected_stars

def test_expanding_helpertable_stars(query_returns):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    stars = result_sim_helptab
    expected_stars = expanded_helptab_stars(stars.copy(),result_sim['provider']['provider_bibcode'][0])

    stars.remove_row(-1) # removing planet object row
    result_sim['h_link']= Table({
        'main_id': np.array(['star1','star3', 'planet1'],dtype=object),
        'parent_main_id': np.array(['system1','system1', 'star3'],dtype=object),
        'h_link_ref': np.array(['ref_hlink','ref_hlink','ref_hlink2'],dtype=object),
        'membership': [999999, 100,999999]
    })

    # Execute
    #isue, plx_qual needs to be masked column with null value 'N'
    #why not '?' as I will use later?
    stars['plx_qual'].mask = [False,False,True]
    result_stars = expanding_helpertable_stars(result_sim_helptab, result_sim, stars.copy())

    # Verify


    for null,colname in zip([999999,999999,999999,999999,999999,999999,999999,999999,
                             '','','','','','','',''],
                            ['membership','coo_err_angle','coo_err_maj',
                             'coo_err_min','parent_oid','mag_i_value','mag_j_value','mag_k_value',
                             'mag_i_ref','mag_j_ref','mag_k_ref','plx_qual',
                             'coo_qual','sptype_string','sptype_qual','h_link_ref']):
        result_stars= nullvalues(result_stars, colname, null)
        expected_stars = nullvalues(expected_stars, colname, null)
    assert len(setdiff(result_stars, expected_stars)) == 0

def test_create_object_table(query_returns):
    # Execute
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    stars = result_sim_helptab.copy()
    stars.remove_row(-1)  # removing planet object row

    stars['type'][np.where(stars['main_id'] == 'star1')] = 'st' #stars in multiple systems function

    return_object = create_objects_table(result_sim_helptab, stars)

    # Verify
    expected_object = Table({
        'main_id': np.array(['star1', 'star3', 'system1', 'planet1'],dtype=object),
        'type': np.array(['st', 'sy', 'sy', 'pl'],dtype=object),
        'ids': np.array(['star1|id1','star3|id3', 'system1', 'planet1'],dtype=object)
    })

    assert len(setdiff(return_object, expected_object)) == 0

def test_create_sim_sources_table(query_returns):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    stars = result_sim_helptab.copy()
    expected_stars = expanded_helptab_stars(stars.copy(),result_sim['provider']['provider_bibcode'][0])
    result_sim['h_link'] = Table({
        'main_id': np.array(['star1','star3', 'planet1'],dtype=object),
        'parent_main_id': np.array(['system1','system1', 'star3'],dtype=object),
        'h_link_ref': np.array(['ref_hlink','ref_hlink','ref_hlink2'],dtype=object),
        'membership': [999999, 100,999999]
    })

    # Execute
    return_sources = create_sim_sources_table(expected_stars,result_sim)

    # Verify
    expected_sources = Table({
        'ref': ['ref_coo', 'ref_coo2', 'ref_sptype', 'ref_sptype2',
                         'ref_plx','ref_plx2','ref_hlink','ref_hlink2',
                         result_sim['provider']['provider_bibcode'][0]],
        'provider_name': ['SIMBAD' for j in range(9)]
    })
    assert len(setdiff(return_sources, expected_sources)) == 0

    #need to add another single star to test stuff properly as I currently only have binary_flag true tones


