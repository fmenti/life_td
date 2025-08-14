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

    membership = MaskedColumn([None,None,100], name='membership', mask=[True, True,False])
    t0 = Table({
        'main_id': ['star1','star2','star3'],
        'coo_ra': [160.,165.,100.],
        'coo_dec': [35.,45.,30.],
        'coo_err_angle': [90.,17.,None],
        'coo_err_maj': [0.02,180.,None],
        'coo_err_min': [0.01,160.,None],
        'oid': [1,2,3],
        'coo_ref': ['ref_coo','ref_coo','ref_coo2'],
        'coo_qual': ['A','A','C'],
        'sptype_string': ['M7.0V','','dM4'],
        'sptype_qual': ['C','','B'],
        'sptype_ref': ['ref_sptype','','ref_sptype2'],
        'plx_err': [0.05,0.9,3],
        'plx_value': [190.,230.,375.],
        'plx_ref': ['ref_plx','ref_plx','ref_plx2'],
        'plx_qual': ['A','A','C'],
        'membership': membership,
        'parent_oid': [4,None,4],
        'h_link_ref': ['ref_hlink','','ref_hlink'],
        'otypes': ['*|**','BD?','*|**'],
        'ids': ['star1|id1','star2|id2','star3|id3'],
        'mag_i_value': [None,None,7.],
        'mag_j_value': [8.,None,6.],
        'mag_k_value': [9.,None,8.]
    })
    parents_without_plx = t0[:0].copy()
    parents_without_plx.add_row(['system1',None,None,None,None,None,4,
                                 '','','','','',None,None,'','',np.ma.masked,None,
                                 '','**','system1',None,None,None])
    children_without_plx = t0[:0].copy()
    children_without_plx.add_row(['planet1',None,None,None,None,None,5,
                                 '','','','','',None,None,'','',np.ma.masked,3,
                                 'ref_hlink2','Pl','planet1',None,None,None])

    ident = Table({
        'main_id': ['star1','star1','star3','star3','system1','planet1'],
        'oid': [1,1,3,3,4,5],
        'id': ['star1','id1','star3','id3','system1','planet1'],
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
    assert len(result_stars) == 3# stars and systems but removed planets

def test_create_ident_table(query_returns):
    # Execute
    result_sim_helptab, result_sim, t0, return_ident = query_returns

    # Verify
    ref = result_sim['provider']['provider_bibcode'][0]
    expected_ident = Table({
        'main_id': ['star1', 'star1','star3', 'star3','system1', 'planet1'],
        'id_ref': [ref, ref, ref, ref, ref, ref],
        'id': ['id1', 'star1','id3', 'star3','system1', 'planet1'],
    })

    assert len(setdiff(return_ident, expected_ident)) == 0

def test_create_h_link_table(query_returns,monkeypatch):
    # Data
    result_sim_helptab, result_sim, t0, return_ident = query_returns
    result_stars = creating_helpertable_stars(result_sim_helptab, result_sim)
    main_id_table = Table({
        'main_id': ['star1','star3', 'planet1'],
        'parent_main_id': ['system1','system1', 'star3'],
        'h_link_ref': ['ref_hlink','ref_hlink','ref_hlink2'],
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
        'main_id': ['star1','star3', 'planet1'],
        'parent_main_id': ['system1','system1', 'star3'],
        'h_link_ref': ['ref_hlink','ref_hlink','ref_hlink2'],
        'membership': [999999, 100,999999]
    })

    # Verify
    assert len(setdiff(result_h_link, expected_h_link)) == 0

    # this stuff below here is for further testing functions

    minimal_simbad = empty_dict.copy()

    minimal_simbad['objects'] = Table({
        'type': ['st'],
        'ids': ['id1|star1'],
        'main_id': ['star1']
    })
    minimal_simbad['sources'] = Table({
        'ref': ['ref1', 'ref2'],
        'provider_name': ['SIMBAD', 'SIMBAD']
    })

