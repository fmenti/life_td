from provider.sdb import *
from astropy.table import Table
from astropy import io

def test_create_sdb_helpertable():
    distance_cut_in_pc = 5.
    data=io.votable.parse_single_table(
        Path().additional_data+"sdb_30pc_09_02_2024.xml").to_table()
    sdb_disks = create_sdb_helpertable(distance_cut_in_pc,'priv. comm.',data)
    columns = ['main_id', 'id', 'sdbid', 'name', '_r', 'raj2000', 'dej2000', 'pmra', 'pmde', 'sdb_host_main_id', 'sp_type', 'sp_bibcode', 'plx_value', 'plx_err', 'plx_bibcode', 'type_short', 'type_long', 'star_comp_no', 'teff', 'e_teff', 'plx_arcsec', 'e_plx_arcsec', 'lstar', 'e_lstar', 'rstar', 'e_rstar', 'logg', 'e_logg', 'mh', 'e_mh', 'a_v', 'e_a_v', 'lstar_1pc', 'e_lstar_1pc', 'disk_r_comp_no', 'temp', 'e_temp', 'ldisk_lstar', 'e_ldisk_lstar', 'rdisk_bb', 'e_rdisk_bb', 'ldisk_1pc', 'e_ldisk_1pc', 'lam0', 'e_lam0', 'beta', 'e_beta', 'dmin', 'e_dmin', 'q', 'e_q', 'model_comps', 'parameters', 'evidence', 'chisq', 'dof', 'model_mtime', 'url', 'phot_url', 'mnest_url', 'model_url', 'type', 'disks_ref']
    plx_in_mas_cut=1000./distance_cut_in_pc
    
    assert sdb_disks.colnames==columns
    assert len(sdb_disks) > 0
    assert sdb_disks['type'][0] == 'di'
    assert sdb_disks['disks_ref'][0] == 'priv. comm.'
    assert max(sdb_disks['plx_value']) > plx_in_mas_cut

def test_provider_sdb():
    distance_cut_in_pc=5.
    data=io.votable.parse_single_table(
        Path().additional_data+"sdb_30pc_09_02_2024.xml").to_table()
    sdb=provider_sdb(distance_cut_in_pc,data)
    
    ref=np.array(['priv. comm.'])
    provider_name=np.array(['Grant Kennedy Disks'])
    example_table=Table((ref,provider_name),names=('ref','provider_name' ))
        
    assert sdb['sources']==example_table


