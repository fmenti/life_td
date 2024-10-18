from provider.sdb import *
from astropy.table import Table

def test_create_gk_helpertable():
    distance_cut_in_pc = 5.
    gk_disks = create_gk_helpertable(distance_cut_in_pc,'priv. comm.')
    columns = ['main_id', 'id', 'sdbid', 'name', '_r', 'raj2000', 'dej2000', 'pmra', 'pmde', 'gk_host_main_id', 'sp_type', 'sp_bibcode', 'plx_value', 'plx_err', 'plx_bibcode', 'type_short', 'type_long', 'star_comp_no', 'teff', 'e_teff', 'plx_arcsec', 'e_plx_arcsec', 'lstar', 'e_lstar', 'rstar', 'e_rstar', 'logg', 'e_logg', 'mh', 'e_mh', 'a_v', 'e_a_v', 'lstar_1pc', 'e_lstar_1pc', 'disk_r_comp_no', 'temp', 'e_temp', 'ldisk_lstar', 'e_ldisk_lstar', 'rdisk_bb', 'e_rdisk_bb', 'ldisk_1pc', 'e_ldisk_1pc', 'lam0', 'e_lam0', 'beta', 'e_beta', 'dmin', 'e_dmin', 'q', 'e_q', 'model_comps', 'parameters', 'evidence', 'chisq', 'dof', 'model_mtime', 'url', 'phot_url', 'mnest_url', 'model_url', 'type', 'disks_ref']
    plx_in_mas_cut=1000./distance_cut_in_pc
    
    assert gk_disks.colnames==columns
    assert len(gk_disks) > 0
    assert gk_disks['type'][0] == 'di'
    assert gk_disks['disks_ref'][0] == 'priv. comm.'
    assert max(gk_disks['plx_value']) > plx_in_mas_cut

def test_provider_gk():
    distance_cut_in_pc=5.
    
    gk=provider_gk(distance_cut_in_pc)
    
    ref=np.array(['priv. comm.'])
    provider_name=np.array(['Grant Kennedy Disks'])
    example_table=Table((ref,provider_name),names=('ref','provider_name' ))
        
    assert gk['sources']==example_table

