from provider.simbad import *
from sdata import table_names




def test_provider_simbad():
    # data
    distance_cut_in_pc = 5
    empty_tables =  ['best_h_link','planet_basic','disk_basic',
                     'mes_mass_pl','mes_teff_st','mes_radius_st',
                     'mes_mass_st','mes_sep_ang']

    # function
    sim = provider_simbad(distance_cut_in_pc)

    # assert
    for table_name in table_names:
        assert type(sim[table_name]) == type(Table())
        if table_name in empty_tables:
            assert len(sim[table_name]) == 0
        else:
            assert len(sim[table_name]) > 0
