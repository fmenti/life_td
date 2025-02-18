from provider.gaia import *
from astropy.table import Table

def test_assign_quality():
    gaia_mes_teff_st_spec = Table(data=[
                       ['3 Cnc', 'TIC 72090501', '6 Lyn'],
                     [3060,4250,4000],
                     ['00000000000000000000000000000990099999995',
                      '00000000000000000000000009900001099999995',
                      '09009090000090090900900009900991099129995']],
                 names=['main_id','teff_gspspec','flags_gspspec'],
                 dtype=[object, float,object])
    gaia_mes_teff_st_spec=assign_quality(gaia_mes_teff_st_spec)

    assert gaia_mes_teff_st_spec['teff_st_qual'][0]=='B'
    assert gaia_mes_teff_st_spec['teff_st_qual'][2] == 'C'
