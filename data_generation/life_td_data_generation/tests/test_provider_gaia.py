from provider.gaia import *
from astropy.table import Table, MaskedColumn


def test_create_mes_teff_st_table():
    # data
    gaia_ref='2022arXiv220800211G'
    phot = MaskedColumn(data = [4269.0712890625,4000.154541015625,4235.400390625,4235.4004],
                         name = 'teff_gspphot', mask=[True,False,False,False])#all but last made up)
    spec = MaskedColumn(data=[3060, 4250, 4000,4108.0],
                        name='teff_gspspec', mask=[False, True, False, False])
    gaia_helptab = Table(data=[
        ['3 Cnc', 'TIC 72090501', '6 Lyn','HD 166348'],
        spec,
        ['00000000000000000000000000000990099999995',
         '00000000000000000000000009900001099999995',
         '09009090000090090900900009900991099129995',
         '00100131000019999999999999999999999999999'],
        phot,
        [gaia_ref,gaia_ref,gaia_ref,gaia_ref]],
        names=['main_id', 'teff_gspspec', 'flags_gspspec','teff_gspphot', 'ref'],
        dtype=[object, float, object, float,object ])

    # function
    gaia = create_mes_teff_st_table(gaia_helptab)

    # assert
    assert len(gaia) == 6 #because two masked objects
    assert gaia['teff_st_ref'][np.where(gaia['teff_st_value']==4108.0)]== '2022arXiv220800211G GSP-Spec'
    assert gaia['teff_st_qual'][np.where(gaia['teff_st_value']==4108.0)] == 'D'
    assert gaia.colnames == ['main_id', 'teff_st_value', 'teff_st_qual',  'teff_st_ref']

