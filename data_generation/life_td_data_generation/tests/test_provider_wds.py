from provider.wds import *
from astropy.table import Table
from provider.utils import query


def test_wds_queryable():
    wds_provider = Table()
    wds_provider['provider_url'] = ["http://tapvizier.u-strasbg.fr/TAPVizieR/tap"]
    adql_query = ["""SELECT top 10 *
                FROM "B/wds/wds" """]
    wds = query(wds_provider['provider_url'][0], adql_query[0])
    assert type(wds) == type(Table())


def test_wds_columns_queryable():
    wds_provider = Table()
    wds_provider['provider_url'] = ["http://tapvizier.u-strasbg.fr/TAPVizieR/tap"]
    adql_query = ["""SELECT TOP 10
                    WDS as wds_name, Comp as wds_comp,
                    sep1 as wds_sep1, sep2 as wds_sep2, 
                    Obs1 as wds_obs1, Obs2 as wds_obs2
                    FROM "B/wds/wds" """]
    wds = query(wds_provider['provider_url'][0], adql_query[0])
    assert wds.colnames == ['wds_name', 'wds_comp', 'wds_sep1', 'wds_sep2', 'wds_obs1', 'wds_obs2']
