from provider.wds import *
from provider.utils import query

def test_wds_queryable():
    wds_provider=Table()
    wds_provider['provider_url']=["http://tapvizier.u-strasbg.fr/TAPVizieR/tap"]
    adql_query=["""SELECT top 10 *
                FROM "B/wds/wds" """]
    wds=query(wds_provider['provider_url'][0],adql_query[0])
    assert type(wds)==type(Table())


