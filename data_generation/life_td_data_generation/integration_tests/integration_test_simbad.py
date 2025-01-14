from provider.simbad import *
from provider.utils import query

def test_simbad_queryable():
    simbad_provider=Table()
    simbad_provider['provider_url']=["http://simbad.u-strasbg.fr:80/simbad/sim-tap"]
    adql_query=["""SELECT top 10 *
                FROM basic"""]
    simbad=query(simbad_provider['provider_url'][0],adql_query[0])
    assert type(simbad)==type(Table())


