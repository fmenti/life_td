from provider.simbad import *
from astropy.table import Table
from provider.utils import query


def test_sim_queryable():
    sim_provider=Table()
    sim_provider['provider_url']=[
            "http://simbad.u-strasbg.fr:80/simbad/sim-tap"]
    adql_query=["""SELECT top 10 *
                FROM basic """]
    sim=query(sim_provider['provider_url'][0],adql_query[0])
    assert type(sim)==type(Table())
