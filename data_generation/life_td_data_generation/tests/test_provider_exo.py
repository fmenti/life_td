from provider.exo import *
from astropy.table import Table
from provider.utils import query


def test_exo_queryable():
    exo_provider=Table()
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    adql_query=["""SELECT top 10 *
                FROM exomercat.exomercat """]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
    print(exo.colnames)
    assert type(exo)==type(Table())

def test_exo_columns_queryable():  
    exo_provider=Table()
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    adql_query=["""SELECT TOP 10
                    main_id, host,
                    binary,letter
                    FROM exomercat.exomercat """]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
    assert exo.colnames==['main_id', 'host', 'binary', 'letter']