from provider.exo2 import *
from utils.io import load
from astropy.table import setdiff
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
                    main_id, host, exomercat_name,
                    binary,letter
                    FROM exomercat.exomercat"""]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
    assert exo.colnames==['main_id', 'host', 'exomercat_name', 'binary', 'letter']

def test_whole_exo_queryable():
    exo_provider=Table()
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    adql_query=["""SELECT *
                FROM exomercat.exomercat """]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
    print(exo.colnames)
    assert type(exo)==type(Table())

def test_query_or_load_exomercat():
    exo, exo_helptab = query_or_load_exomercat()

    assert exo['provider']['provider_url'][0]=="http://archives.ia2.inaf.it/vo/tap/projects"
    assert len(exo_helptab)>0

    
def test_create_exo_helpertable():    
    exo, exo_helptab = create_exo_helpertable()
    [removed_objects] = load(['exomercat_removed_objects'])
    not_correctly_removed_objects = np.isin(exo_helptab['exomercat_name'],
                                           removed_objects['exomercat_name'])
    assert exo['provider']['provider_url'][0]=="http://archives.ia2.inaf.it/vo/tap/projects"
    assert len(exo_helptab['exomercat_name'][np.where(not_correctly_removed_objects)]) == 0


