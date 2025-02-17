from provider.exo import *
from utils.io import load
from astropy.table import setdiff
from provider.utils import query

def test_exo_queryable():
    exo_provider=Table()
    exo_provider['provider_url']=["http://archives.ia2.inaf.it/vo/tap/projects"]
    adql_query=["""SELECT top 10 *
                FROM exomercat.exomercat """]
    exo=query(exo_provider['provider_url'][0],adql_query[0])
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

def bestmass(masstable):
    for i,flag in enumerate(masstable['mass_pl_sini_flag']):
        #if bestmass is mass and not msini measurement
        if masstable['bestmass_provenance'][i] == 'Mass':
            #if measurement is sini measurement
            if flag == 'True':
                masstable['mass_pl_qual'][i]=lower_quality(masstable['mass_pl_qual'][i])
        elif masstable['bestmass_provenance'][i] == 'Msini':
            if flag == 'False':
                masstable['mass_pl_qual'][i]=lower_quality(masstable['mass_pl_qual'][i])
    #is it possible, that I arrive at a different quality thingy than exomercat because I already assign some quality? this step is only needed, if both would end up with same quality
    return masstable
    
# def test_bestmass():    
#     # data
#     # load  exo_mes_mass_pl
#     [exo_mes_mass_pl] = load(['exo_mes_mass_pl'])
#     # group by planet_main_id
#     grouped_mass_pl = exo_mes_mass_pl.group_by('main_id')
#     ind=grouped_mass_pl.groups.indices
#     cols=['main_id','mass_pl_value','mass_pl_qual','mass_pl_sini_flag','bestmass_provenance']
#     for i in range(10):#range(len(ind)-1):
#         for j in range(ind[i],ind[i+1]):
#             print(grouped_mass_pl[cols][j])
#     #there are some with same qual now. those with siniflag true don't have quals yet
    
    
#     # assert
#     # exomercat bestmass has better qual in exo_mes_mass_pl for same object
#     assert False




