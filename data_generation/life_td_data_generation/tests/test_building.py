from building import *
from astropy.table import Table, MaskedColumn

# def test_idsjoin_no_mask():
#     cat = Table(data=[['* 61 Cyg b|61 Cyg b', '', 'GCRV 13273|LTT 16180'],
#                       ['61 Cyg b|markmuster', 'HD 201091|LSPM J2106+3844N', 'USNO-B1.0 1287-00443364']],
#                 names=['ids1','ids2'], 
#                 dtype=[object, object])
#     cat=idsjoin(cat,'ids1','ids2')
    
    
#     for ID in cat['ids']:
#         assert ID.split('|').count('')==0
    
#     assert cat['ids'][0].split('|').count('61 Cyg b')==1
#     assert cat['ids'][0].split('|').count('* 61 Cyg b')==1
#     assert cat['ids'][0].split('|').count('markmuster')==1
    
# def test_idsjoin_with_mask():
#     a=MaskedColumn(data=['* 61 Cyg b|61 Cyg b', '', 'GCRV 13273|LTT 16180'],name='ids1',mask=[True,False,False])
#     cat = Table(data=[a,
#                       ['61 Cyg b|markmuster', 'HD 201091|LSPM J2106+3844N', 'USNO-B1.0 1287-00443364']],
#                 names=['ids1','ids2'], 
#                 dtype=[object, object])
#     cat=idsjoin(cat,'ids1','ids2')
    
    
#     for ID in cat['ids']:
#         assert ID.split('|').count('')==0
    
#     assert cat['ids'][0].split('|').count('61 Cyg b')==1
#     assert cat['ids'][0].split('|').count('markmuster')==1

def test_best_para_id():
    mes_table=Table(data=[[1,2,3],
                          ['* 61 Cyg b', '61 Cyg b', 'GCRV 13273'],
                       ['2000A&AS..143....9W', '1925AnHar.100...17C','2001AJ....122.3466M' ]],
                 names=['object_idref','id','id_source_idref'], 
                 dtype=[int,object,int])
    best_para_table=best_para_id(mes_table)
    
    assert False