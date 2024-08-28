# from building import *
# from astropy.table import Table, MaskedColumn

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

