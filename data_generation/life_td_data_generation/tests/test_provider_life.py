from provider.life import *
from astropy.table import Table

def test_assign_lum_type():
    sptype=np.array(['M5.0IV','G2VI'])
    main_id=np.array(['G 227-48B','* mu. Dra C'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #no lum class given main sequence gets V assigned
    assert result['class_lum'][0]=='IV'
    assert result['class_lum'][1]=='VI'

def test_no_lum_class_in_sptype_string_but_V_assumed():
    sptype=np.array(['M3.51','M3:'])
    main_id=np.array(['G 227-48B','* mu. Dra C'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #no lum class given main sequence gets V assigned
    assert result['class_lum'][0]=='V'
    assert result['class_lum'][1]=='V'

def test_deal_with_minus_in_sptype():
    sptype=np.array(['K2-IIIbCa-1','K2-VbCa-1','K2-III','K-2V','K-2-V'])
    main_id=np.array(['* alf Ari','* alf AriV','* alf Arishort','test','test2'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    assert result['class_lum'][0]=='III'
    assert result['class_lum'][1]=='V'
    assert result['class_lum'][2]=='III'
    
    assert result['class_temp'][3]=='K'
    assert result['class_lum'][3]=='V'
    assert result['class_temp_nr'][3]=='2'
    assert result['class_temp'][4]=='K'
    assert result['class_lum'][4]=='V'
    assert result['class_temp_nr'][4]=='2'

    
def test_deal_with_leading_d_sptype():
    sptype=np.array(['dM3.51','dM3:','dM5.0'])
    main_id=np.array(['G 227-48B','* mu. Dra C','FBS 1415+456'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #no lum class given main sequence gets V assigned
    assert result['class_lum'][0]=='V'
    assert result['class_lum'][1]=='V'
    #dwarf d gets V assigned
    assert result['class_lum'][2]=='V'

    
def test_assign_null_values():
    sptype=np.array(['','D2.3V'])
    main_id=np.array(['test','test2'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #empty entry gets '?' assigned
    assert result['class_temp'][0]=='?'
    assert result['class_temp_nr'][0]=='?'
    assert result['class_lum'][0]=='?'
    assert result['class_ref'][0]=='?'
    
    #non main sequence
    assert result['class_temp'][1]=='?'
    assert result['class_temp_nr'][1]=='?'
    assert result['class_lum'][1]=='?'
    assert result['class_ref'][1]=='?'
    
    