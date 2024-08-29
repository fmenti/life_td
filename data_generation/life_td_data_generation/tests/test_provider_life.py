from provider.life import *
from astropy.table import Table


def test_no_lum_class_in_sptype_string_but_V_assumed():
    sptype=np.array(['M3.51','M3:','dM5.0','M5.0IV','','D2.3V','G2VI'])
    main_id=np.array(['G 227-48B','* mu. Dra C','FBS 1415+456','test','test2','test3','test4'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #no lum class given main sequence gets V assigned
    assert result['class_lum'][0]=='V'
    assert result['class_lum'][1]=='V'

    
def test_deal_with_leading_d_sptype():
    sptype=np.array(['M3.51','M3:','dM5.0','M5.0IV','','D2.3V','G2VI'])
    main_id=np.array(['G 227-48B','* mu. Dra C','FBS 1415+456','test','test2','test3','test4'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #no lum class given main sequence gets V assigned
    assert result['class_lum'][0]=='V'
    assert result['class_lum'][1]=='V'
    #dwarf d gets V assigned
    assert result['class_lum'][2]=='V'
    
    #empty entry gets '?' assigned
    assert result['class_lum'][4]=='?'
    #non main sequence V gets '?' assigned
    assert result['class_lum'][5]=='?'
    
def test_assign_null_values():
    sptype=np.array(['M3.51','M3:','dM5.0','M5.0IV','','D2.3V','G2VI'])
    main_id=np.array(['G 227-48B','* mu. Dra C','FBS 1415+456','test','test2','test3','test4'])
    example_table=Table((main_id,sptype),names=('main_id','sptype_string' ))
    
    result=sptype_string_to_class(example_table,'fake_ref')
    
    #empty entry gets '?' assigned
    assert result['class_temp'][4]=='?'
    assert result['class_temp_nr'][4]=='?'
    assert result['class_lum'][4]=='?'
    assert result['class_ref'][4]=='?'
    
    #non main sequence
    assert result['class_temp'][5]=='?'
    assert result['class_temp_nr'][5]=='?'
    assert result['class_lum'][5]=='?'
    assert result['class_ref'][5]=='?'
    
    