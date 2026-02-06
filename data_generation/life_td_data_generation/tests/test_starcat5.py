from catalog.starcat5 import (
    choose_service,
)

def test_choose_service():
    assert choose_service("")=="http://localhost:8080/tap"
    assert choose_service("heid") == "http://dc.zah.uni-heidelberg.de/tap"
    assert choose_service("gvo") == "http://dc.g-vo.org/tap"





