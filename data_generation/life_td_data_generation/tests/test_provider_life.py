import numpy as np
from astropy.table import Table
from provider.life import (
    assign_null_values,
    deal_with_leading_d_sptype,
    extract_lum_class,
    initiate_columns,
    modeled_param,
    spec,
    sptype_string_to_class,
)


def test_extract_lum_class():
    sptype = ["M5.0IV", "G2VI"]

    assert extract_lum_class(4, sptype[0]) == "IV"
    assert extract_lum_class(2, sptype[1]) == "VI"


def test_assign_null_values():
    sptype = np.array(["", "D2.3V"])
    main_id = np.array(["test", "test2"])
    table = Table(
        [main_id, sptype],
        names=["main_id", "sptype_string"],
        dtype=[object, object],
    )
    example_table = initiate_columns(
        table,
        ["class_temp", "class_temp_nr", "class_lum", "class_ref"],
        [object, object, object, object],
        [True, True, True, True],
    )

    result = assign_null_values(example_table, 0)

    # empty entry gets '?' assigned
    assert result["class_temp"][0] == "?"
    assert result["class_temp_nr"][0] == "?"
    assert result["class_lum"][0] == "?"
    assert result["class_ref"][0] == "?"


def test_deal_with_leading_d_sptype():
    sptype = np.array(["M3.51", "dM3:", "dM5.0"])
    main_id = np.array(["G 227-48B", "* mu. Dra C", "FBS 1415+456"])
    table = Table(
        (main_id, sptype),
        names=("main_id", "sptype_string"),
        dtype=[object, object],
    )
    example_table = initiate_columns(
        table,
        ["class_temp", "class_temp_nr", "class_lum", "class_ref"],
        [object, object, object, object],
        [True, True, True, True],
    )

    # no lum class given main sequence gets V assigned
    assert deal_with_leading_d_sptype(example_table, 0) == "M3.51"
    assert example_table["class_lum"][0] != "V"  # because no d

    assert deal_with_leading_d_sptype(example_table, 1) == "M3:"
    assert example_table["class_lum"][1] == "V"

    assert deal_with_leading_d_sptype(example_table, 2) == "M5.0"
    assert example_table["class_lum"][2] == "V"


def test_assign_lum_type():
    sptype = np.array(["M5.0IV", "G2VI"])
    main_id = np.array(["G 227-48B", "* mu. Dra C"])
    example_table = Table(
        (main_id, sptype),
        names=("main_id", "sptype_string"),
        dtype=[object, object],
    )

    result = sptype_string_to_class(example_table, "fake_ref")

    # no lum class given main sequence gets V assigned
    assert result["class_lum"][0] == "IV"
    assert result["class_lum"][1] == "VI"


def test_no_lum_class_in_sptype_string_but_v_assumed():
    sptype = np.array(["M3.51", "M3:"])
    main_id = np.array(["G 227-48B", "* mu. Dra C"])
    example_table = Table(
        (main_id, sptype),
        names=("main_id", "sptype_string"),
        dtype=[object, object],
    )

    result = sptype_string_to_class(example_table, "fake_ref")

    # no lum class given main sequence gets V assigned
    assert result["class_lum"][0] == "V"
    assert result["class_lum"][1] == "V"


def test_deal_with_minus_in_sptype():
    sptype = np.array(["K2-IIIbCa-1", "K2-VbCa-1", "K2-III", "K-2V", "K-2-V"])
    main_id = np.array(
        ["* alf Ari", "* alf AriV", "* alf Arishort", "test", "test2"]
    )
    example_table = Table(
        (main_id, sptype),
        names=("main_id", "sptype_string"),
        dtype=[object, object],
    )

    result = sptype_string_to_class(example_table, "fake_ref")

    assert result["class_lum"][0] == "III"
    assert result["class_lum"][1] == "V"
    assert result["class_lum"][2] == "III"

    assert result["class_temp"][3] == "K"
    assert result["class_lum"][3] == "V"
    assert result["class_temp_nr"][3] == "2"
    assert result["class_temp"][4] == "K"
    assert result["class_lum"][4] == "V"
    assert result["class_temp_nr"][4] == "2"


def test_modeled_param():
    mp = modeled_param()  # create model table as votable

    assert (
        mp["Teff"][np.where(mp["SpT"] == "M5V")] == 3060
    )  # issue M5V not M5.0V
    assert mp["Teff"][np.where(mp["SpT"] == "M3.5V")] == 3270


def test_assign_teff():
    temp = np.array(['M', 'M', 'M', 'M', 'M','K'])
    temp_nr = np.array(['5.0', '3.5', '5.5', '6.5', '4.0','2.5'])
    lum = np.array(['V', 'V', 'V', 'V', 'V','V'])
    main_id = np.array(['test', 'test2', 'test3', 'test4', 'test5','test6'])
    cat = Table((main_id, temp, temp_nr, lum),
                names=('main_id', 'class_temp', 'class_temp_nr', 'class_lum'),
                dtype=[object, object, object, object])

    result = spec(cat)

    # assert
    assert result['mod_Teff'][
               np.where(result['class_temp_nr'] == '2.5')] == 5100 # should use the one from 2.0
    assert result['mod_Teff'][np.where(result['class_temp_nr'] == '3.5')] == 3270
    assert result['mod_Teff'][np.where(result['class_temp_nr'] == '4.0')] == 3210
    assert result['mod_Teff'][np.where(result['class_temp_nr'] == '5.0')] == 3060
    assert result['mod_Teff'][np.where(result['class_temp_nr'] == '5.5')] == 2930
    assert result['mod_Teff'][np.where(result['class_temp_nr'] == '6.5')] == 2740


def test_spec():
    main_id = np.array(["test", "test2", "test3", "test4", "test5"])
    temp = np.array(["M", "M", "M", "M", "M"])
    temp_nr = np.array(["0", "3.5", "5.5", "6.5", "4"])
    lum = np.array(["V", "V", "V", "V", "V"])

    cat = Table(
        (main_id, temp, temp_nr, lum),
        names=("main_id", "class_temp", "class_temp_nr", "class_lum"),
        dtype=[object, object, object, object],
    )
    # somehow within 30 pc I have one element that has int in there
    cat["class_temp_nr"][np.where(cat["class_temp_nr"] == "0")] = int("0")
    result = spec(cat)

    assert (
        result["mod_Teff"][np.where(result["class_temp_nr"] == "6.5")] == 2740
    )
