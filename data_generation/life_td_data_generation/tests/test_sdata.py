from sdata import paras_dict, para_columns_dict


def test_paras_dict():
    p_d = paras_dict

    assert p_d["sources"] == []


def test_best_parameters_ingestion_loop():
    parameter_groups = para_columns_dict["star_basic"]

    plx_colnames = ["plx_value", "plx_err", "plx_qual",
                    "plx_source_idref",
                    "plx_ref"]

    for para in parameter_groups:
        print(para)
        print(parameter_groups[para]["columns"])
        if para == "plx":
            for i in range(len(parameter_groups[para]["columns"])):
                assert parameter_groups[para]["columns"][i] == plx_colnames[i]

    assert False
