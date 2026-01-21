import importlib #reloading external functions after modification

#self created modules
import life_td as ltd
importlib.reload(ltd)#reload module after changing it





def test_warnings():

    provider_tables_dict, database_tables= ltd.create_life_td(5)

    def metachecking_prov(prov_tables_dict):
        for prov_name, prov_tables in prov_tables_dict.items():
            print(f"testing {prov_name}...")
            for table_name, table in prov_tables.items():
                print(f"testing {table_name}...")
                assert table.meta == {}

    def metachecking_cat(cat):
        for table_name, table in cat.items():
            print(f"testing {table_name}...")
            assert table.meta == {}

    metachecking_prov(provider_tables_dict)
    metachecking_cat(database_tables)
