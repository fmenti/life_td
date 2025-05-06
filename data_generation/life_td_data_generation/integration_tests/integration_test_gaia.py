from provider.gaia import *
from provider.utils import query


def test_gaia_queryable():
    gaia_provider=Table()
    gaia_provider['provider_url']=["https://gea.esac.esa.int/tap-server/tap"]
    adql_query=["""SELECT top 10 *
                FROM gaiadr3.gaia_source """]
    gaia=query(gaia_provider['provider_url'][0],adql_query[0])
    assert type(gaia)==type(Table())

def test_gaia_columns_queryable():  
    gaia_provider=Table()
    gaia_provider['provider_url']=["https://gea.esac.esa.int/tap-server/tap"]
    adql_query=["""SELECT TOP 10
                    s.source_id ,p.mass_flame, p.radius_flame,
                    p.teff_gspphot, p.teff_gspspec, m.nss_solution_type, p.age_flame
                    FROM gaiadr3.gaia_source as s
                    JOIN gaiadr3.astrophysical_parameters as p ON s.source_id=p.source_id
                    LEFT JOIN gaiadr3.nss_two_body_orbit as m ON s.source_id=m.source_id"""]
    gaia=query(gaia_provider['provider_url'][0],adql_query[0])
    assert gaia.colnames==['source_id', 'mass_flame', 'radius_flame', 'teff_gspphot',
                           'teff_gspspec', 'nss_solution_type', 'age_flame']

