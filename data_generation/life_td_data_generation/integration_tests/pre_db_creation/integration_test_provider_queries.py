from astropy.table import Table
from provider.utils import query


def test_wds_queryable():
    wds_provider = Table()
    wds_provider["provider_url"] = [
        "http://tapvizier.u-strasbg.fr/TAPVizieR/tap"
    ]
    adql_query = [
        """SELECT top 10 *
                FROM "B/wds/wds" """
    ]
    wds = query(wds_provider["provider_url"][0], adql_query[0])
    assert type(wds) == type(Table())


def test_simbad_queryable():
    simbad_provider = Table()
    simbad_provider["provider_url"] = [
        "http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    ]
    adql_query = [
        """SELECT top 10 *
                FROM basic"""
    ]
    simbad = query(simbad_provider["provider_url"][0], adql_query[0])
    assert type(simbad) == type(Table())


def test_gaia_queryable():
    gaia_provider = Table()
    gaia_provider["provider_url"] = ["https://gea.esac.esa.int/tap-server/tap"]
    adql_query = [
        """SELECT top 10 *
                FROM gaiadr3.gaia_source """
    ]
    gaia = query(gaia_provider["provider_url"][0], adql_query[0])
    assert type(gaia) == type(Table())


def test_gaia_columns_queryable():
    gaia_provider = Table()
    gaia_provider["provider_url"] = ["https://gea.esac.esa.int/tap-server/tap"]
    adql_query = [
        """SELECT TOP 10
                    s.source_id ,p.mass_flame, p.radius_flame,
                    p.teff_gspphot, p.teff_gspspec, m.nss_solution_type, p.age_flame
                    FROM gaiadr3.gaia_source as s
                    JOIN gaiadr3.astrophysical_parameters as p ON s.source_id=p.source_id
                    LEFT JOIN gaiadr3.nss_two_body_orbit as m ON s.source_id=m.source_id"""
    ]
    gaia = query(gaia_provider["provider_url"][0], adql_query[0])
    assert gaia.colnames == [
        "source_id",
        "mass_flame",
        "radius_flame",
        "teff_gspphot",
        "teff_gspspec",
        "nss_solution_type",
        "age_flame",
    ]


def test_exo_queryable():
    exo_provider = Table()
    exo_provider["provider_url"] = [
        "http://archives.ia2.inaf.it/vo/tap/projects"
    ]
    adql_query = [
        """SELECT top 10 *
                FROM exomercat.exomercat """
    ]
    exo = query(exo_provider["provider_url"][0], adql_query[0])
    assert type(exo) == type(Table())


def test_exo_columns_queryable():
    exo_provider = Table()
    exo_provider["provider_url"] = [
        "http://archives.ia2.inaf.it/vo/tap/projects"
    ]
    adql_query = [
        """SELECT TOP 10
                    main_id, host, exomercat_name,
                    binary,letter
                    FROM exomercat.exomercat"""
    ]
    exo = query(exo_provider["provider_url"][0], adql_query[0])
    assert exo.colnames == [
        "main_id",
        "host",
        "exomercat_name",
        "binary",
        "letter",
    ]


def test_whole_exo_queryable():
    exo_provider = Table()
    exo_provider["provider_url"] = [
        "http://archives.ia2.inaf.it/vo/tap/projects"
    ]
    adql_query = [
        """SELECT *
                FROM exomercat.exomercat """
    ]
    exo = query(exo_provider["provider_url"][0], adql_query[0])
    assert type(exo) == type(Table())
