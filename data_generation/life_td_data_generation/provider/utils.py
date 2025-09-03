"""
Generates the data for the database for each of the data providers separately.
"""

from datetime import datetime

import numpy as np  # arrays
from astropy.table import (
    Column,
    MaskedColumn,
    Table,
    column,
    join,
    table,
    unique,
    vstack,
)
from pyvo.dal import TAPService

from utils.io import load


def initiate_columns(
    table_obj: Table, columns: list[str], types: list[type], mask: list[bool]
) -> Table:
    """
    Initialize multiple columns on a table with dtype and mask settings.

    :param table_obj: Table to receive the new columns.
    :type table_obj: astropy.table.Table
    :param columns: Column names to create.
    :type columns: list[str]
    :param types: Dtypes for the respective columns.
    :type types: list[type]
    :param mask: Mask flags; True -> MaskedColumn, False -> Column.
    :type mask: list[bool]
    :returns: The same table with the new columns added.
    :rtype: astropy.table.Table
    """
    for name, dtype, use_mask in zip(columns, types, mask):
        if use_mask:
            table_obj[name] = MaskedColumn(dtype=dtype, length=len(table_obj))
        else:
            table_obj[name] = Column(dtype=dtype, length=len(table_obj))
    return table_obj



def create_provider_table(
    provider_name: str,
    provider_url: str,
    provider_bibcode: str,
    provider_access: str = datetime.now().strftime("%Y-%m-%d"),
) -> Table:
    """
    Build a one-row provider metadata table.

    :param provider_name: Human-readable provider name.
    :type provider_name: str
    :param provider_url: TAP/base URL for the provider service.
    :type provider_url: str
    :param provider_bibcode: ADS bibcode or similar reference marker.
    :type provider_bibcode: str
    :param provider_access: Access date (YYYY-MM-DD). Defaults to today.
    :type provider_access: str
    :returns: Table with provider metadata columns.
    :rtype: astropy.table.Table
    """
    print(f"Trying to create {provider_name} tables from {provider_access}...")
    provider_table = Table()
    provider_table["provider_name"] = [provider_name]
    provider_table["provider_url"] = [provider_url]
    provider_table["provider_bibcode"] = [provider_bibcode]
    provider_table["provider_access"] = [provider_access]
    return provider_table


def query(
    link: str, adql_query: str, upload_tables: list[table.Table] = []
) -> table.Table:
    """
    Perform a TAP query against a service.

    If upload tables are provided, they are made available to the service
    as TAP_UPLOAD tables.

    :param link: Service access URL.
    :type link: str
    :param adql_query: Query to execute (ADQL).
    :type adql_query: str
    :param upload_tables: Optional tables to upload for join operations.
    :type upload_tables: list[astropy.table.table.Table]
    :returns: Result table returned by the TAP service.
    :rtype: astropy.table.table.Table
    """
    service = TAPService(link)
    if upload_tables == []:
        result = service.run_async(
            adql_query.format(**locals()), maxrec=1600000
        )
    else:
        tables = {}
        for i, t in enumerate(upload_tables, start=1):
            tables[f"t{i}"] = t
        result = service.run_async(
            adql_query, uploads=tables, timeout=None, maxrec=1600000
        )
    return result.to_table()


def remove_catalog_description(cat: table.Table, no_description: bool) -> table.Table:
    """
    Remove description metadata from all columns of a table.

    Useful when merging catalogs has left descriptions inconsistent.

    :param cat: Catalog to sanitize.
    :type cat: astropy.table.table.Table
    :param no_description: When True, clear all column descriptions.
    :type no_description: bool
    :returns: Catalog with cleared descriptions (if requested).
    :rtype: astropy.table.table.Table
    """
    for col in cat.colnames:
        if no_description:
            cat[col].description = ""
    return cat


def fill_sources_table(
    cat: table.Table,
    ref_columns: list[str],
    provider: str,
    old_sources: table.Table = Table(),
) -> table.Table:
    """
    Create or update the sources table from reference columns of a table.

    Entries are unique. Output has two columns: 'ref' and 'provider_name'.

    :param cat: Table on which to gather references.
    :type cat: astropy.table.table.Table
    :param ref_columns: Column names containing reference strings.
    :type ref_columns: list[str]
    :param provider: Provider name tag for all collected references.
    :type provider: str
    :param old_sources: Previously created sources table to extend.
    :type old_sources: astropy.table.table.Table
    :returns: Sources table containing unique refs and provider labels.
    :rtype: astropy.table.table.Table
    """
    if len(cat) > 0:
        cat_sources = Table()
        cat_reflist: list[object] = []

        for k in range(len(ref_columns)):
            # Skip masked entries if the column is a MaskedColumn.
            if type(cat[ref_columns[k]]) == column.MaskedColumn:
                cat_reflist.extend(
                    cat[ref_columns[k]][
                        np.where(cat[ref_columns[k]].mask == False)
                    ]
                )
            else:
                cat_reflist.extend(cat[ref_columns[k]])

        cat_sources["ref"] = cat_reflist
        cat_sources = unique(cat_sources)

        cat_sources["provider_name"] = [
            provider for _ in range(len(cat_sources))
        ]

        sources = vstack([old_sources, cat_sources])
        sources = unique(sources)
    else:
        sources = old_sources
    return sources


def create_sources_table(
    tables: list[table.Table], ref_columns: list[list[str]], provider_name: str
) -> table.Table:
    """
    Build a sources table from multiple input tables and reference columns.

    :param tables: List of input tables to scan for references.
    :type tables: list[astropy.table.table.Table]
    :param ref_columns: Parallel list of reference-column name lists.
    :type ref_columns: list[list[str]]
    :param provider_name: Provider label to attach to each reference.
    :type provider_name: str
    :returns: Final sources table with unique references.
    :rtype: astropy.table.table.Table
    """
    sources = Table()
    for cat, refs in zip(tables, ref_columns):
        sources = fill_sources_table(cat, refs, provider_name, sources)
    return sources


class OidCreator:
    """
    Create ADQL query builder for fetch_main_id using an 'oid' join column.
    """

    def __init__(self, name: str, colname: str) -> None:
        """
        :param name: Output alias for SIMBAD main_id (e.g. 'sim_main_id').
        :type name: str
        :param colname: Column name in the uploaded table to join on.
        :type colname: str
        """
        self.name = name
        self.colname = colname

    def create_main_id_query(self) -> str:
        """
        Build the ADQL query that joins SIMBAD basic on an 'oid' column.

        :returns: ADQL string with column alias and join.
        :rtype: str
        """
        return (
            "SELECT b.main_id AS "
            + self.name
            + """,t1.*
                FROM basic AS b
                JOIN TAP_UPLOAD.t1 ON b.oid=t1."""
            + self.colname
        )


class IdentifierCreator:
    """
    Create ADQL query builder for fetch_main_id using an identifier column.
    """

    def __init__(self, name: str, colname: str) -> None:
        """
        :param name: Output alias for SIMBAD main_id (e.g. 'sim_main_id').
        :type name: str
        :param colname: Identifier column name in the uploaded table.
        :type colname: str
        """
        self.name = name
        self.colname = colname

    def create_main_id_query(self) -> str:
        """
        Build the ADQL query that joins via SIMBAD 'ident' on identifiers.

        :returns: ADQL string with column alias and identifier join.
        :rtype: str
        """
        return (
            "SELECT b.main_id AS "
            + self.name
            + """,t1.*
                    FROM basic AS b
                    JOIN ident ON ident.oidref=b.oid
                        JOIN TAP_UPLOAD.t1 ON ident.id=t1."""
            + self.colname
        )


def fetch_main_id(
    cat: table.Table, id_creator=OidCreator(name="main_id", colname="oid")
) -> table.Table:
    """
    Joins main_id from simbad to the column colname.

    Returns the whole table cat but without any rows where no simbad
    main_id was found.

    :param cat: Astropy table containing column colname.
    :type cat: astropy.table.table.Table
    :param id_creator: OidCreator or IdentifierCreator object.
    :return: Table with all main SIMBAD identifiers that could be found
        in column "name".
    :rtype: astropy.table.table.Table
    """
    # improvement idea to be performed at one point
    # tbd option to match on position instead of main_id or oid
    # SIMBAD TAP service
    TAP_service = "http://simbad.u-strasbg.fr:80/simbad/sim-tap"
    # performing query using external function
    main_id_query = id_creator.create_main_id_query()
    cat = query(TAP_service, main_id_query, [cat])
    return cat


def distance_cut(cat: table.Table, colname: str, main_id: bool = True):
    """
    Sorts out objects not within the provider_simbad distance cut.

    :param cat: Astropy table to be matched against sim_objects table.
    :type cat: astropy.table.table.Table
    :param str colname: Name of the column to use for the match.
    :return: Table like cat without any objects not found in
        sim_objects.
    :rtype: astropy.table.table.Table
    """
    if main_id:
        [sim] = load(["sim_objects"])
        sim.rename_columns(["main_id", "ids"], ["temp1", "temp2"])
        cat = join(
            cat, sim["temp1", "temp2"], keys_left=colname, keys_right="temp1"
        )
        cat.remove_columns(["temp1", "temp2"])
    else:
        [sim] = load(["sim_ident"])
        sim.rename_columns(["id"], ["temp1"])
        cat = join(
            cat, sim["temp1", "main_id"], keys_left=colname, keys_right="temp1"
        )
        cat.remove_columns(["temp1"])
    return cat


def nullvalues(cat, colname, nullvalue, verbose=False):
    """
    This function fills masked entries of specified column.

    :param cat: Astropy table containing the column colname.
    :type cat: astropy.table.table.Table
    :param str colname: Name of a column.
    :param nullvalue: Value to be placed instead of masked elements.
    :type nullvalue: str or float or bool
    :param verbose: If True prints message if the column is not an
        astropy masked column, defaults to False.
    :type verbose: bool, optional
    :return: Astropy table with masked elements of colname replaced
        by nullvalue.
    :rtype: astropy.table.table.Table
    """
    if type(cat[colname]) == column.MaskedColumn:
        cat[colname].fill_value = nullvalue
        cat[colname] = cat[colname].filled()
    elif verbose:
        print(colname, "is no masked column")
    return cat


def replace_value(cat, column, value, replace_by):
    """
    This function replaces values.

    :param cat: Table containing column specified as colname.
    :type cat: astropy.table.table.Table
    :param str column: Designates column in which to replace the
        entries.
    :param value: Entry to be replaced.
    :type value: str or float or bool
    :param replace_by: Entry to be put in place of param value.
    :type replace_by: str or float or bool
    :return: Table with replaced entries.
    :rtype: astropy.table.table.Table
    """
    cat[column][np.where(cat[column] == value)] = [
        replace_by
        for i in range(len(cat[column][np.where(cat[column] == value)]))
    ]
    return cat


def ids_from_ident(ident, objects):
    """
    Concatenates identifier entries of same main_id object.

    This function extracts the identifiers of common main_id objects in
    column id of table ident using the delimiter | and stores the result
    in the column ids of the table objects.

    :param ident: Table containing the rows main_id and id
    :type ident: astropy.table.table.Table
    :param objects: Table containing the columns main_id and ids
    :type objects: astropy.table.table.Table
    :returns: Filled out table objects.
    :rtype: astropy.table.table.Table
    """
    grouped_ident = ident.group_by("main_id")
    ind = grouped_ident.groups.indices
    for i in range(len(ind) - 1):
        # -1 is needed because else ind[i+1] is out of bonds
        ids = []
        for j in range(ind[i], ind[i + 1]):
            ids.append(grouped_ident["id"][j])
        ids = "|".join(ids)
        objects.add_row([grouped_ident["main_id"][ind[i]], ids])
    return objects


def lower_quality(qual):
    if qual == "A":
        qual = "B"
    elif qual == "B":
        qual = "C"
    elif qual == "C":
        qual = "D"
    elif qual == "D":
        qual = "E"
    return qual
