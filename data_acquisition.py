import pyvo as vo #catalog query
import astropy as ap #table handling

def save(cat,path):
    """
    Saves cat in the location path as .xml file.
    :param cat: astropy table to be saved
    :param path: path to where to save the cat table
    """
    ap.io.votable.writeto(ap.io.votable.from_table(cat), f'{path}.xml')
    return

def idref(cat,oldcolnames,newcolnames):
    """
    For each entry in one of the oldcolnames columns in table cat,
    looks up the corresponding source identifier source_id in the 
    source table and adds it in the column <newcolname>_source_idref. Deletes 
    the columns of oldcolnames.
    :param cat: astropy table
    :param oldcolnames: list of strings naming the columns in cat where the reference entries are
    :param newcolnames: list containing the strings to be added at the front of the 
        colmnname _source_idref.
    :return cat: astropy table
    """

#go over all reference parameter columns
    for k in range(len(oldcolnames)):
        colname=oldcolnames[k]
        param=newcolnames[k]
        #initiate empty columns
        cat[f'{param}_source_idref']=[0 for i in range(len(cat))]
        #go through all param bibcode entries in stars
        for i_cat in range(len(cat)):
            #go through all ref entries in sources
            for i_sources in range(len(sources)):
                if sources['ref'][i_sources]==cat[colname][i_cat]:
                    cat[f'{param}_source_idref'][i_cat]=sources['source_id'][i_sources]
        cat.remove_column(oldcolnames[k])
    return cat

#defining quantities
adql_query=[
"""
SELECT TOP 10 main_id,ra,dec,oid, coo_bibcode, plx_bibcode
FROM basic
WHERE basic.plx_value >=50.
""",
"""
SELECT TOP 20 id, name, bestmass, bestmass_url
FROM exomercat.exomercat
"""]
TAP_service=["http://simbad.u-strasbg.fr:80/simbad/sim-tap",
            "http://archives.ia2.inaf.it/vo/tap/projects"]
catalog_names=['stars','planets']
ref_columns=[['coo_bibcode','plx_bibcode'],['bestmass_url']]
new_columns=[['coord','plx'],['bestmass']]
provider_name=['SIMBAD','Exo-MerCat']
provider_bibcode=['2000A&AS..143....9W','2020A&C....3100370A']

#initializing tables
catalog = [ap.table.Table() for i in range(len(catalog_names))]
catalog_sources=[ap.table.Table() for i in range(len(catalog_names))]

#do the same for the stars and the planets objects
for i in range(len(catalog_names)):
    #Performs a query via TAP on the service. 
    service = vo.dal.TAPService(TAP_service[i])
    result=service.run_async(adql_query[i].format(**locals()), maxrec=160000)
    catalog[i]=result.to_table()
    #creates a source table from the references given in the ref_columns
    catalog_reflist=[]
    for j in ref_columns[i]:
        catalog_reflist.extend(catalog[i][j])
    catalog_sources[i]['ref']=catalog_reflist
    catalog_sources[i]=ap.table.unique(catalog_sources[i])#keeps only unique values
    #attaches service information
    catalog_sources[i]['provider_name']=[provider_name[i]]*len(catalog_sources[i])
    catalog_sources[i]['provider_url']=[TAP_service[i]]*len(catalog_sources[i])
    catalog_sources[i]['provider_bibcode']=[provider_bibcode[i]]*len(catalog_sources[i])
sources=ap.table.vstack([catalog_sources[0],catalog_sources[1]])
sources['source_id']=[i+1 for i in range(len(sources))]#introduces an identifier column
save(sources,'data/sources')
#do the same for the stars and the planets tables
for i in range(len(catalog_names)):
    catalog[i]=idref(catalog[i][:],ref_columns[i],new_columns[i])
    save(catalog[i],'data/'+catalog_names[i])