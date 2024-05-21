import pyvo as vo #catalog query

adql_query="""
    SELECT TOP 10 object_id, main_id FROM life_td.object
    WHERE type='st'"""

def example_query(adql_query=adql_query):
    """
    Performs an example query via TAP on the LIFE Target Database.
    
    :param adql_query: Query, default is main identifier of the first 
        ten stellar objects.
    :type adql_query: str
    :returns: Result of the query.
    :rtype: astropy.table.table.Table   
    """
    
    def query(link,query):
        """
        Performs a query via TAP on the service given.
        
        :param link: Service access URL.
        :type link: str
        :param query: Query to be asked of the external database service
             in ADQL.
        :type query: str
        :returns: Result of the query.
        :rtype: astropy.table.table.Table
        """
        
        #defining the vo service using the given link
        service = vo.dal.TAPService(link)

        result=service.run_async(query.format(**locals()), 
                                                maxrec=160000)
    
        cat=result.to_table()
        return cat
    
    service='http://dc.zah.uni-heidelberg.de/tap'
    cat=query(service,adql_query)
    return cat
