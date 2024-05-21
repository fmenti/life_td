.. _usage:

Usage
=====

Database Access Tutorial
-----------------------------

.. _tutorial:

To access the data from the LIFE database you can either use the 
graphical user interface `TOPCAT <https://www.star.bris.ac.uk/~mbt/
topcat/>`_ or do it in a python environment.

TOPCAT:
Install TOPCAT. Lunch TOPCAT. Click the "VO" button on the top header. Click on "Table Access Protocol (TAP) Query". In the newly opened window choose from the dropdown menue of "All TAP services" the "GAVO DC TAP" one.

Python:
You need to install PyVO. Afterwards just copy paste the code below. It will give you back the main identifiers of the first ten stars in the database. Feel free to adapt the adql_query variable to your liking. If you are unsure of what else you could query have a look at the provided `examples <http://dc.zah.uni-heidelberg.de/life/q/ex/examples>`_.

.. code-block:: python

   import pyvo as vo
   
    adql_query="""
        SELECT TOP 10 object_id, main_id FROM life_td.object
        WHERE type='st'"""
        
    link='http://dc.zah.uni-heidelberg.de/tap'
    
    service = vo.dal.TAPService(link)
    
    result=service.run_async(adql_query.format(**locals()),
                                                maxrec=160000)
                                                
    catalog=result.to_table()


LIFE Star Catalog
-----------------

One of the main use cases of the LIFE Target Database is to enable the creation of the LIFE Star Catalog (LIFE-StarCat). It is the stellar sample of potential targets that is used be LIFEsim provide mission yield estimates.

The creation of the fourth version (LIFE-StarCat4) requires python 3 as well as the packages numpy, pyvo and astropy. Below the detailed explanation of its construction:

The stellar sample was compiled from querying the LIFE database (Menti et al. 2021) for stellar objects within 30 pc. Appart from the new source we tried to keep the criteria for inclusion into the catalog as similar to the ones from the previous LIFE catalog version as possible. This allows us to compare the two versions better before changing the criteria in the next catalog version (LIFE-StarCat5) again. (Among others to include higher order systems since they could have better chances for hosting life but in current version not included due to complexity e.g. centaury system). We started by removed objects with luminosity class I, II and III as we focused on main-sequence stars (or close to). In case no luminosity class was given we assumed the objects were main-sequence objects. In order to assess which objects are members of binary or multiple systems we used the LIFE database's hierarchical link. This feature connects an object with its parent and child objects. We therefore removed objects which were neither single stars nor trivial binaries. Concretely we removed objects where either the parent object of the star is istself a child object of another system, where multiple parent objects were given for the same star or which had more than one companion. As a note of caution this classification in the LIFE database is heavily based on Washington Visual Double Star Catalog (WDS; Mason et al. 2001) which in term only describes visual multiples and therefore contains some objects that are not physically bound. We also removed objects with incomplete information for the stellar parameters (this step also included binaries, where one component was not listed as a main-sequence star) or no separation between the binaries was given. In order to ensure that all stellar components of the remaining binary systems could harbor stable planetary systems between 0.5 and 10 AU we transformed the angular separation into physical and used it together with the stability criterion from Holman & Wiegert (1999) that takes into account the stellar masses, their separation and the eccentricity of the binary orbit (which we assumed to be zero). Systems that did not fulfill this stability criterion were removed. Finally we introduced a new parameter for future mission architecture analysis in which we flag objects that are contained within \pm 45 degrees of the ecliptic plane. This simulates the loss of targets by not full sky field of view. The catalog consists of 3909 stars in total (434 wide binary components and 3475 single stars).


