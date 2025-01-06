.. _modules:

Modules
========

This page contains a list of all available modules and a short description of what they are used for. More details on the design of the database can be found in the :ref:`architecture` section. 

.. note::
    
   The building and provider modules are currently in the status of "working but difficult to understand". Some new functions and classes are in the process of being implemented in an attempt to improve the readability.


life_td module
--------------

The data for the LIFE database is generated using the function :py:func:`life_td.create_life_td`:

It is also the one being called when the module is run as a script. 

During development it can happen that one runs the data generation multiple times e.g. when one is debugging a particular provider. As running the whole data generation takes time especially with larger distance cuts I implemented the function :py:func:`life_td.partial_create`. There one can specify which provider data is generated and which one loaded from the previously generated files and can potentially save a lot of time.

If one only wants to load the data e.g. for analysis reasons the function :py:func:`life_td.load_life_td` can be used.

utils module
------------

The utils module contains funtions of general utility used in multiple places throughout life_td_data_generation. Most often used are the functions for reading and writing of data: :py:func:`utils.save` and :py:func:`utils.load`. The later one uses the function :py:func:`utils.stringtoobject` where string type data is transformed into object type one. This has the advantage of allowing strings of varying length. Without it truncation of entries is a risk.

The next functions are for analyzing catalogs. With the function :py:func:`utils.object_contained` one can find out if a given list of objects is contained in a specified catalog. This function is of particular use if one wants to find out in which step of a script a specific object got lost.

The function :py:func:`utils.compare_catalogs` compares two given catalogs by providing statistics about the amount of common objects and parameter distributions. 

If one wants to look only at objects that are present in both of the two catalogs one is analysing the function :py:func:`utils.create_common` can be used.
   
sdata class
-----------

The sdata class is an attempt at separating the structuring of the data from the building and provider modules. It contains the class :py:class:`sdata.structured_data`
which contains functions specifying table manipulytion, loading and saving.

The second class in the sdata class is :py:class:`sdata.provider` which specifies all the tables, columns and parameter types for the final database tables.
     
provider module
--------------- 

The provider module generates the data for the database for each of the data providers separately. 

The provider module contains all the functions dealing with data acquisition and pre structuring from individual data providers. It contains specific functions for each provider: :py:func:`provider.provider_sim` for SIMBAD, :py:func:`provider.provider_gk` for the sdb from Grant Kennedy, :py:func:`provider.provider_exo` for Exo-Mercat, :py:func:`provider.provider_life` for data adapted by our team, :py:func:`provider.provider_wds` for Washington Double Star Catalog, :py:func:`provider.provider_gaia` for Gaia. In addition it also contains functions utilized by those data provider functions. 

The provider module is called when one of the life_td modules functions needs to generate data.

building module
---------------

The building module combines the data from the individual data providers. This is done by the function :py:func:`building.building`. The module also contains specific functions for merging of individual parameters: :py:func:`building.idsjoin` for identifiers, :py:func:`building.objectmerging` for objects and a general table merging function :py:func:`building.merge_table`. Additionally it also contains functions (:py:func:`building.best_para` :py:func:`building.best_parameters_ingestion`) to create tables only containing the best measurement for each object.

analysis module
---------------

The analysis module contains functions for graphical display of the database tables as well as general overview of the contained data.

The most important function :py:func:`analysis.final_plot` creates two plots for visualization of the spectral distribution. 

