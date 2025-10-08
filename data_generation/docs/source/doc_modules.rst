.. _modules:

Modules
========

This page contains a list of all available modules and a short description of what they are used for. More details on the design of the database can be found in the :ref:`architecture` section.

.. note::
    This page is slightly outdated. If you are interested in an up to date version please contact me and I will increase the priority for updating this documentation section.


life_td module
--------------

The data for the LIFE database is generated using the function :py:func:`life_td.create_life_td`:

It is also the one being called when the module is run as a script.

During development it can happen that one runs the data generation multiple times e.g. when one is debugging a particular provider. As running the whole data generation takes time especially with larger distance cuts I implemented the function :py:func:`life_td.partial_create`. There one can specify which provider data is generated and which one loaded from the previously generated files and can potentially save a lot of time.

If one only wants to load the data e.g. for analysis reasons the function :py:func:`life_td.load_life_td` can be used.

utils subpackage
----------------

The utils subpackage contains funtions of general utility used in multiple places throughout life_td_data_generation. Most often used are the functions in the io module for reading and writing of data: :py:func:`utils.io.save` and :py:func:`utils.io.load`. The later one uses the function :py:func:`utils.io.stringtoobject` where string type data is transformed into object type one. This has the advantage of allowing strings of varying length. Without it truncation of entries is a risk.

The other module is for analyzing catalogs. With the function :py:func:`utils.catalog_comparison.object_contained` one can find out if a given list of objects is contained in a specified catalog. This function is of particular use if one wants to find out in which step of a script a specific object got lost.

The function :py:func:`utils.catalog_comparison.compare_catalogs` compares two given catalogs by providing statistics about the amount of common objects and parameter distributions. 

If one wants to look only at objects that are present in both of the two catalogs one is analyzing the function :py:func:`utils.catalog_comparison.create_common` can be used.
   
sdata module
------------

The sdata module is an attempt at separating the structuring of the data from the building and provider modules. It contains different dictionaries which can be used as template for data containers.This is useful for example by having the tables, columns and parameter types for the final database tables already specified.
     
building module
---------------

The building module combines the data from the individual data providers. This is done by the function :py:func:`building.building`. The module also contains specific functions for merging of individual parameters: :py:func:`building.idsjoin` for identifiers, :py:func:`building.objectmerging` for objects and a general table merging function :py:func:`building.merge_table`. Additionally it also contains functions (:py:func:`building.best_para` :py:func:`building.best_parameters_ingestion`) to create tables only containing the best measurement for each object.

The building module is called by the life_td module.

provider subpackage
-------------------

The provider subpackage generates the data for the database for each of the data providers separately. 

The provider subpackage contains all the modules dealing with data acquisition and pre structuring from individual data providers. Those provider modules contain specific functions for each provider: :py:func:`provider.provider_simbad` for SIMBAD, :py:func:`provider.provider_sdb` for the sdb from Grant Kennedy, :py:func:`provider.provider_exo` for Exo-Mercat, :py:func:`provider.provider_life` for data adapted by our team, :py:func:`provider.provider_wds` for Washington Double Star Catalog, :py:func:`provider.provider_gaia` for Gaia. In addition the provider subpackage also contains a utils module used by the provider modules. 

The provider modules are called when the building module functions need to generate data.

analysis module
---------------

The analysis module contains functions for graphical display of the database tables as well as general overview of the contained data. Currently some of those functions are being transformed into unit and integration tests. As a consequence the analysis module is currently not working but will be repaired soon.

The most important function :py:func:`analysis.final_plot` creates two plots for visualization of the spectral distribution.

