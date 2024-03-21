.. life_td documentation master file, created by
   sphinx-quickstart on Wed Mar 20 14:47:58 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to life_td's documentation!
===================================

**life_td** stands for the LIFE Target Database 
<https://dc.zah.uni-heidelberg.de/voidoi/q/lp/custom/10.21938/ke_e6lzO_jjX_vzvVIcwZA> and contains information useful
for the planned LIFE mission <https://life-space-mission.com/> (mid-ir, nulling
interferometer in space). It characterizes possible
target systems including information about stellar,
planetary and disk properties. The data itself is mainly
a collection from different other catalogs. This is the documentation about 
the python data generation code for life_td.

Check out the :doc:`usage` section for further information.

The main goal of this documentation is to make both the code used to create the
LIFE Target Database as well as the database itself follow the FAIR Principle.

FAIR
----

**Findability** Metadata

**Accessibility** Free and open source. The whole code can be accessed through 
GitHub <https://github.com/fmenti/life_td>. Once run the code produces all the 
data that goes into the database. The database itself can be accessed through 
topcat or python (Check out the :ref:`LIFE Target Database Tutorial <tutorial>`).

**Interoperability** Python, VO standards

**Reuse** provenance

.. note::
    
    This project is under active development.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   usage



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
