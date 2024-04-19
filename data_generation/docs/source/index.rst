.. life_td documentation master file, created by
   sphinx-quickstart on Wed Mar 20 14:47:58 2024.
   You can adapt this file completely to your liking, but it should at 
   least contain the root `toctree` directive.

.. _index:

Welcome!
========


This is the documentation of **life_td**, a Python package that generates the data for the `LIFE Target Database <https://
dc.zah.uni-heidelberg.de/voidoi/q/lp/custom/10.21938/
ke_e6lzO_jjX_vzvVIcwZA>`_. The database contains information useful for the 
planned `LIFE mission <https://life-space-mission.com/>`_ (mid-ir, 
nulling interferometer in space). It characterizes possible target 
systems including information about stellar, planetary and disk 
properties. The data itself is mainly a collection from different other 
databases and catalogs. 

In this documentation, we explain the structure of the code base and provide some guides to help you get started if you want to reproduce our results or use our code for your own experiments.

We have tried our best to follow best practices for software development and document all important steps. However, at the end of the day, it is still research code. Therefore, if you do get stuck with anything, please feel free to reach out to us!

This applies, of course, also if you find any bugs or problems (either in the code or the docs)!

Check out the :doc:`usage` section for further information on how to use the database.


.. note::
    
    This project is under active development.

.. toctree::
   :maxdepth: 2
   :caption: Getting started
   :hidden:

   installation
..   tutorials/first_example.ipynb

.. toctree::
   :maxdepth: 2
   :caption: Documentation
   :hidden:
   
   introduction
   architecture
   creation
   usage
   modules

.. toctree::
   :maxdepth: 2
   :caption: About
   :hidden:
   
   about

