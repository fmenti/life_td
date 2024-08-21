.. _introduction:

Introduction
============

Documentation Motivation
------------------------

The main goal of this documentation is to make both the code used to 
create the LIFE Target Database as well as the database itself follow 
the FAIR Principle.

FAIR
^^^^

TBD distinguish between FAIR data and FAIR software

**Findability** Contains richly annotated metadata to make clear what 
the data describes.

**Accessibility** Free and open source. The whole code can be accessed 
through `GitHub <https://github.com/fmenti/life_td>`_. Once run the 
code produces all the data that goes into the database. The database 
itself can be accessed through topcat or python (Check out the 
:ref:`LIFE Target Database Tutorial <tutorials>`).

**Interoperability** Independent of operating system. To access the 
database either python or Topcat is needed. Following the 
`Virtual Observatory <https://ivoa.net/>`_ standards.

**Reuse** Provenance in form of sources and code provided to enable 
reproducing everything.

Database Motivation
-------------------

Broader context - The LIFE initiative
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This project is part of the larger `LIFE <https://life-space-mission.com/>`_ initiative standing for Large Interferometer For Exoplanets. A brief overview is given here: 
LIFE was initiated in 2017 and officially kicked-off in 2018 with the goal to develop the science, the technology and a roadmap for an ambitious space mission that will allow humankind to detect and characterize the atmospheres of hundreds of nearby extrasolar planets including dozens that are similar to Earth. Thanks to NASA’s Kepler and TESS missions and dedicated, long-term exoplanet searches from the ground, we know that small, rocky exoplanets are ubiquitous in the Milky Way and also in the immediate Solar neighbourhood. Detecting the nearest exoplanets, scrutinizing their atmospheres and searching for habitable conditions and indications of biological activity is a cornerstone of 21st century astrophysics.

LIFE Target Database
^^^^^^^^^^^^^^^^^^^^

The Target Database for LIFE and ancillary observations is a
repository for relevant data about potential stellar target systems.
Our primary objectives include offering easy database access and
extracting input target catalogs to predict exoplanet detection
yields. Simultaneously, we aim to determine important missing
observations, offer data to derive mission design decisions,
and provide context for the analysis of data obtained by LIFE.

LIFE samples
^^^^^^^^^^^^

+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+--------------------------------+------+
| Sampel                | Goal                                                                                                                                                             | Assumption                                                                                          | Object Types                   | #    |
+=======================+==================================================================================================================================================================+=====================================================================================================+================================+======+
| Target Database       | As much data as possible on potential targets                                                                                                                    | 30pc, no single brown or white dwarfs                                                               | Stars, Systems, Planets, Disks | ~10⁴ |
+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+--------------------------------+------+
| Star-Cat              | As much data as needed for LIFEsim                                                                                                                               | mature single stars + wide binaries (stable planet orbit in hz)                                     | Stars                          | ~10³ |
+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+--------------------------------+------+
| LIFEsim output sample | Show for which mission parameters e.g. S/N what stellar objects (e.g. amount, to what distance) would be retrieved                                               | highest simulated observation yield                                                                 | Stars                          | ~10² |
+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+--------------------------------+------+
| Golden Targets        | targets that allow detailed characterization in relatively short observation times for oportunity of potential "time-resolved" experiments e.g. seasonal changes | best a priori known candidates e.g. proxima, Teegarden's star/planets, stars also observed with HWO | Stars                          | ~10¹ |
+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+--------------------------------+------+
| Final Targets         | best stars for LIFE mission                                                                                                                                      | golden targets + some of LIFEsim output sample                                                      | Stars                          | ~10² |
+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------+--------------------------------+------+




Complementary databases and catalogs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TBD (HOSTS SPORES, Starchive, HOSTS, HPIC)

