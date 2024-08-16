.. _architecture:

Architecture
============

.. _architecture_intro:

Introduction
------------

Let's start with a bit of terminology used when dealing with databases. A `database` is a collection of data. They support electronic storage and maniputlation of data. `Data` itself are facts related to any object in consideration. For example name, age, temperature. A picture, image, file, etc. can also be considered data. The LIFE Target Database is a `relational database`. That means the data is stored in tables linked to each other with relations. Those tables and relations can be visualized using a `data model`. The one for `life_td` is explained in the section below. 


Use Cases
---------

(TBD add tree for use cases and implementation consequences)

We considered the following five use cases:

A: Catalog extraction for main science objective of LIFE
B: Easy access for fellow scientists
C: Low maintenace for database administrator to save costs.
D: Identify important missing data for future observation proposals.
E: Flexible database design enabling easy expansion in form of adding parameters (e.g. metalicity), providers (e.g. TESS) or adapt to other use cases that might come up in the future (e.g. planned non detection information).
F: Provide context for analysis of data obtained by LIFE.

Which have the following consequences on the database features:

A1: 30pc cut -> max distance at which we can in reasonable time observe desired planet
A2: Main stellar parameters (name, position, distance, spectral type, effective temperature, mass, radius)
A3: Multiplicity information to predict stable planetary orbits as well as take into account lower planet occurence rates in multiples.
A4: Disk information to predict observability. Exozodi as noise source.
A5: Planet information to reduce detection time where planets already known or ruled out from nondetection or Hot Jupiter relation to habitable zone planets.
A6: Best parameters of data collected from different providers.

B1: VO compatibility
B2: Example queries to mitigate ADQL knowledge deficites.
B3: tutorials to help create catalog for own science project

C1: Data ingestion as automated as possible. This is achieved by prefering VO compatible data providers over other databases. Those are in turn prefered over modeled data and as a last resort literature data is collected.
C2: GAVO published, they take case that server runs smoothly.
C3: Everything public and documented to have good knowledge transfer in case of administrator change.

E1: Careful analysis of not only short term use cases but also mid and long term ones. 
E2: Designing and reiterating of the data model. 




.. _architecture_data_model:

Data Model
----------

We used the Unified Modeling Language (UML) to create our data model. Each table in the database is represented as a box called a class. The columns of the table are called attributes. 

.. image:: classdiagramexplanation.png

.. image:: db_data_providers.png

.. image:: data_model_2024.png

TBD: insert example measurement tables

TBD insert image from input over package to database

