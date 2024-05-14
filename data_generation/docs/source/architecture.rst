.. _architecture:

Architecture
============

.. _architecture_intro:

Introduction
------------

Let's start with a bit of terminology used when dealing with databases. A `database` is a collection of data. They support electronic storage and maniputlation of data. `Data` itself are facts related to any object in consideration. For example name, age, temperature. A picture, image, file, etc. can also be considered data. The LIFE Target Database is a `relational database`. That means the data is stored in tables linked to each other with relations. Those tables and relations can be visualized using a `data model`. The one for `life_td` is explained in the section below. 

.. _architecture_data_model:

Data Model
----------

We used the Unified Modeling Language (UML) to create our data model. Each table in the database is represented as a box called a class. The columns of the table are called attributes. 

.. image:: classdiagramexplanation.png

.. image:: db_data_providers.png

.. image:: data_model_2024.png

TBD: insert example measurement tables

TBD insert image from input over package to database

