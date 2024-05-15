.. _installation:

Installation
============

This short guide will walk you through the required steps to set up and install `life_td_data_generation`.

The code in this repository is organized as a Python package named `life_td_data_generation` together with a set of scripts that use the functions and classes of the package. The code was written for `Python <https://www.python.org>`_ 3.10 and above; earlier versions will likely require some small modifications.
To get started, you will need to install the package.
For this, we *strongly* recommend you to use a `virtual environment <https://virtualenv.pypa.io/en/latest/>`_. 

Virtual Environment
-------------------
First install ``virtualenv``, for example with the apt (linux) package manager:

.. code-block:: console

    $ apt install virtualenv
    
Then create a virtual environment:

.. code-block:: console

    $ virtualenv -p python3 folder_name

And activate the environment with:

.. code-block:: console

    $ source folder_name/bin/activate

A virtual environment can be deactivated with:

.. code-block:: console

    $ deactivate

.. important::
   Make sure to adjust the path where the virtual environment is installed and activated.

Installing the life_td package
------------------------------

Once you have set up a suitable virtualenv, clone this repository and install `life_td_data_generation` as a Python package:

.. code-block:: console

    git clone git@github.com:fmenti/life_td.git ;
    cd life_td/data_generation ;
    pip install .
    
Once a local copy of the repository exists, new commits can be pulled from Github with:

.. code-block:: console

    $ git pull origin main
    
Testing life_td_data_generation
-------------------------------

The installation can be tested by running the life_td module as a script:

.. code-block:: console

    $ python life_td.py
    
