.. _installation:

Installation
============

This short guide will walk you through the required steps to set up and install `life_td_data_generation`.


Installing the life_td package
------------------------------

The code in this repository is organized as a Python package named `life_td_data_generation` together with a set of scripts that use the functions and classes of the package.
To get started, you will need to install the package.
For this, we *strongly* recommend you to use a `virtual environment <https://virtualenv.pypa.io/en/latest/>`_. 

.. note:: 

    Attention: The code was written for **Python 3.10 and above**; earlier versions will likely require some small modifications.

Once you have set up a suitable virtualenv, clone this repository and install `life_td_data_generation` as a Python package:

.. code-block:: bash

    git clone git@github.com:fmenti/life_td.git ;
    cd life_td/data_generation ;
    pip install .
