This project creates the LIFE target database (life_td).

This project is structured as follows:
q.rd: Resource descriptor used by DaCHS to create the life_td from the data.
data: Folder containing the input data (as votable xml files) for the life_td.
data_generation: Folder containing the code (as python files) used to create the data.

How and when to run the data_generation files:
To update the life_td (for bug fixes, new data, new features,...) the developer runs 

import life_td as ltd
sim, gk, wds, exo, life, gaia, database_tables= ltd.create_life_td(30)

in an environment with preinstalled numpy, astropy and pyvo.
