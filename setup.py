"""
Setup script to install life_td as a Python package.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from os.path import join, dirname
from setuptools import find_packages, setup


# -----------------------------------------------------------------------------
# RUN setup() FUNCTION
# -----------------------------------------------------------------------------

# Get version from VERSION file
with open(join(dirname(__file__), \
    "data_generation/life_td_data_generation/VERSION")) as version_file:
    version = version_file.read().strip()

# Run setup()
setup(
    name='life_td_data_generation',
    version=version,
    description='life_td_data_generation: Data generation for LIFE Target Database',
    url='https://github.com/fmenti/life_td',
    install_requires=[
        'astropy ~= 7.0.0',
        'numpy ~= 2.2.2',
        'pyvo ~= 1.6',
        'matplotlib ~= 3.10.0',
    ],
    extras_required={
        'docs': [
            'furo',
            'sphinx',
        ]
    },
    packages=find_packages(),
    zip_safe=False,
)
