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
with open(join(dirname(__file__), "life_td/VERSION")) as version_file:
    version = version_file.read().strip()

# Run setup()
setup(
    name='life_td_data_generation',
    version=version,
    description='life_td_data_generation: Data generation for LIFE Target Database',
    url='https://github.com/fmenti/life_td',
    install_requires=[
        'astropy ~= 5.3.0',
        'numpy ~= 1.24.3',
        'pyvo ~= 1.4.1',
        'matplotlib ~= 3.6.3',
    ],
        'docs': [
            'furo',
            'sphinx',
        ]
    },
    packages=find_packages(),
    zip_safe=False,
)
