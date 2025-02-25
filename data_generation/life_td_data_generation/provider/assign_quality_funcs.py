""" 
Generates the data for the database for each of the data providers separately. 
"""

import numpy as np  #arrays
from astropy.table import  MaskedColumn
from provider.utils import lower_quality


def teff_st_spec_assign_quality(gaia_mes_teff_st_spec):
    interval = 41 * 9 / 5.
    gaia_mes_teff_st_spec['teff_st_qual'] = ['?' for j in range(len(gaia_mes_teff_st_spec))]
    for i, flag in enumerate(gaia_mes_teff_st_spec['flags_gspspec']):
        summed = 0
        for j in flag:
            summed += int(j)
        if summed in range(0, int(interval) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'A'
        elif summed in range(int(interval) + 1, int(interval * 2) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'B'
        elif summed in range(int(interval * 2) + 1, int(interval * 3) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'C'
        elif summed in range(int(interval * 3) + 1, int(interval * 4) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'D'
        elif summed in range(int(interval * 4) + 1, int(interval * 5) + 1):
            gaia_mes_teff_st_spec['teff_st_qual'][i] = 'E'
    return gaia_mes_teff_st_spec

def assign_quality_elementwise(exo_helptab, para, i):
    qual = 'B'
    if exo_helptab[para + '_max'][i] == 1e+20:
        qual = lower_quality(qual)
    if exo_helptab[para + '_min'][i] == 1e+20:
        qual = lower_quality(qual)
    return qual

def exo_assign_quality(exo_helptab):
    for para in ['mass', 'msini']:
        exo_helptab[para + '_pl_qual'] = MaskedColumn(dtype=object, length=len(exo_helptab))
        exo_helptab[para + '_pl_qual'] = ['?' for j in range(len(exo_helptab))]
        for i in range(len(exo_helptab)):
            exo_helptab[para + '_pl_qual'][i] = assign_quality_elementwise(exo_helptab, para, i)
    return exo_helptab

def assign_quality(table, column='', special_mode=''):
    # Define a mapping of special modes to their respective handler functions
    mode_handlers = {
        'teff_st_spec': lambda t, _: teff_st_spec_assign_quality(t),
        'exo': lambda t, _: exo_assign_quality(t),
        'gaia_binary': lambda t, col: assign_gaia_binary_quality(t, col),
        'wds_sep1': lambda t, col: assign_wds_sep_quality(t, col, mode='wds_sep1'),
        'wds_sep2': lambda t, col: assign_wds_sep_quality(t, col, mode='wds_sep2'),
    }

    # Check if special_mode exists in handlers
    if special_mode in mode_handlers:
        return mode_handlers[special_mode](table, column)

    # Handle default or special edge cases
    quality = _default_quality(column, special_mode)
    table[column] = [quality] * len(table)
    return table


# Helper function for Gaia binary case
def assign_gaia_binary_quality(table, column):
    table[column] = [
        'B' if table['binary_flag'][j] == 'True' else 'E'
        for j in range(len(table))
    ]
    return table


# Helper function for WDS separation cases
def assign_wds_sep_quality(table, column, mode):
    if column not in table.colnames:
        table[column] = [''] * len(table)  # Initialize column if it doesn't exist
    if mode == 'wds_sep1':
        table[column][:] = [
            'C' if not isinstance(j, np.ma.core.MaskedConstant) else 'E'
            for j in table['sep_ang_obs_date']
        ]
    elif mode == 'wds_sep2':
        table[column][:] = [
            'B' if not isinstance(j, np.ma.core.MaskedConstant) else 'E'
            for j in table['sep_ang_obs_date']
        ]
    return table



# Encapsulate default/assumption logic for fallback cases
def _default_quality(column, special_mode):
    if column == 'coo_gal_qual':
        # Quality transformation logic for ra/dec to gal
        return '?'
    elif special_mode in ['teff_st_phot', 'radius_st_flame', 'mass_st_flame']:
        return 'B'
    elif special_mode in ['model', 'wds_binary']:
        return 'C'
    elif special_mode == 'sim_binary':
        return 'D'
    return '?'  # Default to unknown quality flag


