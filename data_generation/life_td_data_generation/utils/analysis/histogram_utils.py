import numpy as np
from enum import Enum
from dataclasses import dataclass

@dataclass
class Plotparas:
    width = 0.15
    color = ["tab:blue", "tab:cyan", "tab:green", "tab:olive"]


class SpectralType(Enum):
    O = "O"
    B = "B"
    A = "A"
    F = "F"
    G = "G"
    K = "K"
    M = "M"

def spectral_type_histogram(spectypes):
    """
    Makes a histogram of the spectral type distribution.

    :param spectypes: array containing a spectral type for each star
    :type spectypes:
    :returns: array containig the number of stars for
        each spectral type in spec.
    :rtype: np.array
    """
    specdist = np.zeros(len(SpectralType))
    for spectype in spectypes:
        if spectype not in ["", "nan"]:
            for j, spectraltype in enumerate(SpectralType):
                if spectype[0] == spectraltype.name:
                    specdist[j] += 1
    specdist = specdist.astype(int)
    return specdist

def x_position(x, n_samples, sample_index):
    """
    Calculates x position for plots so that not all samples plotted over each other.
    """
    width_of_samples = n_samples * Plotparas.width
    sample_location = sample_index * Plotparas.width
    # if I change sign in front of sample_location result is that catalogs get shown in order backwards
    return x - width_of_samples / 2 + sample_location


