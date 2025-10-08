import matplotlib.pyplot as plt
import numpy as np  # arrays
from utils.analysis.analysis import (
    distancecut_subplot,
    full_subplot,
    tight_plot,
)
from utils.analysis.histogram_utils import (
    SpectralType,
    spectral_type_histogram,
    x_position,
)
from utils.io import Path


def spectral_type_histogram_catalog_comparison(
    stellar_catalogs,
    labels,
    path=Path().plot + "sthcc.png",
    distance_colname="dist",
    spectral_type_colname="spec",
):
    """
    Creates the figure for the RNAAS article.

    Spectral type distribution of catalogs with shading for amount below 20pc.
    """
    spectral_type_samples = [
        stellar_cat[spectral_type_colname] for stellar_cat in stellar_catalogs
    ]

    # ititializes the data containers
    spec = [spectraltype.name for spectraltype in SpectralType]
    # initializes array to contain spectral distribution of total samples
    specdist_total = np.zeros((len(spectral_type_samples), len(spec)))
    # initializes array to contain spectral distribution of sub samples
    specdist_sub = np.empty_like(specdist_total)

    plt.figure(figsize=(8, 4))
    x = np.arange(len(spec))

    for i in range(len(spectral_type_samples)):
        sample_total = stellar_catalogs[i][spectral_type_colname]
        specdist_total[i] = spectral_type_histogram(sample_total)

        x_pos = x_position(x, len(spectral_type_samples), i)

        full_subplot(x_pos, specdist_total[i], i, labels)

        # now make same for only within 20pc sample
        sample_sub = sample_total[
            np.where(stellar_catalogs[i][distance_colname] < 20.0)
        ]
        specdist_sub[i] = spectral_type_histogram(sample_sub)

        distancecut_subplot(x_pos, specdist_sub[i], i, "//")

        sample_sub = sample_total[
            np.where(stellar_catalogs[i][distance_colname] < 10.0)
        ]
        specdist_sub[i] = spectral_type_histogram(sample_sub)

        distancecut_subplot(x_pos, specdist_sub[i], i, "////")

    # creating a single label for the tree hatch barplots
    plt.bar(
        [1],
        [0],
        log=True,
        edgecolor="black",
        hatch="//",
        label=["d < 20pc"],
        facecolor="none",
    )
    plt.bar(
        [1],
        [0],
        log=True,
        edgecolor="black",
        hatch="////",
        label=["d < 10pc"],
        facecolor="none",
    )

    xt, yt = tight_plot(x, spec)
    plt.xticks(xt, yt)
    plt.title("Spectral type distribution")
    plt.ylabel("Number of stars")
    plt.xlabel("Spectral types")
    plt.legend(loc="upper left")
    plt.savefig(path, dpi=300)
    return
