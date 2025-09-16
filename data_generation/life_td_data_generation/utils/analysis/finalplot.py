import matplotlib.pyplot as plt
import numpy as np  # arrays
from typing import Iterable, List, Sequence, Tuple  # typing helpers
from numpy.typing import ArrayLike  # NumPy array protocol types

# self created modules
from utils.io import Path
from utils.analysis.analysis import (
    spectral_type_histogram,
    x_position,
    SpectralType,
)


def get_distance_cut(
    distance_samples: Sequence[ArrayLike],
) -> float:
    """
    Compute the maximum distance across all distance samples.

    Each element of ``distance_samples`` is expected to be an array-like
    (e.g., numpy array or column) of distances in parsec.

    :param distance_samples: Sequence of per-catalog distance arrays.
    :type distance_samples: Sequence[ArrayLike]  # list of 1D numeric arrays
    :returns: Maximum distance over all samples.
    :rtype: float
    """
    distance_cut = 0.0
    for sample in distance_samples:
        # Ensure sample is array-like; relies on numpy max for speed/robustness
        try:
            sample_max = float(np.max(sample))
        except ValueError:
            # Empty sequence: skip
            continue
        if sample_max > distance_cut:
            distance_cut = sample_max
    return distance_cut


def define_xticks(distance_cut: float) -> Tuple[List[str], float]:
    """
    Prepare x-ticks and bin step size for the distance-binned subplot.

    Bins are chosen as 5 pc for cuts <= 35 pc, and 10 pc for cuts
    <= 50 pc. Distances larger than 50 pc are not supported here.

    :param distance_cut: Maximum distance in parsec.
    :type distance_cut: float
    :returns: Tuple of x-tick labels and bin step size in parsec.
    :rtype: Tuple[List[str], float]
    """
    if distance_cut <= 35.0:
        xticks_total = ["0-5", "5-10", "10-15", "15-20", "20-25", "25-30"]
        stepsize = 5.0
    elif distance_cut <= 50.0:
        xticks_total = ["0-10", "10-20", "20-30", "30-40", "40-50"]
        stepsize = 10.0
    else:
        # Preserve original behavior (print) instead of raising
        print("Error: distance bigger than 50pc not programmed.")
        # Fallback to 10 pc bins up to 50 for graceful degradation
        xticks_total = ["0-10", "10-20", "20-30", "30-40", "40-50"]
        stepsize = 10.0

    # Limit number of bins according to cutoff
    xticks = xticks_total[: int(distance_cut / stepsize)]
    return xticks, stepsize


def myround(x: float, base: int = 5) -> float:
    """
    Round a value to the nearest multiple of a given base.

    :param x: Value to round.
    :type x: float
    :param base: Base to which ``x`` is rounded.
    :type base: int
    :returns: Rounded value.
    :rtype: float
    """
    return base * round(x / base)


def _extract_samples(
    stars: Sequence[object],
) -> Tuple[List[ArrayLike], List[ArrayLike]]:
    """
    Extract spectral type and distance columns from star tables.

    This function expects that each element in ``stars`` behaves like an
    astropy Table (or similar) with two columns: first the spectral type
    identifier and second the distance in parsec. We keep this flexible
    to avoid a hard dependency here.

    :param stars: Tables with at least two columns:
                  spectral type and distance [pc].
    :type stars: Sequence[object]  # astropy-like tables
    :returns: Two lists containing the spectral types and distances.
    :rtype: Tuple[List[ArrayLike], List[ArrayLike]]
    """
    spectral_type_samples = [
        stellar_cat[stellar_cat.colnames[0]] for stellar_cat in stars
    ]
    distance_samples = [
        stellar_cat[stellar_cat.colnames[1]] for stellar_cat in stars
    ]
    return spectral_type_samples, distance_samples


def _compute_distance_binned_distributions(
    spectral_type_samples: Sequence[ArrayLike],
    distance_samples: Sequence[ArrayLike],
    xticks: Sequence[str],
    stepsize: float,
    spec_labels: Sequence[str],
) -> np.ndarray:
    """
    Compute spectral-type histograms per distance bin for each sample.

    :param spectral_type_samples: Per-sample spectral type arrays.
    :type spectral_type_samples: Sequence[ArrayLike]
    :param distance_samples: Per-sample distance arrays [pc].
    :type distance_samples: Sequence[ArrayLike]
    :param xticks: Text labels for the distance bins (e.g., "0-5").
    :type xticks: Sequence[str]
    :param stepsize: Bin width in parsec corresponding to ``xticks``.
    :type stepsize: float
    :param spec_labels: Ordered list of spectral types (e.g., O..M).
    :type spec_labels: Sequence[str]
    :returns: 3D array with shape (nsamples, nbins, nspec).
    :rtype: np.ndarray
    """
    nsamples = len(spectral_type_samples)
    nbins = len(xticks)
    nspec = len(spec_labels)

    sub_specdist = np.zeros((nsamples, nbins, nspec))
    for i in range(nsamples):
        distances_i = np.asarray(distance_samples[i])
        sptypes_i = np.asarray(spectral_type_samples[i])
        sample_max = np.max(distances_i) if distances_i.size else 0.0

        for j in range(nbins):
            low = j * stepsize
            high = (j + 1) * stepsize
            if low > sample_max:
                # No stars in bins beyond the sample's max distance
                sub_specdist[i][j] = [0.0] * nspec
            else:
                mask = (distances_i > low) & (distances_i < high)
                sub_specdist[i][j] = spectral_type_histogram(sptypes_i[mask])
    return sub_specdist


def starcat_distribution_plot(
    stars: Sequence[object],
    labels: Sequence[str],
    path: str = f"{Path().plot}final_plot.png",
) -> None:
    """
    Plot spectral-type distributions (overall and distance-binned).

    This generates a figure with two subplots:
    #. Left: spectral type histogram for each sample.
    #. Right: spectral type histogram per distance bin for each sample.

    Note:
        This plot is mainly retained for backward compatibility. For
        future work consider using a more structured visualization
        pipeline.

    :param stars: Sequence of tables (e.g., astropy Tables). Each table
                  must have its first column as spectral type and its
                  second as distance [pc].
    :type stars: Sequence[object]  # astropy-like tables
    :param labels: Legend labels corresponding to ``stars``.
    :type labels: Sequence[str]
    :param path: File path to save the PNG figure.
    :type path: str
    :returns: None
    :rtype: None
    """
    # Extract arrays
    spectral_type_samples, distance_samples = _extract_samples(stars)

    # Determine distance range and bins
    distance_cut_in_pc = myround(get_distance_cut(distance_samples))
    xticks, stepsize = define_xticks(distance_cut_in_pc)

    # Shared configuration
    width = 0.15  # width of the bars
    colors = ["tab:blue", "tab:orange", "tab:green"]  # limited palette

    # Initialize spectral labels and containers
    spec = [spectraltype.name for spectraltype in SpectralType]
    nsamples = len(spectral_type_samples)
    specdist = np.zeros((nsamples, len(spec)))

    # Prepare figure
    fig, ((ax1, ax2)) = plt.subplots(
        1,
        2,
        sharey="row",
        figsize=[16, 5],
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    # First subplot: overall spectral type distribution
    x = np.arange(len(spec))

    for i in range(nsamples):
        specdist[i] = spectral_type_histogram(spectral_type_samples[i])
        # Remove O and B for tighter plot: use [2:]
        ax1.bar(
            x_position(x[2:], nsamples, i),
            specdist[i][2:],
            width,
            align="center",
            label=labels[i],
            color=colors[i % len(colors)],
        )

    ax1.set_yscale("log")
    plt.sca(ax1)
    plt.xticks(x[2:], spec[2:])
    ax1.set_title("Spectral type distribution")
    ax1.set_ylabel("Number of stars")
    ax1.set_xlabel("Spectral types")
    ax1.legend(loc="upper left")

    # Second subplot: spectral type distribution vs. distance
    index = np.arange(len(xticks))  # number of bars per x-tick label

    sub_specdist = _compute_distance_binned_distributions(
        spectral_type_samples=spectral_type_samples,
        distance_samples=distance_samples,
        xticks=xticks,
        stepsize=stepsize,
        spec_labels=spec,
    )

    # Plot per spectral type, skipping O and B as above
    for i in range(nsamples):
        for l in np.arange(len(spec))[2:]:
            ax2.bar(
                x_position(
                    (len(xticks) + 1) * index,
                    len(stars) + (len(stars) + 2) * (len(spec) - 2),
                    (len(stars) + 2) * l + i,
                ),
                tuple(sub_specdist[i][:, l]),
                width,
                color=colors[i % len(colors)],
            )

    ax2.set_xlabel("Distance [pc]")
    plt.sca(ax2)
    plt.xticks((len(xticks) + 1) * index, xticks)
    ax2.set_title("Spectral type and distance distribution")

    plt.savefig(path, dpi=300)
    return
