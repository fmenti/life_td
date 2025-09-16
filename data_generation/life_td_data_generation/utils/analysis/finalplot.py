import matplotlib.pyplot as plt
import numpy as np  # arrays


# self created modules
from utils.io import Path
from utils.analysis.analysis import spectral_type_histogram,x_position,SpectralType

def get_distance_cut(distance_samples):
    """
    Gets the maximum distance of the sample of catalogs.
    """
    distance_cut = 0
    for i in range(len(distance_samples)):
        if max(distance_samples[i]) > distance_cut:
            distance_cut = max(distance_samples[i])
    return distance_cut

def define_xticks(distance_cut: float):
    """
    Prepares the xticks of the second subfigure in the final_plot.

    :param float distance_cut: Max distance of sample in parsec.
    :returns:
    :rtype:
    """
    if distance_cut <= 35.0:
        xticks_total = ["0-5", "5-10", "10-15", "15-20", "20-25", "25-30"]
        stepsize = 5.0
    elif distance_cut <= 50.0:
        xticks_total = ["0-10", "10-20", "20-30", "30-40", "40-50"]
        stepsize = 10.0
    else:
        print("Error: distance bigger than 50pc not programmed.")

    xticks = xticks_total[: int(distance_cut / stepsize)]
    return xticks, stepsize


def myround(x, base=5):
    return base * round(x / base)

def starcat_distribution_plot(stars, labels, path=f"{Path().plot}final_plot.png"):
    """
    Plots spectral distribution in two subplots.

    This plot is a bit of a mess, in the future use the function
    spectral_type_histogram_catalog_comparison instead.
    Makes plot with two subfigures, first a histogram of spectral
    type and then one of spectral type for each distance sample of
    0-5, 5-10, 10-15 and 15-20pc.

    :param stars: list of astropy.table.table.Table of shape (,2)
        with columns first spectral type and then distance in pc.
    :param labels: list containing the labels for the plot
    :param path: location to save the plot
    """
    # tbd: replace stars with two variables to prevent having [0] hardcoded stuff
    # list of arrays of spectral type info: spectral_type_samples
    # list of arrays of spectral type info: distance_samples

    spectral_type_samples = [
        stellar_cat[stellar_cat.colnames[0]] for stellar_cat in stars
    ]
    distance_samples = [
        stellar_cat[stellar_cat.colnames[1]] for stellar_cat in stars
    ]

    distance_cut_in_pc = myround(get_distance_cut(distance_samples))
    xticks, stepsize = define_xticks(distance_cut_in_pc)

    # variables for both subfigures
    width = 0.15  # with of the bars
    color = ["tab:blue", "tab:orange", "tab:green"]

    # ititializes the data containers
    spec = [spectraltype.name for spectraltype in SpectralType]
    specdist = np.zeros((len(spectral_type_samples), len(spec)))

    # prepares subfigure display
    fig, ((ax1, ax2)) = plt.subplots(
        1,
        2,
        sharey="row",
        figsize=[16, 5],
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    # variables for first subfigure
    x = np.arange(len(spec))

    # first subfigure
    # remove spectraltypes that are not present (O and B) for tighter plot
    # => x -> x[2:], specdist[i] -> specdist[i][2:]

    # first_sub_figure(ax1)
    for i in range(len(spectral_type_samples)):
        specdist[i] = spectral_type_histogram(spectral_type_samples[i])
        # shift x so that not all samples plotted over each other
        ax1.bar(
            x_position(x[2:], len(spectral_type_samples), i),
            specdist[i][2:],
            width,
            align="center",
            label=labels[i],
            color=color[i],
        )

    ax1.set_yscale("log")
    # problem beim zweiten plot: min() arg is an empty sequence
    plt.sca(ax1)
    plt.xticks(x[2:], spec[2:])
    ax1.set_title("Spectral type distribution")
    ax1.set_ylabel("Number of stars")
    ax1.set_xlabel("Spectral types")
    ax1.legend(loc="upper left")

    # second subfigure
    index = np.arange(len(xticks))  # wie viele bars pro label es haben wird
    # n = np.arange(7)#wieviele labels
    sub_specdist = np.zeros(
        (len(spectral_type_samples), len(xticks), len(spec))
    )
    # def distance_binned_spectral_distributions(spectral_type_samples,distance_samples):

    # return sub_specdist

    for i, spectral_type_sample in enumerate(spectral_type_samples):
        # what am I doing here? I create the distance binned spectral distributions
        for j in range(len(xticks)):
            # here have an error because lifestarcat only goes up to 20pc
            if j * stepsize > max(distance_samples[i]):
                sub_specdist[i][j] = [0.0] * 7
            else:
                sub_specdist[i][j] = spectral_type_histogram(
                    spectral_type_samples[i][
                        np.where(
                            (distance_samples[i] > j * stepsize)
                            * (distance_samples[i] < (j + 1) * stepsize)
                        )
                    ]
                )
        # here I put in the axis bars
        for l in np.arange(len(spec))[2:]:
            # print(sub_specdist[i][:,l],
            # 5*index-((n+(n+2)*(s-2))*width)+((n+2)*l+i)*width)
            # problem if in one bin all 0
            # make if clause for it

            ax2.bar(
                x_position(
                    (len(xticks) + 1) * index,
                    len(stars) + (len(stars) + 2) * (len(spec) - 2),
                    (len(stars) + 2) * l + i,
                ),
                tuple(sub_specdist[i][:, l]),
                width,
                color=color[i],
            )

    ax2.set_xlabel("Distance [pc]")
    plt.sca(ax2)
    plt.xticks((len(xticks) + 1) * index, (xticks))  # (xticks_name))
    ax2.set_title("Spectral type and distance distribution")
    plt.savefig(path, dpi=300)
    return
