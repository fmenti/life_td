import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.stats import norm
from utils.analysis.analysis import get_data, database_para_temp_class_plot_prep
from utils.io import Path, load

# not showing plots
# show_kw=True #<--- added
show_kw = False

if show_kw:  # <--- added
    plt.ion()  # <--- added

with open(Path().data + "distance_cut.txt") as h:
    distance_cut = float(h.readlines()[0])


def model_exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c


def plot_data_and_fit(title, data, p0,fit=True, bins = 10):
    """
    Plot histogram of data with an exponential decay fit.

    :param str title: Plot title
    :param array-like data: Data to plot and fit
    :param list p0: Initial parameters [a, b, c] for exponential decay fit: a * exp(-b * x) + c
    """
    # Setup plot
    plt.figure()
    plt.title(title)

    # Create histogram
    bin_heights, bin_borders, _ = plt.hist(
        data, bins=bins, density=True, label="pdf of data"
    )
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

    if fit:
        # Fit exponential decay model
        popt, _ = curve_fit(model_exp_decay, bin_centers, bin_heights, p0=p0)
        a_opt, b_opt, c_opt = popt

        # Create fitted curve
        x_model = np.linspace(bin_centers[0], max(bin_borders), 100)
        y_model = model_exp_decay(x_model, a_opt, b_opt, c_opt)


        # Plot fitted curve
        plt.plot(x_model, y_model, color="r", label="fit")

    # Add labels and legend
    plt.ylabel("amount of objects")
    plt.legend(loc="upper right")
    plt.show()


def ravsdec(x_label, y_label, x, y):
    plt.figure()
    ra = np.linspace(0, 360)
    plt.scatter(x, y, s=2)
    # plt.scatter(within45deg['coo_ra'],within45deg['coo_dec'],s=2)
    # ecliptic plane in equatorial coordinates
    if x_label == "coo_ra":
        ecliptic = (23.4) * np.sin(2 * np.pi * ra / 360)
        plt.plot(ra, ecliptic, color="k",label='ecliptic')
        plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title("coordinate distribution")
    # plt.savefig('plots/quasi_aitoff', dpi=300)
    plt.show()
    return


def norm_fit(data, title):
    bins = 10
    plt.hist(
        data, bins, density=True, label="pdf of data"
    )  # , alpha=0.6, color='g')
    # do I need to keep density=True here or can I use non normalized display?

    mu, std = norm.fit(data)
    y_fit = norm.pdf(np.sort(data), mu, std)

    plt.title(title)
    plt.plot(np.sort(data), y_fit, color="r", label="fit")
    plt.legend()
    plt.show()
    return mu


def test_data_makes_sense_main_id():
    # experimental numbers
    dist_dist = Table(
        data=[
            [5, 10, 30],
            [420, 1728, 22107],# does not scale beyond
            [205, 931, 15869],
            [101, 413, 4625],
            [113, 366, 1464],
            [6, 18, 149], # does not scale beyond
        ],
        names=["dist", "total", "st", "sy", "pl", "di"],
        dtype=[float, float, float, float, float, float],
    )

    [table] = load(["objects"], location=Path().data)
    data = table["main_id", "type"]

    spec = np.array(["System", "Star", "Exoplanet", "Disk"])

    plt.figure()
    plt.title(f"Object type distribution up to {distance_cut} pc")
    plt.xlabel("Number of objects")
    plt.hist(data["type"], histtype="bar", log=True, orientation="horizontal")
    plt.yticks(np.arange(4), spec)
    # plt.savefig(Path().plot+path, dpi=300)
    plt.show()

    number_total = len(data)
    number_st = len(data[np.where(data["type"] == "st")])
    number_sy = len(data[np.where(data["type"] == "sy")])
    number_pl = len(data[np.where(data["type"] == "pl")])
    number_di = len(data[np.where(data["type"] == "di")])

    print('distance: ', distance_cut)
    print('total: ', number_total)
    print('stars: ', number_st)
    print('systems: ', number_sy)
    print('planets: ', number_pl)
    print('disks: ', number_di)

    total = number_total ** (1 / 3) / distance_cut
    st = number_st ** (1 / 3) / distance_cut
    sy = number_sy ** (1 / 3) / distance_cut
    pl = number_pl ** (1 / 3) / distance_cut
    di = number_di ** (1 / 3) / distance_cut

    assert total < 2 and total > 0.5
    assert st < 1.5 and st > 0.5
    assert sy < 1.5 and sy > 0.4
    if distance_cut <= 30.:
        assert pl < 1.5 and pl > 0.3
        assert di < 0.4 and di > 0.1
    else:
        print('behavior <30pc not scalable for disks and planets')
    # tbd make this via analytical and not just experimental numbers


def test_data_makes_sense_mass_st():
    # data
    data = get_data("mes_mass_st", "mass_st_value")

    plot_data_and_fit("Stellar Mass", data, [1, 8, 0],
                      bins = 100
)

    #add one function that plots the masses binned by temperature class
    #ms_tempclass = np.array(["O", "B", "A", "F", "G", "K", "M"])
    #should have a lot already programmed in the analysis file

    # assert
    assert max(data) < 60  # O3V
    assert min(data) > 0.074  # brown dwarf

def test_data_makes_sense_mass_st_class():
    # data
    [table] = load(["star_basic"], location=Path().data)

    data = database_para_temp_class_plot_prep(table, "mass_st_value",
                                   "Stellar Mass [Msun]")

    # would also be nice to do it for all masses not just bestmass. need a join for that

    temp_class_list = np.array(["O", "B", "A", "F", "G", "K", "M"])
    maxlist=[60.,18.,2.7,1.75,1.10,0.92,0.6]
    minlist=[17.,2.,1.5,1.,0.85,0.55,0.074]
    # wait, those are ms numbers.
    # -> masses shouldn't change much with age just radii

    # issues, gaia mass vs simbad sptype?
    # mass A star too low -> find object in db and have closer look
    #bad qual sptype and gaia flame mass -> does gaia offer sptype?

    for tempclass,maxvalue,minvalue in zip(temp_class_list,maxlist,minlist):
        testcase = data["mass_st_value"][np.where(data["class_temp"] == tempclass)]
        if len(testcase)>0:
            assert max(testcase) < maxvalue
            assert min(testcase) > minvalue

def test_data_makes_sense_teff_st_class():
    # data
    [table] = load(["star_basic"], location=Path().data)

    data = database_para_temp_class_plot_prep(table, "teff_st_value",
                                   "Stellar Temperature [K]")

    # would also be nice to do it for all masses not just bestmass. need a join for that

    temp_class_list = np.array(["O", "B", "A", "F", "G", "K", "M"])
    maxlist=[45000.,31500.,10000.,7300,5950.,5300.,3900.]
    minlist=[31500.,10000.,7300.,5950.,5300.,3900.,2300.]
    # wait, those are ms numbers.
    # even around 10 MS M stars currently have too high temps from gaia-gspec

    for tempclass,maxvalue,minvalue in zip(temp_class_list,maxlist,minlist):
        testcase = data["teff_st_value"][np.where(data["class_temp"] == tempclass)]
        if len(testcase)>0:
            assert max(testcase) < maxvalue
            assert min(testcase) > minvalue


def test_data_makes_sense_temp_st():
    # want a scatter plot x axis distance and y axis temperature
    # not sure what to assert
    # not sure how to get data -> what to do about fill values?
    # data
    # loading the correct table
    [table] = load(["star_basic"], location=Path().data)
    # exctracting the correct columns
    arr = table["dist_st_value", "teff_st_value"]
    arr2 = arr[np.where(arr["dist_st_value"] != 1e20)]
    data = arr2[np.where(arr2["teff_st_value"] != 1e20)]

    # plt.figure()
    fig, ax = plt.subplots(
        figsize=(9, 6)
    )  # subplots so that I can overplot old version?

    ax.scatter(data["dist_st_value"], data["teff_st_value"], s=2)
    ax.set_yscale("log")

    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel("Temperature [K]")

    database_para_temp_class_plot_prep(table, "teff_st_value",
                                   "Stellar Temperature [K]")

    # assert
    assert max(data["teff_st_value"]) < 45000  # O3V
    assert min(data["teff_st_value"]) > 2300  # brown dwarf

def test_data_makes_sense_radius_st():
    # want a scatter plot x axis distance and y axis temperature
    # not sure what to assert
    # not sure how to get data -> what to do about fill values?
    # data
    # loading the correct table
    [table] = load(["star_basic"], location=Path().data)
    # exctracting the correct columns
    arr = table["dist_st_value", "radius_st_value"]
    arr2 = arr[np.where(arr["dist_st_value"] != 1e20)]
    data = arr2[np.where(arr2["radius_st_value"] != 1e20)]

    # plt.figure()
    fig, ax = plt.subplots(
        figsize=(9, 6)
    )  # subplots so that I can overplot old version?

    ax.scatter(data["dist_st_value"], data["radius_st_value"], s=2)
    ax.set_yscale("log")

    ax.set_xlabel("Distance [pc]")
    ax.set_ylabel("Radius [Rsun]")

    database_para_temp_class_plot_prep(table, "radius_st_value",
                                       "Stellar Radius [Rsun]")

    # assert
    assert min(data["radius_st_value"]) > 0.095  # brown dwarf


def test_data_makes_sense_mass_pl():
    # data
    data = get_data("mes_mass_pl", "mass_pl_value")

    plot_data_and_fit("Planetary Mass", data, [1, 0.1, 0.1])

    # assert
    assert max(data) < 75  # m star (in jupiter masses)
    assert min(data) > 0
    # not working for distance cut 20
    # one object with 515 Mjup which would actually be a m star
    # looks like they have 0 values in there


def test_data_makes_sense_coo():
    # data
    data_ra = get_data("star_basic", "coo_ra")
    data_dec = get_data("star_basic", "coo_dec")

    # make sky plot
    ravsdec("coo_ra", "coo_dec", data_ra, data_dec)

    assert min(data_ra) > 0
    assert max(data_ra) < 360
    assert min(data_dec) > -90
    assert max(data_dec) < 90


def test_data_makes_sense_coo_gal():
    # data
    data_l = get_data("star_basic", "coo_gal_l")
    data_b = get_data("star_basic", "coo_gal_b")

    # make sky plot
    ravsdec("coo_gal_l", "coo_gal_b", data_l, data_b)

    assert min(data_l) > 0
    assert max(data_l) < 360
    assert min(data_b) > -90
    assert max(data_b) < 90


def test_data_makes_sense_mag_i():
    # data
    data = get_data("star_basic", "mag_i_value")

    mu = norm_fit(data, "Stellar i-band magnitude")

    assert mu < 11 and mu > 7


def test_data_makes_sense_mag_j():
    # data
    data = get_data("star_basic", "mag_j_value")

    mu = norm_fit(data, "Stellar J-band magnitude")

    assert mu < 11 and mu > 7


def test_data_makes_sense_mag_k():
    # data
    data = get_data("star_basic", "mag_k_value")

    mu = norm_fit(data, "Stellar K-band magnitude")

    assert mu < 11 and mu > 6  # failed at 5pc


def test_data_makes_sense_plx():
    # data
    data = get_data("star_basic", "plx_value")

    plot_data_and_fit("Stellar Parallax", data, [1, 0.1, 0.1],fit = False)

    # assert
    assert min(data) > 1/(distance_cut+distance_cut/10.)*1000
    assert min(data) < 1000  # 1 pc equivalent marcsec
# not working, doin it differently in simbad:     plx_cut = plx_in_mas_cut - plx_in_mas_cut / 10.0
# maybe instead do a sharp cut and deal with boundary objects later?

def test_data_makes_sense_dist_st():
    # data
    data = get_data("star_basic", "dist_st_value")

    plot_data_and_fit("Stellar Distance", data, [1, -0.5, 0])

    # assert
    assert min(data) < distance_cut+distance_cut/10.
    assert min(data) > 1

def test_data_makes_sense_dist_rad():
    # data
    data = get_data("disk_basic", "rad_value")
    data=data.astype(float)

    plot_data_and_fit("Disk Radius", data, [1, 0.1, 0.1])

    #stuff not working yet because vertical lines and no other parameter in disk
    # assert
    #assert min(data) < distance_cut+distance_cut/10.
    #assert min(data) > 1
