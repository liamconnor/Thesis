import numpy
import matplotlib
matplotlib.use("PDF")
import pylab

# ----------------------------------------------------------------------
# PARAMETERS
# ----------------------------------------------------------------------

# I/O
data_file = "data/results.dat"
plot_file = "../fig4b.pdf"

# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

def ReadResults(data_file):
    """ 
    Read the results from the bootstrap analysis.
    """

    mass, m, merr, b, berr = numpy.genfromtxt(data_file, dtype="float64", comments="#", usecols=(0,5,6,7,8), unpack=True)

    return mass, m, merr, b, berr

def PlotResults(plot_file, mass, m, merr, b, berr):
    """
    Plot the slope and intercept as a function of halo mass.
    """

    #
    # Get ranges
    #

    xmin  = 12.0 
    xmax  = 15.75
    ymin1 = -0.0025
    ymax1 =  0.008
    ymin2 = -0.0105
    ymax2 = -0.0045 

    ms = 4.0 

    xlabel  = r"${\rm log}M\ [M_\odot]$"
    ylabel1 = r"${\rm Slope}$"
    ylabel2 = r"${\rm Intercept}$"

    #
    # Make the plot
    #

    fig_width_pt  = 300.
    inches_per_pt = 1. / 72.27
    golden_mean   = (numpy.sqrt(5.) - 1.) / 2.
    fig_width     = fig_width_pt * inches_per_pt
    fig_height    = fig_width #* golden_mean
    fig_size      = [fig_width, fig_height]
    params        = {'backend': 'ps', 'axes.labelsize': 13, 'text.fontsize': 10, \
                    'legend.fontsize': 9, 'xtick.labelsize': 12, 'ytick.labelsize': 12, \
                    'text.usetex': True, 'figure.figsize': fig_size}

    pylab.rcParams.update(params)
    pylab.figure(1)
    pylab.clf()

    ax1 = pylab.subplot(211)

    pylab.hlines(0., xmin, xmax, linestyles="dashed", colors="k")
    pylab.errorbar(mass, m, fmt="o", yerr=merr, color="b", ms=ms)
    pylab.ylabel(ylabel1)

    ax2 = pylab.subplot(212)

    pylab.hlines(0., xmin, xmax, linestyles="dashed", colors="k")
    pylab.errorbar(mass, b, fmt="o", yerr=berr, color="b", ms=ms)
    pylab.ylabel(ylabel2)
    pylab.xlabel(xlabel)

    #
    # Adjust plot
    #

    ax1.set_xlim((xmin, xmax))
    ax2.set_xlim((xmin, xmax))
    ax1.set_ylim((ymin1, ymax1))
    ax2.set_ylim((ymin2, ymax2))

    pylab.setp(ax1.get_xticklabels(), visible=False)

    pylab.subplots_adjust(left=0.20, right=0.97, bottom=0.12, top=0.97, hspace=0.)

    #
    # Save the plot
    #

    pylab.savefig(plot_file)
    pylab.close()

# ----------------------------------------------------------------------
# MAIN 
# ----------------------------------------------------------------------

#
# Read the relevant data
#

mass, m, merr, b, berr = ReadResults(data_file)

#
# Plot the data
#

PlotResults(plot_file, mass, m, merr, b, berr)

