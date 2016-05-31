import numpy
import matplotlib
matplotlib.use("PDF")
import pylab

# ----------------------------------------------------------------------
# PARAMETERS
# ----------------------------------------------------------------------

# I/O
data_file = "data/tdmdnu06.dat"
plot_file = "../fig4a.pdf"

# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

def ReadData(data_file):
    """
    Returns the scatter points and best-fit line.
    """

    # Get scatter points
    x, y = numpy.genfromtxt(data_file, dtype="float32", comments="#", unpack=True)

    # Get line
    fr = open(data_file, "r")
    lines = fr.readlines()
    fr.close()

    for i in xrange(len(lines)):

        if "# slope  =" in lines[i]:

            cline = str.split(lines[i])
            m     = float(cline[-1])

        if "# inter  =" in lines[i]:

            cline = str.split(lines[i])
            b     = float(cline[-1])

    return x, y, m, b

def PlotDMvsDNU(plot_file, x, y, m, b):
    """
    Plots the relevant data.
    """

    #
    # Plot parameters
    #

    xmin = -0.15 
    xmax =  0.15
    ymin = -0.02
    ymax =  0.01

    xlabel = r"$\delta_\nu^{\rm rec} - \langle \delta_\nu^{\rm rec}(\delta_c^{\rm rec}) \rangle$"
    ylabel = r"$(M^{\rm TN}-M^{\rm TC})/M^{\rm TC}$"

    #
    # Line to plot
    #

    xfit = numpy.linspace(xmin, xmax, 100)
    yfit = m*xfit + b

    #
    # Initialize plot
    #

    fig_width_pt  = 300.
    inches_per_pt = 1. / 72.27
    golden_mean   = (numpy.sqrt(5.) - 1.) / 2.
    fig_width     = fig_width_pt * inches_per_pt
    fig_height    = fig_width
    fig_size      = [fig_width, fig_height]
    params        = {'backend': 'ps', 'axes.labelsize': 13, 'text.fontsize': 11, \
                    'legend.fontsize': 10, 'xtick.labelsize': 11, 'ytick.labelsize': 11, \
                    'text.usetex': True, 'figure.figsize': fig_size}

    pylab.rcParams.update(params)
    pylab.figure(1)
    pylab.clf()

    #
    # Make the plot
    #

    pylab.scatter(x, y, c="0.6", marker="o", lw=0, s=2.0, edgecolor="none")
    pylab.hlines(0., xmin, xmax, colors="k", linestyles="dashed", lw=1.3)
    pylab.plot(xfit, yfit, "r-", lw=1.3)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    #
    # Adjust the plot
    #

    pylab.subplots_adjust(left=0.20, right=0.97, bottom=0.12, top=0.97, hspace=0.)

    pylab.axes().set_xlim((xmin, xmax))
    pylab.axes().set_ylim((ymin, ymax))

    #
    # Save the plot
    #

    pylab.savefig(plot_file)
    pylab.close()
    
# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------

#
# Read the points and best-fit line
#

x, y, m, b = ReadData(data_file)

#
# Plot the data
#

PlotDMvsDNU(plot_file, x, y, m, b)

