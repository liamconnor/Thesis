import numpy
import matplotlib
matplotlib.use("PDF")
import pylab

# ----------------------------------------------------------------------
# PARAMETERS
# ----------------------------------------------------------------------

# I/O
data_sim  = "data/Dipole_sim.txt"
data_lin  = "data/HongMing.txt"
plot_file = "../fig3.pdf"

# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

def ReadCorrelation(data_file): 
    """
    Returns r and xi(r).
    """

    # Read the data
    k, d = numpy.genfromtxt(data_file, dtype="float32", usecols=(0,1), unpack=True)

    # Delete any trouble points
    ib = numpy.where((k <= 0.) | (d <= 0.))
    k  = numpy.delete(k, ib)
    d  = numpy.delete(d, ib)
    
    return k, d

def ReadMBPTCorrelation(df):
    """
    Returns the dipole correlation function from Hong-Ming's paper.
    """

    r, d = numpy.genfromtxt(df, dtype="float32", usecols=(0,3), unpack=True)

    return r, d

def PlotXi(plot_file, xsim, ysim, xlin, ylin):
    """
    Plot the simulated and MBT dipole correlation functions. 
    """

    #
    # Plot parameters
    #

    xmin = 1.e-1
    xmax = 1.e2
    ymin = 1.e-4
    ymax = 1.e-1

    xlabel = r"$r\ [{\rm Mpc}/h]$"
    ylabel = r"$\xi_{c\nu1}$"

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

    pylab.loglog(xsim, ysim, "k-", label=r"${\rm TianSmall}$")
    pylab.loglog(xlin, ylin, "b:", label=r"${\rm Zhu\ et\ al.\ (2014b)}$")

    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    lg = pylab.legend(loc="upper left", fancybox=True, numpoints=1)
    lg.draw_frame(False)

    #
    # Adjust the plot
    #

    pylab.subplots_adjust(left=0.16, right=0.97, bottom=0.12, top=0.97, hspace=0.)

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
# Read the relevant data
#

xsim, ysim = ReadCorrelation(data_sim)
xlin, ylin = ReadMBPTCorrelation(data_lin)

#
# Plot the data
#

PlotXi(plot_file, xsim, ysim, xlin, ylin) 


