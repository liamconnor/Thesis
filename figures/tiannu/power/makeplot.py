import numpy
import matplotlib
matplotlib.use("PDF")
import pylab

# ----------------------------------------------------------------------
# PARAMETERS
# ----------------------------------------------------------------------

# I/O
file_dmdm = "data/ngp_density_delta0.010_dmdm.dat"
file_nunu = "data/ngp_density_delta0.010_nunu.dat"
file_haha = "data/halopower_g0g1.dat"
file_cold = "data/ngp_density_delta0.010_dmdm_cold.dat"
file_tran = "data/ith2_mnu0p05_z0p05_tk.dat"
plot_file = "../fig2.pdf"

# ----------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------

def ReadPower(data_file):
    """
    Returns k and Delta^2(k).
    """

    # Read the data
    k, d, g = numpy.genfromtxt(data_file, dtype="float32", usecols=(1,2,3), unpack=True)

    # Delete any trouble points
    ib = numpy.where((k <= 0.) | (d <= 0.))
    k  = numpy.delete(k, ib)
    d  = numpy.delete(d, ib)
    g  = numpy.delete(g, ib)   
 
    return k, d, g

def ComputeLinearPower(tfile):
    """
    Compute the linear power spectra for CDM and neutrinos.
    """

    # Cosmological parameters from the simulation

    As = 2.109524866209489e-09
    h  = 0.67
    ns = 0.96
    k0 = 0.05

    # Read the transfer function

    k, tdm, tnu = numpy.genfromtxt(tfile, dtype="float64", usecols=(0,1,5), unpack=True)

    # Convert to dimensionless power spectra

    pdm = 2.*numpy.pi**2 * As * h**ns * k0 * (k/k0)**ns * tdm**2
    pnu = 2.*numpy.pi**2 * As * h**ns * k0 * (k/k0)**ns * tnu**2

    pdm = k**3 * pdm / 2. / numpy.pi**2
    pnu = k**3 * pnu / 2. / numpy.pi**2

    pdm *= h**3
    pnu *= h**3

    return k, pdm, pnu    

def PlotPower(plot_file, xdm, ydm, xnu, ynu, xha, yha, xl, cl, nl, xcc, ycc, gnu):
    """
    Plot power spectra and compare to linear/non-linear models.
    """

    assert numpy.array_equal(xdm, xcc)

    #
    # Compute poisson noise
    #

    L  = 1200.
    nc = 13824.
    ni = (L/nc)**3.
    pn = xnu**3 * ni / 2.0 / numpy.pi**2

    #
    # Get total matter power spectra
    #

    omegac = 0.32
    omegan = 0.05/93.14/0.67**2
    omegam = omegac + omegan

    dmatcd = ycc.copy()
    dmatnu = omegac/omegam*ydm + omegan/omegam*((ynu-pn)/gnu)
    pdiff  = dmatnu/dmatcd

    #
    # Plot parameters
    #

    xlabel  = r"$k\ [h/{\rm Mpc}]$"
    ylabel  = r"$\Delta^2(k)$"
    ylabel2 = r"$P_{mm}^{\rm TN}/P_{mm}^{\rm TC}$"

    xmin  = 1.e-3
    xmax  = 1.e2
    ymin  = 2.e-6
    ymax  = 1.e3  
    ymin2 = 0.988
    ymax2 = 0.997 

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

    ax1 = pylab.subplot(211)

    pylab.loglog(xl, cl, "k:", label=r"$cc\ ({\rm Linear})$")
    pylab.loglog(xdm, ydm, "k-", label=r"$cc$")
    pylab.loglog(xl, nl, "r:", label=r"$\nu\nu\ ({\rm Linear})$")
    pylab.loglog(xnu, (ynu-pn)/gnu, "r-", label=r"$\nu\nu$")
    pylab.ylabel(ylabel)

    lg = pylab.legend(loc="upper left", fancybox=True, numpoints=1)
    lg.draw_frame(False)

    ax2 = pylab.subplot(212)

    pylab.semilogx(xcc, pdiff)

    pylab.ylabel(ylabel2)
    pylab.xlabel(xlabel)

    #
    # Adjust the plot
    #

    ax1.set_xlim((xmin, xmax))
    ax2.set_xlim((xmin, xmax))
    ax1.set_ylim((ymin, ymax))
    ax2.set_ylim((ymin2, ymax2))

    pylab.setp(ax1.get_xticklabels(), visible=False)

    pylab.subplots_adjust(left=0.18, right=0.97, bottom=0.12, top=0.97, hspace=0.)

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

xdm, ydm, gdm = ReadPower(file_dmdm)
xnu, ynu, gnu = ReadPower(file_nunu)
xha, yha, gha = ReadPower(file_haha)
xcc, ycc, gcc = ReadPower(file_cold)

#
# Compute linear power spectra
#

xl, cl, nl = ComputeLinearPower(file_tran)

#
# Plot the data
#

PlotPower(plot_file, xdm, ydm, xnu, ynu, xha, yha, xl, cl, nl, xcc, ycc, gnu)


