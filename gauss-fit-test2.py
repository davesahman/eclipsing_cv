#!/usr/bin/env python3

# Program to fit gaussian to light curve- first attempt!
# D Sahman Nov 2018
#=======================================================

import numpy as np
from numpy import exp, linspace,random
import scipy
from scipy.optimize import curve_fit
import glob
import matplotlib
from matplotlib import pyplot as plt
from hipercam.hlog import Hlog
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import gaussian_fwhm_to_sigma

# Set up gaussian model
def gaussian(x, amp, cen, wid):
    return amp * exp (-(x-cen)**2/(2*wid**2))


def est_fwhm(x, y):
    """
    Estimate FWHM of a Gaussian
    """
    half_max = 0.5*y.max()
    within_halfmax = y > half_max
    x_within_halfmax = x[within_halfmax]
    return x_within_halfmax.max() - x_within_halfmax.min()
    

def get_lc(fname, ccdnam, targ_ap, comp_aps):
    """
    Extract relative lightcurve from ULTRACAM logfile

    Calculates the relative lightcurve from an 
    ULTRACAM logfile, taking into account error
    codes and returning valid data only.

    Parameters
    ----------
    fname: string
        Name of logfile
    ccdnam: string
        CCD name ('1', '2' or '3' for ULTRACAM)
    targ_ap: string
        The aperture of the target
    comp_aps: string or list of strings
        The aperture (or apertures of the comparison stars)

    Returns
    -------
    t, y, ye: np.ndarray
        MJD, relative counts and error
    """
    lf = Hlog.from_ulog(fname)

    ts_targ = lf.tseries(ccdnam, targ_ap)
    # make sure that comp_aps is a list
    comp_aps = list(comp_aps)
    # make comparison time series from first comp ap
    ts_comp = lf.tseries(ccdnam, comp_aps[0])
    # sum any other comparisons onto same tseries
    for comp_ap in comp_aps[1:]:
        ts_comp += lf.tseries(ccdnam, comp_ap)
    # make relative lightcurve
    rel_ts = ts_targ/ts_comp

    # now mask out error codes 6, 8, 9, 10, 11
    mask = ((rel_ts.mask != 6) & (rel_ts.mask != 8) &
            (rel_ts.mask != 9) & (rel_ts.mask != 10) &
            (rel_ts.mask != 11))

    # return masked data
    return rel_ts.t[mask], rel_ts.y[mask], rel_ts.ye[mask]

    
# set print option to print to 12 decimal places
# ==============================================
np.set_printoptions(precision=12)

# TODO: get fname, ccdname, etc from command line arguments
fname = 'oycar_1.log'
ccdnam = '2'
targ_ap = '1'
comp_aps = '2'
mjd, counts, counts_err = get_lc(fname, ccdnam,targ_ap, comp_aps)

# Calculate derivative of the signal
# ==================================
smooth_counts = convolve(counts, Box1DKernel(10))
df = smooth_counts[1:] - smooth_counts[:-1]
dt = mjd[1:] - mjd[:-1]
dfdt = df/dt
mid_mjd = 0.5 * mjd[:-1] + 0.5 * mjd[1:]

# trim the data to exclude noisy end points
mid_mjd = mid_mjd[9:-11]
dfdt = dfdt[9:-11]

# subtract off integer part of MJD
# otherwise floating point errs in minimisation
tfloor = int(mid_mjd.min())
mid_mjd -= tfloor

# Smooth the derivative
dfdt = convolve(dfdt, Box1DKernel(10))

# Fit two Gaussians to the max and min values of dfdt_ts
#-------------------------------------------------------
# find min and max values of dfdt_ts
min_arg = np.argmin(dfdt)
t_min = mid_mjd[min_arg]
max_arg = np.argmax(dfdt)
t_max = mid_mjd[max_arg]

# get two gaussians around peak to fit
t1 = mid_mjd[min_arg-20:min_arg+20]
t2 = mid_mjd[max_arg-20:max_arg+20]
g1 = np.fabs(dfdt[min_arg-20:min_arg+20])
g2 = np.fabs(dfdt[max_arg-20:max_arg+20])

# First gaussian fit
# Initial guesses
amp = -1 * dfdt[min_arg]
cen = t_min
wid = est_fwhm(t1, g1) * gaussian_fwhm_to_sigma
init_vals = [amp ,cen, wid]

best_vals, covar = curve_fit(gaussian, t1, g1, p0=init_vals)
print("Start of eclipse")
print("================")
print("init_vals = ", init_vals)
print ("best vals = ", best_vals)
print(" amp = %.3f +/- %.3f" % (best_vals[0], np.sqrt(covar[0,0])))
print(" cen = %.3f +/- %.3f" % (best_vals[1], np.sqrt(covar[1,1])))
print(" wid = %.3f +/- %.3f" % (best_vals[2], np.sqrt(covar[2,2])))

# Second Gaussain fit

# Initial guesses
amp2 = dfdt[max_arg]
cen2 = t_max
wid2 = est_fwhm(t2, g2) * gaussian_fwhm_to_sigma
print("cen2 = ", cen2)
init_vals2 = [amp2 ,cen2, wid2]

best_vals2, covar2 = curve_fit(gaussian, t2, g2, p0=init_vals2)
print("End of eclipse")
print("==============")
print("init_vals2 = ", init_vals2)
print ("best vals2 = ", best_vals2)
print(" amp2 = %.3f +/- %.3f" % (best_vals2[0], np.sqrt(covar2[0,0])))
print(" cen2 = %.3f +/- %.3f" % (best_vals2[1], np.sqrt(covar2[1,1])))
print(" wid2 = %.3f +/- %.3f" % (best_vals2[2], np.sqrt(covar2[2,2])))

fig, axs = plt.subplots(2,2,figsize=(10,7))
axs[0,0].plot(mjd, counts, linestyle='-', color='red')
ax = axs[0, 0]
# ax.set_xlabel("Time (days)")
ax.set_title(" Light curve - Red CCD")
axs[0, 1].plot(mid_mjd, dfdt, linestyle='-', color ='red')
ax = axs[0, 1]
# ax.set_xlabel("Time (days)")
ax.set_title("Smoothed Derivative")
axs[1, 0].plot(t1, gaussian(t1, *best_vals), label="fit")
axs[1, 0].plot(t1, g1, label="data")
axs[1, 0].legend()
ax = axs[1, 0]
ax.set_xlabel("Time (MJD - {})".format(tfloor))
ax.set_title("First Gaussian Fit")
axs[1, 1].plot(t2, gaussian(t2, *best_vals2), label="fit")
axs[1, 1].plot(t2, g2, label="data")
axs[1, 1].legend()
ax = axs[1, 1]
ax.set_xlabel("Time (days)")
ax.set_title("Second Gaussian Fit")
plt.show()

# Calculate Eclipse Time
mid_cen = tfloor + (best_vals[1]+best_vals2[1])/2
#mid_cen_err = np.sqrt(covar[1,1]) if np.sqrt(covar[1,1]) > np.sqrt(covar2[1,1]) else np.sqrt(covar2[1,1])
mid_cen_err = max(np.sqrt(covar[1,1]),np.sqrt(covar2[1,1]))  
print("Eclipse Time = ",mid_cen, "+/-", mid_cen_err)

