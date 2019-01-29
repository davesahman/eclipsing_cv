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

# set print option to print to 12 decimal places
# ==============================================
np.set_printoptions(precision=12)

# The lines below are alternative file load methods
# =================================================

# name, mjd, flag, expose, ccd, fwhm, beta, naper, x, y, xm, ym, exm, eym, counts, sigma, sky, nsky, nrej, worst, errorflag, naper2, x2, y2, xm2, ym2, exm2, eym2, counts2, sigma2, sky2, nsky2, nrej2, worst2, errorflag2 = np.loadtxt('oycar_1.log', unpack=True)
# data = np.loadtxt('oycar_1.log', usecols=[0,1,4,7,14,15,20,27,28,34])

# Load file
# =========

data = np.loadtxt('oycar_1.log')

a = data.shape
rows = a[0]
cols = a[1]

# Remove exposures with flag = 0
# ==============================
indx0 = []
for i in range(rows):
    if data[i,2] == 0:
        indx0.append(i)

data = np.delete(data,indx0, axis=0)

a = data.shape
rows = a[0]
cols = a[1]


# Remove error flags with values of 6, 8, 9, 10, 11
# =================================================
indx6 = []
indx8 = []
indx9 = []
indx10 = []
indx11 = []

for i in range(rows):
    if data[i,20] == 6 or data[i,34] == 6:
        indx6.append(i)
    if data[i,20] == 8 or data[i,34] == 8:
        indx8.append(i)
    if data[i,20] == 9 or data[i,34] == 9:
        indx9.append(i)      
    if data[i,20] == 10 or data[i,34] == 10:
        indx10.append(i)   
    if data[i,20] == 11 or data[i,34] == 11:
        indx11.append(i)  

data = np.delete(data,indx6, axis=0)
data = np.delete(data,indx8, axis=0)
data = np.delete(data,indx9, axis=0)
data = np.delete(data,indx10, axis=0)
data = np.delete(data,indx11, axis=0)

# Sort into types of ccd - I assume there are only 3
# ==================================================

indx1 = []
indx2 = []
indx3 = []
for i in range(rows):
    if data[i,4] == 1:
        indx1.append(i)
    if data[i,4] == 2:
        indx2.append(i)
    if data[i,4] == 3:
        indx3.append(i)

data1 = data[indx1]
data2 = data[indx2]
data3 = data[indx3]

# Note this dataset has the blue frames co-added into
# every third frame. There is no need to swap MJD from the second 
# frame into third frame, because the third frame already has a MJD 
# with mid-point of the three frames.
# So all I need to do is remove all frames with zero counts.
# ========================================================= 

a3 = data3.shape
rows3 = a3[0]
cols3 = a3[1]

indx4 = []
for i in range(rows3):
    if data3[i,20] == 12 or data3[i,34] == 12:
        continue
    indx4.append(i)

data4 = data3[indx4]    


# Calculate derivative of the signal
# ==================================
from astropy.convolution import convolve, Box1DKernel

mjd1 = data1[:,1:2]
mjd1 = np.squeeze(mjd1[:])
counts1 = data1[:,14:15]
counts1 = np.squeeze(counts1[:])
mjd2 = data2[:,1:2]
counts2 = data2[:,14:15]
mjd4 = data4[:,1:2]
counts4 = data4[:,14:15]
smooth_counts1 = convolve(counts1, Box1DKernel(10))
df1 = smooth_counts1[1:] - smooth_counts1[:-1]
dt1 = mjd1[1:] - mjd1[:-1]
dfdt = df1/dt1

a = np.zeros(shape=(df1.shape))
mid_mjd1 = 0.5 * mjd1[:-1] + 0.5 * mjd1[1:]

# trim the data to exclude noisy end points
mid_mjd1_t = mid_mjd1[9:-11]
dfdt_t = dfdt[9:-11]

# Smooth the derivative
dfdt_ts = convolve(dfdt_t, Box1DKernel(10))

# Fit two Gaussians to the max and min values of dfdt_ts
#-------------------------------------------------------

# find min and max values of dfdt_ts
min = np.amin(dfdt_ts,axis=None)
min_arg = np.argmin(dfdt_ts)
t_min = mid_mjd1_t[min_arg]
max1 = np.amax(dfdt_ts,axis=None)
max_arg = np.argmax(dfdt_ts)
t_max = mid_mjd1_t[max_arg]
t1 = mid_mjd1_t[min_arg-20:min_arg+20]
t1_1 = (t1-58144)*86400

t2 = mid_mjd1_t[max_arg-20:max_arg+20]
t2_1 = (t2-58144)*86400
g1 = np.fabs(dfdt_ts[min_arg-20:min_arg+20])
g2 = np.fabs(dfdt_ts[max_arg-20:max_arg+20])

# First gaussian fit

# Initial guesses
amp = -1 *min
cen = (t_min-58144)*86400
wid = 12.
init_vals = [amp ,cen, wid]

# Set up gaussian model

def gaussian(x, amp, cen, wid):
    return amp * exp (-(x-cen)**2/(2*wid**2)) 

best_vals, covar = curve_fit(gaussian, t1_1, g1, p0=init_vals)
print("Start of eclipse")
print("================")
print("init_vals = ", init_vals)
print ("best vals = ", best_vals)
print(" amp = %.3f +/- %.3f" % (best_vals[0], np.sqrt(covar[0,0])))
print(" cen = %.3f +/- %.3f" % (best_vals[1], np.sqrt(covar[1,1])))
print(" wid = %.3f +/- %.3f" % (best_vals[2], np.sqrt(covar[2,2])))

# Second Gaussain fit

# Initial guesses
amp2 = max1
cen2 = (t_max-58144)*86400
wid2 = 12
print("cen2 = ",cen2)
init_vals2 = [amp2 ,cen2, wid2]

best_vals2, covar2 = curve_fit(gaussian, t2_1, g2, p0=init_vals2)
print("End of eclipse")
print("==============")
print("init_vals2 = ", init_vals2)
print ("best vals2 = ", best_vals2)
print(" amp2 = %.3f +/- %.3f" % (best_vals2[0], np.sqrt(covar2[0,0])))
print(" cen2 = %.3f +/- %.3f" % (best_vals2[1], np.sqrt(covar2[1,1])))
print(" wid2 = %.3f +/- %.3f" % (best_vals2[2], np.sqrt(covar2[2,2])))

fig,axs = plt.subplots(2,2,figsize=(10,7))
mjd1 = data1[:,1:2]
counts1 = data1[:,14:15]
axs[0,0].plot(mjd1,counts1, linestyle='-', color='red')
ax =axs[0,0]
# ax.set_xlabel("Time (days)")
ax.set_title(" Light curve - Red CCD")
axs[0,1].plot(mid_mjd1_t,dfdt_ts, linestyle='-', color ='red')
ax =axs[0,1]
# ax.set_xlabel("Time (days)")
ax.set_title("Smoothed Derivative")
axs[1,0].plot(t1_1,gaussian(t1_1,*best_vals),label="fit")
axs[1,0].plot(t1_1,g1,label="data")
axs[1,0].legend()
ax =axs[1,0]
ax.set_xlabel("Time (secs)")
ax.set_title("First Gaussian Fit")
axs[1,1].plot(t2_1,gaussian(t2_1,*best_vals2),label="fit")
axs[1,1].plot(t2_1,g2,label="data")
axs[1,1].legend()
ax =axs[1,1]
ax.set_xlabel("Time (secs)")
ax.set_title("Second Gaussian Fit")
plt.show()

# Calculate Eclipse Time

mid_cen = (best_vals[1]+best_vals2[1])/2
#mid_cen_err = np.sqrt(covar[1,1]) if np.sqrt(covar[1,1]) > np.sqrt(covar2[1,1]) else np.sqrt(covar2[1,1])
mid_cen_err = max(np.sqrt(covar[1,1]),np.sqrt(covar2[1,1])) 
mid_cen_d = (mid_cen/86400)+58144 
mid_cen_err_d = (mid_cen_err/86400)
print("Eclipse Time = ",mid_cen_d, "+/-", mid_cen_err_d)

