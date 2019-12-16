#!/usr/bin/env python3

# Program to fit gaussian to light curve- first attempt!
# D Sahman Nov 2018
#=======================================================

import numpy as np
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
print ('Please enter filename: ')
filename=input()
data = np.loadtxt(filename)

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

# print ('indx6  ', indx6)
# print ('indx8  ', indx8)
# print ('indx9  ', indx9)
# print ('indx10  ', indx10)
# print ('indx11  ', indx11)

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

# Plot all three CCDs
# ===================

fig, axarr = plt.subplots(3, sharex=True)
mjd1 = data1[:,1:2]
counts1 = data1[:,14:15]
mjd2 = data2[:,1:2]
counts2 = data2[:,14:15]
mjd4 = data4[:,1:2]
counts4 = data4[:,14:15]
axarr[0].plot(mjd1,counts1, linestyle='-', color='red')
axarr[0].set_ylabel('counts')
axarr[1].plot(mjd2,counts2, linestyle='-', color='green')
axarr[1].set_ylabel('counts')
axarr[2].plot(mjd4,counts4, linestyle='-', color='blue')
axarr[2].set_ylabel('counts')
axarr[2].set_xlabel('MJD')
plt.show()
