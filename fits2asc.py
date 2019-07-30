#!/usr/bin/env python3

# PLots lightcurve from ucam/uspec files
# User needs to add file names
# on command line when invoking.
# Copied from Athena July 2019 by D Sahman
# ===============================

import pyfits
import numpy
import pylab
import sys

def getData(file):
    data =pyfits.getdata(file)
    mjd  =data['MJD']
    c1   =data['Counts_1']
    c2   =data['Counts_2']
    e1   =data['Sigma_1']
    e2   =data['Sigma_2']
    eflag1 =data['Eflag_1'] == 0
    eflag2 =data['Eflag_2'] == 0
    flux =c1/c2
    fluxerr =(((e1/c1)**2+(e2/c2)**2)**0.5)*flux
    check =eflag1 & eflag2
    return (mjd[check], flux[check], fluxerr[check])

files = sys.argv[1:]
mjdlist = []
fluxlist = []
fluxerrlist = []

for file in files:
    mjd, flux, fluxerr = getData(file)
    mjdlist.append(mjd)
    fluxlist.append(flux)
    fluxerrlist.append(fluxerr)

mjd_all = numpy.concatenate(mjdlist)
flux_all = numpy.concatenate(fluxlist)
fluxerr_all = numpy.concatenate(fluxerrlist)

pylab.plot(mjd_all, flux_all, 'r.')
pylab.errorbar(mjd_all, flux_all, yerr=fluxerr_all, color="red", fmt='.')
pylab.xlabel("MJD")
pylab.ylabel("Flux")
pylab.show()