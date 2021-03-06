#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:16:58 2017
Input values for modelling and plotting Geoelectric fields using
@author: joan
"""
import numpy as np
import os

###############################################################################
# 1) Input variables
###############################################################################

run_parallel_version = False
parallel_pool_size = 5

# 1.1) Main Paths

# main_path = 'C:/Users/u11236/Desktop/Compute_GeoElectric_Fields/'
main_path = os.path.dirname(os.path.abspath(__file__))


# 1.1.2) Folder with input magnetic time series Joan's
# storm = '22-23_06_2015'

# AWAGS magnetic storm data
storm = '01-02_12_1989'

# VIC64R Auslamp 21-22-12-2016
# storm = '21-22-12-2016'


# 1.2) Modelling approaches based on Campanya et al., 2018
# 1.2.1) (1) for Approach #1 and (2) for Approach #2
mode = 1


###############################################################################
# Select the fields depending if we are with Approach #1 or Approach #2
###############################################################################
fname = ""
if mode == 1:
    fname = "SECS_"
if mode == 2:
    fname = "dSECS_"

# 1.2.2) if 1 will not compute interpolated magnetic fields
avoid_comp_secs = 0

# 1.3) Periods of interest
# 1.3.1) maximum period to analyse (seconds)
# hi = 10 ** 4.2   # Joan's
hi = 10 ** 4.0   # VIC64R
# hi = 10 ** 10.0

# 1.3.2) minimum period to analyse (seconds)
low = 10 ** 2
# no bandpass filter
# low = 0
# 1.4) Sampling rate (seconds)
samp = 60.

earthrad, ionorad = 6371000.0, 6471000.0  # 6371000.0, 6481000.0
Samp_P_day = int(86400 / samp)  # Number of samples pr day !!!!


# 1.5 Time series properties
# 1.5.1) Starting point for the analysis
mint = 1000

# 1.5.2) End point for the analysis
maxt = -800

# 1.6) Area of interest for SECS interpolation
# Define grid: uniform in lat and long Ireland
#secswest, secseast, secssouth, secsnorth = -15, 15, 43, 65

# SA
#secswest, secseast, secssouth, secsnorth = 122, 145, -43, -20

# SA extended, only do 2 telluric sites, memory error, but 365 SECS is ok. input 25 sites
secswest, secseast, secssouth, secsnorth = 122, 155, -43, -20

# reduce the inpt sites to 23 (ROM, CND CNB excluded) still memeory error
#secswest, secseast, secssouth, secsnorth = 122, 151, -43, -20


# SA and VIC, memory error 25 sites, not working for this case
#secswest, secseast, secssouth, secsnorth = 120, 160, -45, -17

# 1.7) Only for Approach #2
# 1.7.1) Ref. magnetic site - regional signal (Approach #2)
reg_ref = 'CLF'

# 1.7.2) Ref. magnetic sites to compute e_fields (Approach #2)

rmf = ['HAD', 'BIR']

###############################################################################
# 2) Additional inputs
###############################################################################

# No need to modify them if following the suggested structure and parameters
# from Campanya et al., 2018

# 2.1) Errors
# 2.1.1) Error floor for the Non-plane wave approximation
e_nvpwa = 1 / np.sqrt(10.0 ** (10.0 / 10.0))

# 2.1.2) Error floor for the MT and quasi-MT tensor relationships
ef_tf = 10e-2

# 2.1.3)  Error floor for H tensor relationship
ef_h = 2e-2

# 2.2) Statistics for error propagation
# stat = 1000
stat = 500

# 2.3) Paths of interest
# 2.3.1) Folder with data from a particular geomagnetic storm
data_path = os.path.join(main_path, 'in', 'data', storm)

# 2.3.2) Folder with magnetic fields time series
mag_path = os.path.join( data_path, 'B')

# 2.3.3) Folder with electric field time series
e_path = os.path.join(data_path, 'E')

# 2.3.4) Folder with input parameters
in_path = os.path.join(main_path, 'in')

# 2.3.5) Folder with output parameters
out_path = os.path.join(main_path, 'out')

# 2.3.6) Folder with electromagnetic transfer functions
tf_path = os.path.join(main_path, 'in', 'data', 'TF')

# 2.3.7) Folder with inputs - outputs for SECS interpolation
secs_path = os.path.join(out_path, 'SECS')

# 2.4) Files with sites of interest
# 2.4.1) Magnetic observatories

obs_f = 'Observatories.dat'

# 2.4.2) Sites where to calculate the electric fields
sit_f = 'sites_interest.dat'

# 2.5) Time series properties
# 2.5.1) length of the time series

len_ts = sorted(os.listdir(mag_path))
len_path = os.path.join(mag_path, str(len_ts[0]))
fileHandler = open(len_path, 'r')
len_val = np.loadtxt(fileHandler)
fileHandler.close()

length = len_val.shape[0]
