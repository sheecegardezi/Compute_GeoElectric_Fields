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


# 1.1) Main Paths
main_path = 'C:/Users/u12089/Desktop/Compute_GeoElectric_Fields/'
# 1.1.2) Folder with input magnetic time series
storm = '22-23_06_2015'

# 1.2) Modelling approaches based on Campanya et al., 2018
# 1.2.1) (1) for Approach #1 and (2) for Approach #2
mode = 1

# 1.2.2) if 1 will not compute interpolated magnetic fields
avoid_comp_secs = 0

# 1.3) Periods of interest
# 1.3.1) maximum period to analyse (seconds)
hi = 10 ** 4.2

# 1.3.2) minimum period to analyse (seconds)
low = 10 ** 2

# 1.4) Sampling rate (seconds)
samp = 60.

# 1.5 Time series properties
# 1.5.1) Starting point for the analysis
mint = 1000

# 1.5.2) End point for the analysis
maxt = -800

# 1.6) Area of interest for SECS interpolation
secswest, secseast, secssouth, secsnorth = -15, 15, 43, 65

# 1.7) Only for Approach #2
# 1.7.1) Ref. magnetic site - regional signal (Approach #2)
reg_ref = 'CLF'

# 1.7.2) Ref. magnetic sites to compute e_fields (Approach #2)
# rmf = ['HAD', 'BIR']
rmf = ['HAD', 'BIR']
###############################################################################
# 2) Additional inputs
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
stat = 1000

# 2.3) Paths of interest
# 2.3.1) Folder with data from a particluar geomagnetic storm
data_path = main_path + 'in/data/' + storm + '/'

# 2.3.2) Folder with magnetic fields time series
mag_path = data_path + 'B/'

# 2.3.3) Folder with electric field time series
e_path = data_path + 'E/'

# 2.3.4) Folder with input parameters
in_path = main_path + 'in/'

# 2.3.5) Folder with output parameters
out_path = main_path + 'out/'

# 2.3.6) Folder with electromagnetic transfer functions
tf_path = main_path + 'in/data/TF/'

# 2.3.7) Folder with inputs - outputs for SECS interpolation
secs_path = out_path + 'SECS/'

# 2.4) Files with sites of interest
# 2.4.1) Magnetic observatories
obs_f = 'Observatories.dat'

# 2.4.2) Sites where to calculate the electric fields
sit_f = 'sites_interest.dat'

# 2.5) Time series properties
# 2.5.1) length of the time series
len_ts = sorted(os.listdir(mag_path))
len_path = mag_path + str(len_ts[0])
f = open(len_path, 'r')
len_val = np.loadtxt(f)
f.close()

length = len_val.shape[0]
