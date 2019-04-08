# #!/usr/bin/env python2
# # -*- coding: utf-8 -*-
# """
# # Code for computing the Spherical Elementary Current Systems (SECS) to
# # interpolate magnetic field variations.
# # Original Author: C.D.Beggan (adapted from AJmcK) [ciar@bgs.ac.uk]
# # Ported from MATLAB to Python by Sean Blake (Trinity College Dublin) [blakese@tcd.ie]
# # Modified by Joan Campanya (Trinity College Dublin) [joan.campanya@tcd.ie]
# # Date: 31-Aug-2018
# """
#
#
# import matplotlib.dates as mdates
# from multiprocessing import Pool
#
# from constants import in_path, samp, length, secswest, secssouth, secseast, secsnorth, sit_f, out_path, storm, mode
#
# import numpy as np
# from scipy.interpolate import griddata
# from scr.secs_pre import secsmatrix_XYonly
# from scr.EM_modelling import lats_secs, lons_secs
# # exec(open('secs_pre.py').read())
#
#
# in_path = in_path
# Samp = samp
# Leng = length
#  ################################################################################
# # DATA INPUT
# ################################################################################
# # Declare the constants
# earthrad, ionorad = 6371000.0, 6471000.0    #6371000.0, 6481000.0
# Samp_P_day =int( 86400 / Samp) # Number of samples pr day !!!!
# # Define grid: uniform in lat and long
# #secswest, secseast, secssouth, secsnorth = -13, 10, 48, 63
# #secswest, secseast, secssouth, secsnorth = -15, 15, 40, 67  #!!!!!
#
# lonmin, latmin = secswest, secssouth
#
# dlon = 0.5 # spacing !!!!
#
# lonpts = int((secseast - secswest)*(1/dlon)+1)
# latpts = int((secsnorth - secssouth)*(1/dlon)+1)
#
#
# ncol = latpts * lonpts * 2
# numb_of_days =int( Leng * Samp / 86400)   # Days recording !!!!
#
# print ("Forming Grid:", lonpts, "x", latpts, "[lon x lat grid]")
# lon, lat = np.meshgrid(np.arange(lonmin, lonmin+dlon*lonpts-dlon + 0.1, dlon),
#     np.arange(latmin, latmin + dlon*latpts-dlon+0.1, dlon))
#
# # Folder which contains all of your inputs
# datafolder =  str(out_path)
#
# # site data (name, lat, lon) are in file: .../datafolder/sites.txt
# sitesfile = str(in_path) + "Observatories.dat"
# sitenames = np.loadtxt(sitesfile, usecols = (0,), unpack = True, skiprows =0, dtype = str)
# sitelat, sitelon = np.loadtxt(sitesfile, usecols = (1,2), unpack = True, skiprows =0)
#
# nsites = len(sitelat)
# allBmeas = np.zeros((numb_of_days * Samp_P_day, 2*nsites))
# allBxmeas, allBymeas = [], []
#
# # TODO remove dependancy
# if mode == 1:
#     fname = "SECS_"
# if mode == 2:
#     fname = "dSECS_"
#
# # Load data from files- varying Bx and By values, one per minute
# for index, value in enumerate(sitenames):
#     filename = datafolder + "SECS/" + str(fname) + value + storm + ".dat"
#
#
#     X ,Y = np.loadtxt(filename, usecols = (0,1), unpack = True, skiprows = 0)
#
#     X = X[0:numb_of_days*Samp_P_day]
#     Y = Y[0:numb_of_days*Samp_P_day]
#     allBmeas[:,index*2] = X
#     allBmeas[:,index*2 + 1] = Y
#
# ################################################################################
# # Construct the SECS matrix
# Amatrix_ext = np.zeros((nsites*2, int(latpts*lonpts)))
#
# for n in np.arange(0, nsites, 1):
#     nm = (2*n)
#     Tex, Tey = secsmatrix_XYonly(latpts, lonpts, ncol, sitelat[n], sitelon[n],
#       lat, lon, earthrad, ionorad)
#
#     Amatrix_ext[nm] = Tex
#     Amatrix_ext[nm+1] = Tey
#
#     testx = Tex.flatten()
#     testy = Tey.flatten()
#
#     for i in testx:
#         if np.isnan(i) == True:
#             print (n)
#             Tex_backup = Tex
#     for j in testx:
#         if np.isnan(j) == True:
#             Tey_backup = Tey
#             print (n)
#
# # make sure there are no Nan's in the data
# Amatrix_ext = np.nan_to_num(Amatrix_ext)
#
# # Perform a Singular Value Decomposition of the SECS matrix
# Ue, We, Ve = np.linalg.svd(Amatrix_ext, full_matrices = False)
# Ve = Ve.T
#
# svdthresh = We[0]/100.0
# for index, value  in enumerate(We):
#     if value < svdthresh:
#         We[index] = 0
#
# truncate = len(We)
# for i in We:
#     if i < svdthresh:
#         truncate = truncate - i
#
#
# print ("Calculating B throughout the grid")
# TmatrixX = np.zeros((lonpts*latpts, lonpts*latpts))
# TmatrixY = np.zeros((latpts*lonpts, latpts*lonpts))
#
# for m in np.arange(0, int(lonpts), 1):
#     grdlat = lat[:,m]
#     grdlon = lon[:,m]
#     zzz = np.arange(int(m*latpts), int((m*latpts) + latpts), 1)
#
#     Tx, Ty = secsmatrix_XYonly(latpts, lonpts, ncol, grdlat, grdlon, lat,
#         lon, earthrad, ionorad)
#
#     TmatrixX[:, zzz] = Tx
#     TmatrixY[:, zzz] = Ty
#
# TmatrixX = TmatrixX * -1
# TmatrixY = TmatrixY * -1
#
#
# print ("Calculating B-Fields")
#
# # Arbitrarily calculates every 7th minute over the 2 day period
# #for minute in np.arange(0, numb_of_days* Samp_P_day, 1):
# A = np.matrix(Ve[:, range(0, truncate)])
# B = np.matrix(np.diag(1/We[:truncate]))
# C = np.matrix(Ue[:, range(0, truncate)].T)
# Tmatrix = np.vstack((TmatrixX, TmatrixY))
# Tmatrix = np.nan_to_num(Tmatrix)
#
# #for minute in range(0, numb_of_days* Samp_P_day, 1):
# def single_calc(minute):
#     Bmeas = allBmeas[minute]
#
#     D = np.matrix(Bmeas).T
#     Ecurr = A*B*C*D
#     Bfield_XYonly = Tmatrix*Ecurr
#
#     Bx = Bfield_XYonly[0:latpts*lonpts]
#     Bx = np.reshape(Bx, (lonpts, latpts), order = 'C')
#     Bx = Bx.T
#     Bx = np.array(Bx)
#     Bx_1D = Bx.flatten()
#
#     By = Bfield_XYonly[latpts*lonpts:]
#     By = np.reshape(By, (lonpts, latpts), order = 'C')
#     By = By.T
#     By = np.array(By)
#     By_1D = By.flatten()
#
#     # interpolating the data
#     lon_1D = lon.flatten()
#     lat_1D = lat.flatten()
#
#     points = np.column_stack([lon_1D, lat_1D])
#     otherpoints = []
#
#     # 3.3) Open the observatories datafile and select the coordinates of the sites
#     # # TODO REMOVE THIS DEPANDANCY
#     # f = open(in_path + str(sit_f), 'r')
#     # lats_secs, lons_secs = np.loadtxt(f,
#     #                                   usecols=(1, 2),
#     #                                   unpack=True,
#     #                                   skiprows=0)
#     # f.close()
#
#     otherpoints = np.column_stack([lons_secs, lats_secs])
#
#     bx_interp = griddata(points, Bx_1D, (otherpoints), method = 'cubic')
#     by_interp = griddata(points, By_1D, (otherpoints), method = 'cubic')
#
# #    filename = out_path + 'SECS/' + "%08d" % (minute+1) + ".txt"
# #    outt = (bx_interp, by_interp)
# #    np.savetxt(filename, outt)
#
#     print (minute)
#     return(bx_interp, by_interp)
#
#
# # # TODO REMOVE DEPENDANCY
# # f = open(in_path + str(sit_f),'r')
# #
# # lats_secs, lons_secs = np.loadtxt(f,
# #                                 usecols=(1, 2),
# #                                 unpack=True,
# #                                 skiprows = 0)
# # f.close()
#
# output_secs = (np.zeros([numb_of_days* Samp_P_day,2,len(lats_secs) ]))
#
# #######################
# # Non-parallel_version
# ######################
# print("Total runs",numb_of_days* Samp_P_day)
#
# for i in range(0, numb_of_days* Samp_P_day, 1):
#    bx_interp, by_interp = single_calc(i)
#    output_secs[i,:,:] = np.array([bx_interp, by_interp])
# # !!!!!!if using non_parallel comment the lines below !!!!!!!!
#
# # Parallel version
# # def output_secs_par(i):
# #     bx_interp, by_interp = single_calc(i)
# #     output = np.array([bx_interp, by_interp])
# #     return(output)
# #
# #
# # with Pool(60) as p:
# #      output_secs = p.map(output_secs_par, range(0, numb_of_days* Samp_P_day, 1))