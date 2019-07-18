#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
 Algorithm to modelled geoelectric fields at a particular site based on data
from magnetic observatories and geophysical parameters
 to measured data.

 Original author: Joan Campanya i Llovet (Trinity College Dublin) [joan.campanya@tcd.ie]
 08 - 2018, Dublin, Ireland.

 Inputs ("in" folder)
   1) Sites of interest where to measure E fields ('sites_interest.dat')
   2) Magnetometers to be used for for modelling E_fields('Observatories.dat')
   3) Electromagnetic transfer functions (TF folder)
   4) Magnetic time series at the magnetic observatories (B folder within
      the folder for each storm)

 Input values are specified in the input.py file in the "scr" folder.
 The program will read the input from this file and look for the corresponsing
 files in the "in" folder

 Outputs (in "out" folder)
   1) Modelled E fields at the sites of interest
   2) Results inputs and outputs for the Spherical elementary current systems
     (SECS folder)

Altgorithm related to Campanya et al. 2018 publication in Space Weather AGU Journal
"""
#######################################################################

import sys
from multiprocessing import Pool
from scipy.interpolate import griddata
from scr.secs_pre import secsmatrix_XYonly
from scr.functions_EM_modelling import *

import numpy as np
import matplotlib.pyplot as plt

from constants import *

plt.close('all')


def main():
    #######################################################################
    # 2) Read and Compute B fields
    #######################################################################
    print('Read and compute B fields')
    #######################################################################
    # 2.0) Approach #1
    #######################################################################
    if mode == 1:
        # 2.1) Read magnetic field data at the magnetic observatories
        mh_obs, h_obs = read_magnetics(in_path, str(obs_f), mag_path, secs_path, samp, hi, low, length, storm, 1)



    #######################################################################
    # 2.0) Approach #2
    #######################################################################
    elif mode == 2:
        # 2.1) Compute regional magnetics at the magnetic observatories
        comp, std_comp = compute_regional_b_field(in_path, str(obs_f), reg_ref, mag_path, tf_path, samp, hi, low, ef_h,
                                                  stat)

        # 2.2) Compute regional magnetics at the sites used to calculate E fields
        comp_rmf, std_comp_rmf = compute_regional_b_field(in_path, str(sit_f), reg_ref, mag_path, tf_path, samp, hi, low,
                                                          ef_h, stat)

        # 2.3) Calculate differences between modelled and measured magnetic fields
        # 2.3.1) At the magnetic Observatories
        # 2.3.1.1) Read magnetics
        pmh_obs, h_obs = read_magnetics(in_path, str(obs_f), mag_path, secs_path, samp, hi, low, length, storm, 0)

        # 2.3.1.2) Compute diferences
        mh_obs, dh_obs, std_dh_obs = difference_b_fields(pmh_obs, h_obs, comp, std_comp, secs_path, stat, storm)

    else:
        print("Error: Invalid mode selected.")
        sys.exit(0)

    #print("stop here if only want to calculating SECS.")
    #sys.exit(0)
    #######################################################################
    # 3) Interpolate magnetic fields, or differences, between observatories
    #######################################################################


    #######################################################################
    # 3.2) Open de observatories datafile and select the names of the sites
    #######################################################################
    with open(in_path + str(sit_f), 'r') as infile:
        name_secs = np.loadtxt(infile, usecols=(0,), unpack=True, skiprows=0, dtype=str)

    #######################################################################
    # 3.3) Open the observatories datafile and select the coordinates of the sites
    #######################################################################
    with open(in_path + str(sit_f), 'r') as infile:
        lats_secs, lons_secs = np.loadtxt(infile, usecols=(1, 2), unpack=True, skiprows=0)

    print(lats_secs)
    print(lons_secs)

    """
    # Code for computing the Spherical Elementary Current Systems (SECS) to
    # interpolate magnetic field variations.
    # Original Author: C.D.Beggan (adapted from AJmcK) [ciar@bgs.ac.uk]
    # Ported from MATLAB to Python by Sean Blake (Trinity College Dublin) [blakese@tcd.ie]
    # Modified by Joan Campanya (Trinity College Dublin) [joan.campanya@tcd.ie] 
    # Date: 31-Aug-2018
    """

    ################################################################################
    # 3.4 DATA INPUT  SECS routine
    # if avoid_comp_secs != 1:  # In case it was computed already
    ################################################################################
    if avoid_comp_secs != 1:  # In case it was computed already
        lonmin, latmin = secswest, secssouth

        dlon = 0.5  # spacing !!!!

        lonpts = int((secseast - secswest) * (1 / dlon) + 1)
        latpts = int((secsnorth - secssouth) * (1 / dlon) + 1)

        ncol = latpts * lonpts * 2
        numb_of_days = int(length * samp / 86400)  # Days recording !!!!

        print("Forming Grid:", lonpts, "x", latpts, "[lon x lat grid]")
        lon, lat = np.meshgrid(np.arange(lonmin, lonmin + dlon * lonpts - dlon + 0.1, dlon),
                               np.arange(latmin, latmin + dlon * latpts - dlon + 0.1, dlon))

        # Folder which contains all of your inputs
        datafolder = str(out_path)

        # site data (name, lat, lon) are in file: .../datafolder/sites.txt
        sitesfile = str(in_path) + "Observatories.dat"
        sitenames = np.loadtxt(sitesfile, usecols=(0,), unpack=True, skiprows=0, dtype=str)
        sitelat, sitelon = np.loadtxt(sitesfile, usecols=(1, 2), unpack=True, skiprows=0)

        nsites = len(sitelat)
        allBmeas = np.zeros((numb_of_days * Samp_P_day, 2 * nsites))
        allBxmeas, allBymeas = [], []

        # Load data from files- varying Bx and By values, one per minute
        for index, value in enumerate(sitenames):
            filename = datafolder + "SECS/" + str(fname) + value + storm + ".dat"

            X, Y = np.loadtxt(filename, usecols=(0, 1), unpack=True, skiprows=0)

            X = X[0:numb_of_days * Samp_P_day]
            Y = Y[0:numb_of_days * Samp_P_day]
            allBmeas[:, index * 2] = X
            allBmeas[:, index * 2 + 1] = Y

        ################################################################################
        # Construct the SECS matrix
        ################################################################################

        Amatrix_ext = np.zeros((nsites * 2, int(latpts * lonpts)))

        for n in np.arange(0, nsites, 1):
            nm = (2 * n)
            Tex, Tey = secsmatrix_XYonly(latpts, lonpts, ncol, sitelat[n], sitelon[n], lat, lon, earthrad, ionorad)

            Amatrix_ext[nm] = Tex
            Amatrix_ext[nm + 1] = Tey

            # testx = Tex.flatten()
            # testy = Tey.flatten()

            testx = Tex
            testy = Tey

            # for i in testx:
            #     if np.isnan(i) == True:
            #         print(n)
            #         Tex_backup = Tex
            # for j in testx:
            #     if np.isnan(j) == True:
            #         Tey_backup = Tey
            #         print(n)

        # make sure there are no Nan's in the data
        Amatrix_ext = np.nan_to_num(Amatrix_ext)

        # Perform a Singular Value Decomposition of the SECS matrix
        Ue, We, Ve = np.linalg.svd(Amatrix_ext, full_matrices=False)
        Ve = Ve.T

        svdthresh = We[0] / 100.0
        for index, value in enumerate(We):
            if value < svdthresh:
                We[index] = 0

        truncate = len(We)
        for i in We:
            if i < svdthresh:
                truncate = truncate - i

        print("Calculating B throughout the grid")
        TmatrixX = np.zeros((lonpts * latpts, lonpts * latpts))
        TmatrixY = np.zeros((latpts * lonpts, latpts * lonpts))

        for m in np.arange(0, int(lonpts), 1):
            grdlat = lat[:, m]
            grdlon = lon[:, m]
            zzz = np.arange(int(m * latpts), int((m * latpts) + latpts), 1)
            Tx, Ty = secsmatrix_XYonly(latpts, lonpts, ncol, grdlat, grdlon, lat, lon, earthrad, ionorad)

            TmatrixX[:, zzz] = Tx
            TmatrixY[:, zzz] = Ty

        TmatrixX = TmatrixX * -1
        TmatrixY = TmatrixY * -1

        print("Calculating B-Fields")

        # Arbitrarily calculates every 7th minute over the 2 day period
        # for minute in np.arange(0, numb_of_days* Samp_P_day, 1):
        A = np.matrix(Ve[:, range(0, truncate)])
        B = np.matrix(np.diag(1 / We[:truncate]))
        C = np.matrix(Ue[:, range(0, truncate)].T)
        Tmatrix = np.vstack((TmatrixX, TmatrixY))
        np.set_printoptions(precision=4)
        print(Tmatrix.size)
        Tmatrix = np.nan_to_num(Tmatrix)


        # for minute in range(0, numb_of_days* Samp_P_day, 1):
        def single_calc(minute):
            Bmeas = allBmeas[minute]

            D = np.matrix(Bmeas).T
            Ecurr = A * B * C * D
            Bfield_XYonly = Tmatrix * Ecurr

            Bx = Bfield_XYonly[0:latpts * lonpts]
            Bx = np.reshape(Bx, (lonpts, latpts), order='C')
            Bx = Bx.T
            Bx = np.array(Bx)
            Bx_1D = Bx.flatten()

            By = Bfield_XYonly[latpts * lonpts:]
            By = np.reshape(By, (lonpts, latpts), order='C')
            By = By.T
            By = np.array(By)
            By_1D = By.flatten()

            # interpolating the data
            lon_1D = lon.flatten()
            lat_1D = lat.flatten()

            points = np.column_stack([lon_1D, lat_1D])

            otherpoints = np.column_stack([lons_secs, lats_secs])

            bx_interp = griddata(points, Bx_1D, (otherpoints), method='cubic')
            by_interp = griddata(points, By_1D, (otherpoints), method='cubic')

            #    filename = out_path + 'SECS/' + "%08d" % (minute+1) + ".txt"
            #    outt = (bx_interp, by_interp)
            #    np.savetxt(filename, outt)

            print(minute)
            return (bx_interp, by_interp)


        output_secs = []
        if run_parallel_version:
            ################################################################################
            # Parallel version
            ################################################################################
            def output_secs_par(i):
                bx_interp, by_interp = single_calc(i)
                output = np.array([bx_interp, by_interp])
                return output


            with Pool(parallel_pool_size) as p:
                output_secs = p.map(output_secs_par, range(0, numb_of_days * Samp_P_day, 1))
        else:
            ################################################################################
            # Non-parallel_version
            ################################################################################
            output_secs = (np.zeros([numb_of_days * Samp_P_day, 2, len(lats_secs)]))
            for i in range(0, numb_of_days * Samp_P_day, 1):
                bx_interp, by_interp = single_calc(i)
                output_secs[i, :, :] = np.array([bx_interp, by_interp])


        # end if here if avoid_comp_secs != 1:  # In case it was computed already

        # 3.4) Compute SECS
        #if avoid_comp_secs != 1:  # In case it was computed already
        print('interpolate magnetic fields using SECS')
        # exec(open('SECS_interpolation.py').read())
        # obs_bx_secs = np.array(obs_bx_secs)
        # obs_by_secs = np.array(obs_by_secs)
        obs_bx_secs = np.array(output_secs)[:, 0, :]
        obs_by_secs = np.array(output_secs)[:, 1, :]
        # Write results from SECS
        for value, ip in enumerate(name_secs):
            secsx_id = open(secs_path + fname + storm + "magBx" + ip, 'w+')
            np.savetxt(secsx_id, obs_bx_secs[:, value], fmt=['%10.3f'])
            secsx_id.close()
            secsy_id = open(secs_path + fname + storm + "magBy" + ip, 'w+')
            np.savetxt(secsy_id, obs_by_secs[:, value], fmt=['%10.3f'])
            secsy_id.close()
    # stop here if only want to calculating SECS, LJW 2019-04-11
    print("stop here if only want to calculating SECS.")
    sys.exit(0)


    ################################################################################
    # 4) Define new variables for computing the e_fields
    ################################################################################
    print('Defining new variables for computing electric fields')
    print(stat)
    # TODO investigate breaking into a loop
    # 4.1) Read sites where we want to compute the E fields
    e_site, e_lat, e_lon = read_co(in_path + str(sit_f))
    """
    print(e_lat)
    print(e_lon)
    e_site_1 = [1]
    e_lat_1 = [1]

    for i in range(len(e_site)):
        e_site_1[0] = e_site[i]
        e_lat_1[0] = e_lat[i]
        print (i,e_site_1, e_lat_1, len(e_site_1))
        // for loop to compute e field at one site at a time to avoid memeory shortage error - Liejun 2019
        // to do here 
    """

    # 4.2) Define size of several variables needed to compute the E fields
    size_a = ([length, len(e_site), 2])
    size_b = [length, len(e_site), len(rmf), 2]

    secs_e = np.zeros(size_a)
    stdsecs_e = np.zeros(size_a)
    std_error = np.zeros(size_a)
    av_e_fields = np.zeros(size_a)
    std_reg_e = np.zeros(size_a)
    av_reg_e = np.zeros(size_a)
    std_secs_e = np.zeros(size_a)
    av_secs_e = np.zeros(size_a)
    std_loc_e = np.zeros(size_a)
    av_loc_e = np.zeros(size_a)

    c_reg_e = np.zeros(size_b)
    c_std_reg_e = np.zeros(size_b)
    d_rmf_e = np.zeros(size_b)
    std_rmf_e = np.zeros(size_b)

    ################################################################################
    # 4.3) Define variables that will be used for error propagation, including stat
    ################################################################################
    secs_e_fields = np.zeros([length, len(e_site), 2, stat])
    reg_e_fields = np.zeros([length, len(e_site), 2, len(rmf), stat])

    print(stat)
    if mode == 1:
        ################################################################################
        # Approach #1
        ################################################################################
        loc_emt = np.zeros(size_a)
        std_loc_emt = np.zeros(size_a)
        e_fields = np.zeros([length, len(e_site), 2, stat])
        comp_e_fields = np.zeros([length, len(e_site), 2, stat])

    elif mode == 2:
        ################################################################################
        # Approach #2
        ################################################################################
        loc_emt = np.zeros(size_a)
        std_loc_emt = np.zeros(size_a)
        e_fields = np.zeros([length, len(e_site), 2, len(rmf), stat])

        comp_e_fields = np.zeros([length, len(e_site), 2, stat])

        loc_e_fields = np.zeros([length, len(e_site), 2, stat])

        reg_e_fields = np.zeros([length, len(e_site), 2, len(rmf), stat])
    else:
        print("Error: Invalid mode selected.")
        sys.exit(0)

    #######################################################################
    # 5) Compute electric fields at the sites of interest
    ################################################################################
    print('Computing electric fields')
    for v1, ip1 in enumerate(e_site):  # Sites where to calculate the e fields

        ################################################################################
        # 5.1) Read magnetics computed from the interpolation approach
        ################################################################################
        rmfb = read_rmfb(secs_path + fname, ip1, storm)

        ################################################################################
        # 5.2) Compute the error associated with the interpolation of the magnetic
        # fields (Based on Figure 7, Campanya et al. 2018 )
        ################################################################################
        error_bf = error_secs_interpolation(e_lat[v1], e_lon[v1], in_path, mode, obs_f)

        ################################################################################
        # 5.3) Compute the electric fields unsing results from SECS
        ################################################################################
        ex_secs, ey_secs, std_ex_secs, std_ey_secs = compute_e_fields_secs(rmfb, tf_path, ip1, samp, hi, low, error_bf,
                                                                           ef_tf, e_nvpwa, stat)

        sax = np.array([ex_secs]).T
        say = np.array([ey_secs]).T
        sa1 = np.array([std_ex_secs]).T
        sa2 = np.array([std_ey_secs]).T

        ################################################################################
        # 5.4) Define the vectors with electric fields including statistics
        ################################################################################
        secs_e[:, v1, :] = np.hstack((sax, say))
        stdsecs_e[:, v1, :] = np.hstack((sa1, sa2))

        ################################################################################
        # 5.5) Compute mean and std of electric fields
        ################################################################################
        if mode == 1:
            ################################################################################
            # 5.5.1) Following Approach #1
            ################################################################################

            if str(e_site[v1]) in h_obs:
                # 5.5.1.1) If the magnetics were recorded use the local tensor
                val = -99
                for v4, ip4 in enumerate(h_obs):
                    if ip4 == ip1:
                        val = v4

                # Compute the electrics
                (loc_exmt, loc_eymt, std_loc_exmt, std_loc_eymt) = compute_e_fields(mh_obs[:, val, :],
                                                                                    0 * mh_obs[:, val, :], tf_path, ip1,
                                                                                    ip1, samp, hi, low, ef_tf, e_nvpwa,
                                                                                    stat)

                ax = np.array([loc_exmt]).T
                ay = np.array([loc_eymt]).T
                a1 = np.array([std_loc_exmt]).T
                a2 = np.array([std_loc_eymt]).T

                loc_emt[:, v1, :] = np.hstack((ax, ay))
                std_loc_emt[:, v1, :] = np.hstack((a1, a2))

                # Define the e vector using several values within the error bars
                for ip in range(0, stat):
                    e_fields[:, v1, :, ip] = (
                                loc_emt[:, v1, :] + np.random.standard_normal([length, 2]) * std_loc_emt[:, v1, :])


            else:
                ################################################################################
                # 5.5.1.2) If the magnetics were NOT recorded use Approach #1
                ################################################################################

                for ip in range(0, stat):
                    # Define the E vector using several values within the error bars
                    e_fields[:, v1, :, ip] = (
                                secs_e[:, v1, :] + np.random.standard_normal([length, 2]) * stdsecs_e[:, v1, :])

            # Define the E vector
            ce_fields = e_fields

        if mode == 2:
            ################################################################################
            # 5.5.2) Following Approach #2
            ################################################################################

            for v2, ip2 in enumerate(rmf):
                # Sites used as a reference magnetic fields to calculate the e_fields
                # Compute the regional electric field
                (c_reg_ex, c_reg_ey, c_std_reg_ex, c_std_reg_ey) = compute_e_fields(comp_rmf[:, v2, :],
                                                                                    std_comp_rmf[:, v2, :], tf_path, ip1,
                                                                                    ip2, samp, hi, low, ef_tf, e_nvpwa,
                                                                                    stat)

                cax = np.array([c_reg_ex]).T
                cay = np.array([c_reg_ey]).T
                ca1 = np.array([c_std_reg_ex]).T
                ca2 = np.array([c_std_reg_ey]).T

                c_reg_e[:, v1, v2, :] = np.hstack((cax, cay))
                c_std_reg_e[:, v1, v2, :] = np.hstack((ca1, ca2))

            for v3, ip3 in enumerate(rmf):

                # Define the E vector using several values within the error bars
                # (for Approach #2 the results form SECS are the local e_field)
                for ip in range(0, stat):
                    e_fields[:, v1, :, v3, ip] = (
                                c_reg_e[:, v1, v3, :] + np.random.standard_normal([length, 2]) * c_std_reg_e[:, v1, v3, :])

                    secs_e_fields[:, v1, :, ip] = (
                    (secs_e[:, v1, :] + np.random.standard_normal([length, 2]) * stdsecs_e[:, v1, :]))

                    # Computed to create regional and local electric fields
                    loc_e_fields[:, v1, :, ip] = (
                                secs_e[:, v1, :] + np.random.standard_normal([length, 2]) * stdsecs_e[:, v1, :])

                    reg_e_fields[:, v1, :, v3, ip] = (
                                c_reg_e[:, v1, v3, :] + np.random.standard_normal([length, 2]) * c_std_reg_e[:, v1, v3, :])

                if str(e_site[v1]) in h_obs:
                    # If the magnetics were recorded use the magnetics and local tensor
                    val = -99
                    for v4, ip4 in enumerate(h_obs):
                        if ip4 == ip1:
                            val = v4

                    (loc_exmt, loc_eymt,
                     std_loc_exmt, std_loc_eymt) = compute_e_fields(pmh_obs[:, val, :], 0 * pmh_obs[:, val, :], tf_path,
                                                                    ip1, ip1, samp, hi, low, ef_tf, e_nvpwa, stat)

                    ax = np.array([loc_exmt]).T
                    ay = np.array([loc_eymt]).T
                    a1 = np.array([std_loc_exmt]).T
                    a2 = np.array([std_loc_eymt]).T

                    loc_emt[:, v1, :] = np.hstack((ax, ay))
                    std_loc_emt[:, v1, :] = np.hstack((a1, a2))

                    for ip in range(0, stat):
                        e_fields[:, v1, :, v3, ip] = (
                                    loc_emt[:, v1, :] + np.random.standard_normal([length, 2]) * std_loc_emt[:, v1, :])
                        secs_e_fields[:, v1, :, ip] = secs_e[:, v1, :] * 0

            ce_fields = np.average(e_fields, axis=3) + secs_e_fields
            reg_ee = np.average(reg_e_fields, axis=3)
            loc_ee = np.copy(loc_e_fields)

    ################################################################################
    # 6) Calculate the average/representative E_field, and the standard deviation
    ################################################################################
    print('Calculating mean and std of the computed electric fields')
    for v4, ip4 in enumerate(e_site):
        for ip2 in range(0, 2):
            for ip3 in range(0, length):
                std_error[ip3, v4, ip2] = np.std(ce_fields[ip3, v4, ip2, :])
                av_e_fields[ip3, v4, ip2] = np.average(ce_fields[ip3, v4, ip2, :])

                if mode == 2:
                    std_reg_e[ip3, v4, ip2] = np.std(reg_ee[ip3, v4, ip2, :])
                    av_reg_e[ip3, v4, ip2] = np.average(reg_ee[ip3, v4, ip2, :])

                    std_loc_e[ip3, v4, ip2] = np.std(loc_ee[ip3, v4, ip2, :])
                    av_loc_e[ip3, v4, ip2] = np.average(loc_ee[ip3, v4, ip2, :])

    ################################################################################
    # 7 Write E fields in the out folder
    ################################################################################
    print('Saving the electric fields')
    for vi, ppi in enumerate(e_site):
        # Total electric field
        ex_path_id = (out_path + str(ppi) + '_Ex_' + str(storm) + '_Approach_#' + str(mode) + '.dat')

        ey_path_id = (out_path + str(ppi) + '_Ey_' + str(storm) + '_Approach_#' + str(mode) + '.dat')

        ex = np.array([av_e_fields[:, vi, 0], std_error[:, vi, 0]])
        ey = np.array([av_e_fields[:, vi, 1], std_error[:, vi, 1]])

        np.savetxt(ex_path_id, ex.T, fmt=['%15.5f', '%15.5f'])
        np.savetxt(ey_path_id, ey.T, fmt=['%15.5f', '%15.5f'])

        ################################################################################
        # Local and Regional electric fields (only for Approach #2)
        ################################################################################
        if mode == 2:  # Approach #2
            loc_ex_path_id = (out_path + str(ppi) + '_Ex_' + str(storm) + '_Approach_#' + str(mode) + 'local.dat')
            loc_ey_path_id = (out_path + str(ppi) + '_Ey_' + str(storm) + '_Approach_#' + str(mode) + 'local.dat')
            reg_ex_path_id = (out_path + str(ppi) + '_Ex_' + str(storm) + '_Approach_#' + str(mode) + 'regional.dat')
            reg_ey_path_id = (out_path + str(ppi) + '_Ey_' + str(storm) + '_Approach_#' + str(mode) + 'regional.dat')

            reg_ex = np.array([av_reg_e[:, vi, 0], std_reg_e[:, vi, 0]])
            reg_ey = np.array([av_reg_e[:, vi, 1], std_reg_e[:, vi, 1]])
            loc_ex = np.array([av_loc_e[:, vi, 0], std_loc_e[:, vi, 0]])
            loc_ey = np.array([av_loc_e[:, vi, 1], std_loc_e[:, vi, 1]])

            np.savetxt(loc_ex_path_id, loc_ex.T, fmt=['%15.5f', '%15.5f'])
            np.savetxt(loc_ey_path_id, loc_ey.T, fmt=['%15.5f', '%15.5f'])
            np.savetxt(reg_ex_path_id, reg_ex.T, fmt=['%15.5f', '%15.5f'])
            np.savetxt(reg_ey_path_id, reg_ey.T, fmt=['%15.5f', '%15.5f'])


if __name__ == "__main__":
    main()
