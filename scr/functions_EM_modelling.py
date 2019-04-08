#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:16:58 2017

Altgorithm related to Campanya et al. 2018 publication in Space Weather AGU Journal
"""
import numpy as np
import scipy.signal
import seaborn as sns
import pandas as pd
from numpy.fft.fftpack import ifft
import geopy.distance
import matplotlib.pyplot as plt
sns.set()

###########################################################################
def scr_fft(x, y, s):
    """ Calculate Fourier transform for two time series (x, y components)
		
		Parameters
		-----------
		x = time series 1
		y = time series 2
		s = sampling rate (in seconds)

		Returns
		-----------
		perd = periods in seconds
		x_fft = fft of time series 1
		y_fft = fft of time series 2

		-----------------------------------------------------------------
	"""
    
    w = scipy.signal.tukey(x.shape[0], 0.1) 
    x = x * w
    y = y * w
    freq = (np.fft.fftfreq(x.shape[0], d=float(s)))
    for i in range(0,len(freq)):
        if freq[i] == 0:
            freq[i] = 1e-99

    perd = freq ** (-1)
    x_fft = np.fft.fftpack.fft(x)
    y_fft = np.fft.fftpack.fft(y)

    return (perd, x_fft, y_fft)

##########################################################################
def read_co(path):
    
    """ Read name of the sites and coordinates of the site from the input
        files. Latitude and longitude should be in degrees.
		
		Parameters
		-----------
		path = path of the site with the name of the sites and coordinates

		Returns
		-----------
		name = Name of the site
		lat = latitude of the site (in degrees)
		lon = longitude of the site (in degrees)

		-----------------------------------------------------------------
    """
   
    a = pd.read_csv(path, 
                    header = None, 
                    skiprows = None, 
                    sep='\s+'
                    )
    
    a = np.array(a)
    name = a[:,0]
    lat =  a[:,1]
    lon =  a[:,2]
    
    return(name, lat, lon)
    
##########################################################################
def read_variable(n, index, data):
    k = 1
    val = []
    for j in range(0, n):
        if k < n:
            val_1 = np.loadtxt([data[index + j + 1]])
            val = np.append(val, val_1)
            k = len(val)
            val = np.array(val)
    return(val)
    
##########################################################################
def read_rmfb(path, name, storm):
    
    """ Read interpolated magnetic fields from SECS
		
		Parameters
		-----------
		path = path for the SECS output files
         storm = Name of the selected storm  
         name = Name of the site  


		Returns
		-----------
		rmfb = Interpolated magnetic fields from SECS

		-----------------------------------------------------------------
	"""

    f = open(path 
             + str(storm) 
             + "magBx" 
             + name , 
             'r'
             )
    
    rmf_bx = np.loadtxt(f, skiprows = 0)
    f.close()
        
    f = open(path 
             + str(storm) 
             + "magBy" 
             + name ,
             'r'
             )
    
    rmf_by = np.loadtxt(f, skiprows = 0)
    f.close()
    
    rmfb = np.array([rmf_bx, rmf_by]).T
    
    return(rmfb)

##########################################################################
def read_magnetics( in_path, sites, mag_path, secs_path, samp, hi, low, 
                   length, storm, var):
    
    """ Read magnetic time series from input folder
		
		Parameters
		-----------
		in_path = Folder with input parameters
		sites = Name of the site
         mag_path = Folder with magnetic fields time series
         secs_path = Folder with inputs - outputs for SECS interpolation
         samp = Sampling rate (seconds)
         hi = maximum period to analyse (seconds)
         low = minimum period to analyse (seconds)
         length = length of the time series
         storm = Folder with input magnetic time series
         var = if (1) write the magnetic time series

		Returns
		-----------
		mag_s = Magnetic time series
		s = Name of the site

		-----------------------------------------------------------------
	"""

    # Check the magnetic observatories used in the experiment
    s, s_lat, s_lon = read_co(in_path + str(sites))
    mag_s = np.zeros([length,len(s),2])

    # Read the time series for these magnetic fields
    for ip in range (0,len(s)):
        mag = str(s[ip])
         
        x_input = str(mag_path) + str(mag) + "Bx.txt"
        y_input = str(mag_path) + str(mag) + "By.txt"

        f = open(x_input,'r')
        dif_bx = np.loadtxt(f)
        f.close()

        f = open(y_input,'r')
        dif_by = np.loadtxt(f)
        f.close()                       

        # Remove detrend
        dif_bx = scipy.signal.detrend(dif_bx)
        dif_by = scipy.signal.detrend(dif_by)
        
        # compute FFT
        perd, dif_bx_fft, dif_by_fft = scr_fft(dif_bx,dif_by,samp)
        
        # Select periods of interest
        factor = np.ones(perd.shape[0])        
        for i, v in enumerate (perd):
            if (v < low) or (v > hi):
                factor[i] = 0

        bx_cfft = (factor * dif_bx_fft) 
        by_cfft = (factor * dif_by_fft) 
        
        # Compute ifft
        bx_c = np.real(np.array(
                ifft(bx_cfft + np.conj(np.roll(bx_cfft[::-1], 1)))))
        
        by_c = np.real(np.array(
                ifft(by_cfft + np.conj(np.roll(by_cfft[::-1], 1)))))

        # Deffine the vector of magnetic field data
        ax = np.array([bx_c]).T
        ay = np.array([by_c]).T
        mag_s[:,ip,:] = np.hstack((ax, ay))

    # Write the magnetic data to be used by the SECS interpolation algorithm
    if var == 1:
        for ip in range (0, len(s)):
            mag = str(s[ip])
            w_path = (str(secs_path) 
                     + "SECS_" 
                     + str(mag) 
                     + str(storm) 
                     + ".dat")
            
            f_id = open(w_path,'w+')
            np.savetxt(f_id, mag_s[:,ip,:], fmt=['%15.5f', '%15.5f'])
            f_id.close()
             
    return(mag_s, s)

###############################################################################
def compute_regional_b_field(in_path, sites, reg_ref, mag_path, tf_path, 
                              samp, hi, low, ef_h, stat):

    """ Calculate regional b field at the magnetic observatories using
        reference magnetic field and inter-station tensor relationships
		
		Parameters
		-----------
		in_path = Folder with input parameters
		sites = Name of the site
         reg_ref = Reference magnetic field (e.g. CLF in Campanya et al. 2018)
         mag_path = Folder with magnetic fields time series
         tf_path = Folder with electromagnetic tensor relationships
         samp = Sampling rate (seconds)
         hi = Maximum period to analyse (seconds)
         low = Minimum period to analyse (seconds)
         ef_h = Error floor for H tensor relationship
         stat = Statistics for error propagation
         
         
		Returns
		-----------
         c =  Computed magnetic fields
         std_c = Standard devaition for the computed magnetic fields
		-----------------------------------------------------------------
	"""

    # Check the magnetic observatories used in the experiment
    s, s_lat, s_lon = read_co(str(in_path) + str(sites))

    # Read the magnetic time series from Reference magnetic field 
    x_input = str(mag_path)  + str(reg_ref) + "Bx.txt"
    y_input = str(mag_path)  + str(reg_ref) + "By.txt"

    f = open(x_input,'r')
    reg_refbx = np.loadtxt(f)
    f.close()

    f = open(y_input,'r')
    reg_refby = np.loadtxt(f)
    f.close()                       
    
    # Deffine variables
    c = np.zeros([len(reg_refbx),len(s),2])
    std_c = np.zeros([len(reg_refbx),len(s),2])

    # Calculate detrend
    reg_refbx = scipy.signal.detrend(reg_refbx)
    reg_refby = scipy.signal.detrend(reg_refby)
           
    # Compute fft
    perd, reg_refbx_fft, reg_refby_fft = scr_fft(reg_refbx,
                                                 reg_refby,
                                                 samp
                                                 )

    # Calculate spectra at the sites of interest                   
    for ip in range (0,len(s)):
        mag = str(s[ip])
            
        # Read inter-station transfer functions
        filename = tf_path + "B" + str(mag) + "B" + str(reg_ref) + "_s.j"
        
        f = open(filename, 'r')
        data = f.readlines()
        f.close()
        
        # Read the components of the tensors
        for index, line in enumerate(data):
            if line.startswith('ZXX'):
                index_zxx = index
            if line.startswith("ZXY"):
                index_zxy = index
            if line.startswith("ZYX"):
                index_zyx = index
            if line.startswith("ZYY"):
                index_zyy = index
            if line.startswith("TZX"):
                index_tzx = index
            if line.startswith("TZY"):
                index_tzy = index
  
        data_zxx = data[index_zxx + 2:index_zxy]
        zxx = np.loadtxt(data_zxx)
        data_zxy = data[index_zxy + 2:index_zyx]
        zxy = np.loadtxt(data_zxy) 
        data_zyx = data[index_zyx + 2:index_zyy]
        zyx = np.loadtxt(data_zyx) 
        data_zyy = data[index_zyy + 2:index_tzx]
        zyy = np.loadtxt(data_zyy)    
        per_z = zxx[:,0]

        # Deffine variables
        factor = np.ones(perd.shape[0])
        zxx_int=np.zeros([perd.shape[0],3])
        zxy_int=np.zeros([perd.shape[0],3])
        zyx_int=np.zeros([perd.shape[0],3])
        zyy_int=np.zeros([perd.shape[0],3])
        
        # select periods of interest
        for i, v in enumerate (perd):
            if (v < low) or (v > hi):
                factor[i] = 0

        for i in range (0,3):
            zxx_int[:,i] = np.interp(perd, per_z, zxx[:,i+1]) * factor
            zxy_int[:,i] = np.interp(perd, per_z, zxy[:,i+1]) * factor
            zyx_int[:,i] = np.interp(perd, per_z, zyx[:,i+1]) * factor
            zyy_int[:,i] = np.interp(perd, per_z, zyy[:,i+1]) * factor

        # Deffine variables
        bx_calc = np.zeros([reg_refbx.shape[0], stat])
        by_calc = np.zeros([reg_refby.shape[0], stat])
        
        for ik in range(0,len(zxx_int)):
            if zxx_int[ik,2] <= ef_h:
                zxx_int[ik,2] = ef_h    
            if zxy_int[ik,2] <= ef_h:
                zxy_int[ik,2] = ef_h 
            if zyx_int[ik,2] <= ef_h:
                zyx_int[ik,2] = ef_h 
            if zyy_int[ik,2] <= ef_h:
                zyy_int[ik,2] = ef_h 
        
        # Calculate B fields using ITF
        for i in range(0,stat):
            
            bx_1 = ((((zxx_int[:,0] 
                   + np.random.standard_normal() 
                   * zxx_int[:,2])
                   + (zxx_int[:,1] 
                   + np.random.standard_normal() 
                   * zxx_int[:,2])*1j) 
                   * reg_refbx_fft) 
                   + (((zxy_int[:,0] 
                   + np.random.standard_normal()
                   * zxy_int[:,2])
                   + (zxy_int[:,1] 
                   + np.random.standard_normal() 
                   * zxy_int[:,2])*1j) 
                   * reg_refby_fft))
                   
            by_1 = ((((zyx_int[:,0] 
                   + np.random.standard_normal() 
                   * zyx_int[:,2]) 
                   + (zyx_int[:,1] 
                   + np.random.standard_normal() 
                   * zyx_int[:,2])*1j) 
                   * reg_refbx_fft) 
                   + (((zyy_int[:,0] 
                   + np.random.standard_normal() 
                   * zyy_int[:,2]) 
                   + (zyy_int[:,1] 
                   + np.random.standard_normal() 
                   * zyy_int[:,2])*1j) 
                   * reg_refby_fft))
            
            bx_calc[:,i] = ifft(bx_1 + np.conj(np.roll(bx_1[::-1],1))).real 
            by_calc[:,i] = ifft(by_1 + np.conj(np.roll(by_1[::-1],1))).real 
                
        # Define mean value
        me_bx=np.zeros(reg_refbx.shape[0])
        me_by=np.zeros(reg_refbx.shape[0])
        
        for i in range(0,reg_refbx.shape[0]):
            me_bx[i]=bx_calc[i,:].mean()
            me_by[i]=by_calc[i,:].mean()
            
        # Calculate standard deviation
        bx_s=np.zeros(reg_refbx.shape[0])
        by_s=np.zeros(reg_refbx.shape[0])
        
        for i in range(0,reg_refbx.shape[0]):
            bx_s[i]=bx_calc[i,:].std()
            by_s[i]=by_calc[i,:].std()  
        
        # Deffine the output vectors
        tf_bx = np.array([np.copy(me_bx)]).T
        tf_by = np.array([np.copy(me_by)]).T
        std_bx = np.array([np.copy(bx_s)]).T
        std_by = np.array([np.copy(by_s)]).T
        
        # Arrays with 1) computed b fields and 2) standard deviation
        c[:,ip,:] = np.hstack((tf_bx, tf_by))
        std_c[:,ip,:] = np.hstack((std_bx, std_by))

    return (c, std_c)

##########################################################################
def difference_b_fields(ms, s, c, std_c, secs_path, stat, storm):

    """ Calculate diference between b fields. i.e. between measured at the
        magnetic observatories and the computed regional signal
		
		Parameters
		-----------
		ms = Magnetic time series
		s = Name of the site
         c = Computed magnetic fields
         std_c = Standard devaition for the computed magnetic fields
         secs_path = Folder with inputs - outputs for SECS interpolation
         stat = Statistics for error propagation
         storm = Folder with input magnetic time series        
         
		Returns
		-----------
         ms =  Magnetic time series
         mean_ds = Difference between magnetic fields (mean value)
         std_ds = Standard deviation of difference between magnetic fields
		-----------------------------------------------------------------
	"""

    
    ds = np.zeros([ms.shape[0], ms.shape[1], ms.shape[2],stat])
    
    # Calculate the difference between modelled magnetic fields using 
    # inter-station transfer functions and measured magnetic fields
    
    for i in range(0,stat):
        ds[:,:,:,i] = ms - (c + np.random.standard_normal(c.shape)*std_c)
    
    mean_ds = ds.mean(3)
    std_ds = ds.std(3)

    # Write the differences between magnetic fields 
    # To be used for the SECS interpolation method on Approach #2
    for ip in range (0, ms.shape[1]):
        mag = str(s[ip])
        w_path = str(secs_path) + "dSECS_" + str(mag) + str(storm) + ".dat"
        
        f_id = open(w_path,'w+')
        np.savetxt(f_id, mean_ds[:,ip,:], fmt=['%15.5f' , '%15.5f'])
        f_id.close()

    return(ms, mean_ds, std_ds)

##########################################################################
def error_secs_interpolation(e_lat, e_lon, in_path, mode, obs):

    """ Compute the errors caused by SECS interpolation (based on Figure 7, 
        Campanya et al 2018)
		
		Parameters
		-----------
         e_lat = latitude sites to compute electric fields
         e_lon = longitude sites to compute magnetic fields
         in_path = Folder with input parameters
         mode = (1) for Approach #1 and (2) for Approach #2
         obs_f = Name of the file with name and coordinated of the 
                 magnetic observatories

		Returns
		-----------
         error_bf = error associated with the SECS interpolation approach
		-----------------------------------------------------------------
	"""
    
    # Compute the distance of each site to the magnetic observatories

    Obs, Obs_lat, Obs_lon = read_co(in_path + str(obs))
    dist = np.zeros(len(Obs))

    for i in range(0,len(Obs)):
        coo_a = (e_lat, e_lon)
        coo_b = (Obs_lat[i], Obs_lon[i])        
        dist[i] = geopy.distance.vincenty(coo_a, coo_b).km
    
    if mode == 1: # Approach 1   
        snr = -1.75e-2 * dist.min() + 12.81                         

    if mode == 2: # Approach 2
        snr = -1.48e-2 * dist.min() + 9.70
    # Compute the error                             
    error_bf = 1.0/np.sqrt([10**(snr/10)])
    
    return(error_bf)

##########################################################################
def compute_e_fields_secs(sb, tf_path, e_site, samp, hi, low, error_bf,
                          ef_tf, e_nvpwa, stat):

    """ Compute E fields using magnetic time series from SECS
		
		Parameters
		-----------
         sb = magnetic time series from SECS
         tf_path = Folder with electromagnetic tensor relationships
         e_site = Name of the site to compute electric fields
         samp = Sampling rate (seconds)
         hi = Maximum period to analyse (seconds)
         low = Minimum period to analyse (seconds)
         error_bf = Error associated with the SECS interpolation approach
         ef_tf = Error floor for the MT and quasi-MT tensor relationships
         e_ncpwa =  Error floor for the Non-plane wave approximation
         stat = Statistics for error propagation

         
		Returns
		-----------
         tf_ex = electric time series x component (mean value)
         tf_ey = electric time series y component (mean value)
         std_ex = standard deviation of the electric time series x component
         std_ey = standard deviation of the electric time series y component
         
         -----------------------------------------------------------
	"""
    
    s_bx = sb[:,0]*1e-9  # convert to nT
    s_by = sb[:,1]*1e-9  # convert to nT

    # Remove detrend
    s_bx = scipy.signal.detrend(s_bx)
    s_by = scipy.signal.detrend(s_by)
    
    # Tukey window to avoid inestabilities at the edges of the time series
    window = scipy.signal.tukey(s_by.shape[0], 0.1) 
    s_bx = window * s_bx 
    s_by = window * s_by
    
    # get Frequencies / periods
    freqB = np.fft.fftfreq(s_bx.shape[0], d = samp)
    for i in range(0,len(freqB)):
        if freqB[i] == 0:
            freqB[i] = 1e-99

    perB = freqB ** -1
    
    # Compute fft
    s_bx_fft = np.fft.fftpack.fft(s_bx)
    s_by_fft = np.fft.fftpack.fft(s_by)

    
    # Read tensors
    file_format = -1
    try:
        filename = (str(tf_path) 
                    + "E"
                    + str(e_site) 
                    + "B"
                    + str(e_site) 
                    + "_s.j")
        
        print(filename)
        f = open(filename, 'r')
        data = f.readlines()
        f.close()
        print('j file')
        file_format = 0
    except:
        check = 0

    try:
         filename = (str(tf_path) 
                     + "E"
                     + str(e_site) 
                     + "B"
                     + str(e_site)
                     + "_s.edi")
         
         f = open(filename, 'r')
         data = f.readlines()
         f.close()
         print('edi file')
         file_format = 1
    except:
        check = 1

    if file_format == 0:     
        try:
            for index, line in enumerate(data):
                if line.startswith("ZXX"):
                    index_zxx = index
                if line.startswith("ZXY"):
                    index_zxy = index
                if line.startswith("ZYX"):
                    index_zyx = index
                if line.startswith("ZYY"):
                    index_zyy = index
                if line.startswith("TZX"):
                    index_tzx = index
                if line.startswith("TZY"):
                    index_tzy = index
                  
            data_zxx = data[index_zxx + 2 : index_zxy]
            zxx = np.loadtxt(data_zxx)    
            data_zxy = data[index_zxy + 2 : index_zyx]
            zxy = np.loadtxt(data_zxy)     
            data_zyx = data[index_zyx + 2 : index_zyy]
            zyx = np.loadtxt(data_zyx)    
            data_zyy = data[index_zyy + 2 : index_tzx]
            zyy = np.loadtxt(data_zyy)  
            per_z = zxx[:,0]
            
            zxx[:,1:3] = (1) * zxx[:,1:3]
            zxy[:,1:3] = (1) * zxy[:,1:3]
            zyx[:,1:3] = (1) * zyx[:,1:3]
            zyy[:,1:3] = (1) * zyy[:,1:3]

            print(per_z)
        except:
            check = 2     
            
    if file_format == 1:     
        try:
            for index, line in enumerate(data):
                if 'NFREQ' in line[:5]:
                    for j in range(0, len(line)):
                        if line[j] == '=':
                            n_freq = int(line[j+1::])
        
                if 'NPER' in line[:5]:
                    for j in range(0, len(line)):
                        if line[j] == '=':
                            n_freq = int(line[j+1::])
        
                if '>FREQ' in line[:5]:
                    freq = read_variable(n_freq, index, data)
                    
                if '>PERI' in line[:5]:
                    per = read_variable(n_freq, index, data)
                    freq = 1./per
                    
                if '>ZROT' in line[:5]:
                    zrot = read_variable(n_freq, index, data)
        
                if '>ZXXR' in line[:5]:
                    zxxr = read_variable(n_freq, index, data)
                    
                if '>ZXXI' in line[:5]:
                    zxxi = read_variable(n_freq, index, data)
        
                if '>ZXX.V' in line[:6]:
                    zxxv = read_variable(n_freq, index, data)
                    zxxstd = np.sqrt(zxxv)
    
                if '>ZXYR' in line[:5]:
                    zxyr = read_variable(n_freq, index, data)
        
                if '>ZXYI' in line[:5]:
                    zxyi = read_variable(n_freq, index, data)
        
                if '>ZXY.V' in line[:6]:
                    zxyv = read_variable(n_freq, index, data)
                    zxystd = np.sqrt(zxyv)
    
                if '>ZYXR' in line[:5]:
                    zyxr = read_variable(n_freq, index, data)
                if '>ZYXI' in line[:5]:
                    zyxi = read_variable(n_freq, index, data)
                    
                if '>ZYX.V' in line[:6]:
                    zyxv = read_variable(n_freq, index, data)
                    zyxstd = np.sqrt(zyxv)
    
                if '>ZYYR' in line[:5]:
                    zyyr = read_variable(n_freq, index, data)
                    
                if '>ZYYI' in line[:5]:
                    zyyi = read_variable(n_freq, index, data)
                    
                if '>ZYY.V' in line[:6]:
                    zyyv = read_variable(n_freq, index, data)
                    zyystd = np.sqrt(zyyv)
            try:
                periods = 1./freq
            except:
                periods = per    
            
            zxx = np.column_stack([periods, -1*zxxr, -1*zxxi, zxxstd])
            zxy = np.column_stack([periods, -1*zxyr, -1*zxyi, zxystd])
            zyx = np.column_stack([periods, -1*zyxr, -1*zyxi, zyxstd])
            zyy = np.column_stack([periods, -1*zyyr, -1*zyyi, zyystd])
            per_z = zxx[:,0]
    
            if per_z[0] > per_z[1]:
                per_z = per_z[::-1]
                zxx = zxx[::-1]
                zxy = zxy[::-1]
                zyx = zyx[::-1]
                zyy = zyy[::-1]

        except:
            check = 3
            
    if file_format == -1:
        print('Cannot read the MT impedance tensor for site:' + str(e_site))
        print('MT impdeance tensor must be j. or edi. file')

    # Select the periods of interest
    factor = np.ones(perB.shape[0])
    zxx_int=np.zeros([perB.shape[0],3])
    zxy_int=np.zeros([perB.shape[0],3])
    zyx_int=np.zeros([perB.shape[0],3])
    zyy_int=np.zeros([perB.shape[0],3])
    
    for i, v in enumerate (perB):
        if (v < low) or (v > hi):
           factor[i] = 0
    
    for i in range (0,3):
        zxx_int[:,i] = np.interp(perB, per_z, zxx[:,i+1])*factor
        zxy_int[:,i] = np.interp(perB, per_z, zxy[:,i+1])*factor
        zyx_int[:,i] = np.interp(perB, per_z, zyx[:,i+1])*factor
        zyy_int[:,i] = np.interp(perB, per_z, zyy[:,i+1])*factor

    # Deffine Variables
    ex_calc=np.zeros([s_bx.shape[0],stat])
    ey_calc=np.zeros([s_bx.shape[0],stat])
   
    # Deffine Error floor
    zzz_det = np.sqrt(np.abs(((zxy_int[:,0] + zxy_int[:,1]*1j)
              * (zyx_int[:,0] + zyx_int[:,1]*1j)) 
              - ((zxx_int[:,0] + zxx_int[:,1]*1j) 
              * (zyy_int[:,0] + zyy_int[:,1]*1j))))

    for ik in range(0,len(zxx_int)):
        if zxx_int[ik,2] <= zzz_det[ik]*ef_tf:
            zxx_int[ik,2] = zzz_det[ik]*ef_tf 
        if zxy_int[ik,2] <= zzz_det[ik]*ef_tf:
            zxy_int[ik,2] = zzz_det[ik]*ef_tf 
        if zyx_int[ik,2] <= zzz_det[ik]*ef_tf:
            zyx_int[ik,2] = zzz_det[ik]*ef_tf 
        if zyy_int[ik,2] <= zzz_det[ik]*ef_tf:
            zyy_int[ik,2] = zzz_det[ik]*ef_tf 
                       
    # Compute electric fields    
    for i in range(0,stat):
        ex_1 = ((((zxx_int[:,0] + np.random.standard_normal()*zxx_int[:,2])
              + (zxx_int[:,1] + np.random.standard_normal()*zxx_int[:,2])*1j) 
              * (s_bx_fft + s_bx_fft*np.random.standard_normal()*(error_bf))) 
              + (((zxy_int[:,0]+np.random.standard_normal()*zxy_int[:,2])
              + (zxy_int[:,1]+np.random.standard_normal()*zxy_int[:,2])*1j) 
              * (s_by_fft + s_by_fft*np.random.standard_normal()*(error_bf))))
    
        ey_1 = ((((zyx_int[:,0]+np.random.standard_normal()*zyx_int[:,2])
              + (zyx_int[:,1]+np.random.standard_normal()*zyx_int[:,2])*1j) 
              * (s_bx_fft + s_bx_fft*np.random.standard_normal()*(error_bf))) 
              + (((zyy_int[:,0]+np.random.standard_normal()*zyy_int[:,2])
              + (zyy_int[:,1]+np.random.standard_normal()*zyy_int[:,2])*1j) 
              * (s_by_fft + s_by_fft*np.random.standard_normal()*(error_bf))))
        
        # Add error associated with the non-validity of the planewave approx.
        ex_1 = ex_1 + np.random.standard_normal() * ex_1 * e_nvpwa
        ey_1 = ey_1 + np.random.standard_normal() * ey_1 * e_nvpwa

        # Compute ifft
        ex_calc[:,i] = (np.real(ifft(ex_1 + np.conj(np.roll(ex_1[::-1],1)))) 
                       * (1000. / (4*np.pi*1e-7)))
        ey_calc[:,i] = (np.real(ifft(ey_1 + np.conj(np.roll(ey_1[::-1],1)))) 
                       * (1000. / (4*np.pi*1e-7)))

    # Define mean value
    mex=np.zeros(s_bx.shape[0])
    mey=np.zeros(s_bx.shape[0])
    
    for i in range(s_bx.shape[0]):
        mex[i]=ex_calc[i,:].mean()
        mey[i]=ey_calc[i,:].mean()
        
    # Calculate standard deviation (errorbars)
    ex_s=np.zeros(s_bx.shape[0])
    ey_s=np.zeros(s_bx.shape[0])
    
    for i in range(s_bx.shape[0]):
        ex_s[i]=ex_calc[i,:].std()
        ey_s[i]=ey_calc[i,:].std()
    
    # Deffine the outputs
    tf_ex = np.array(np.copy(mex))
    tf_ey = np.array(np.copy(mey))
    std_ex = np.array(np.copy(ex_s))
    std_ey = np.array(np.copy(ey_s))

    return (tf_ex, tf_ey, std_ex, std_ey)


#######################################################################
def compute_e_fields(sb,std_sb, tf_path, e_site, b_site, samp, hi, low,
                     ef_tf, e_nvpwa, stat):

    """ Compute E fields using magnetic time series from SECS
		
		Parameters
		-----------
         sb = magnetic time series
         std_sb = standard deviation of the magnetic time series
         tf_path = Folder with electromagnetic tensor relationships
         e_site = Name of the site to compute electric fields
         b_site = Name of the site for which the magnetic fields where
                  used to compute electric fields
         samp = Sampling rate (seconds)
         hi = Maximum period to analyse (seconds)
         low = Minimum period to analyse (seconds)
         ef_tf = Error floor for the MT and quasi-MT tensor relationships
         e_ncpwa =  Error floor for the Non-plane wave approximation
         stat = Statistics for error propagation

         
		Returns
		-----------
         tf_ex = electric time series x component (mean value)
         tf_ey = electric time series y component (mean value)
         std_ex = standard deviation of the electric time series x component
         std_ey = standard deviation of the electric time series y component
         
         -----------------------------------------------------------
	"""
    ###################################################################
    # Read magnetic data and convert from nT to T
    s_bx = sb[:,0]*1e-9  
    s_by = sb[:,1]*1e-9  
    std_s_bx = std_sb[:,0]*1e-9 
    std_s_by = std_sb[:,1]*1e-9 
                    
    ###################################################################
    # Compute electric fields
    s_bx = scipy.signal.detrend(s_bx)   
    s_by = scipy.signal.detrend(s_by)  
    
    # Tukey window
    window = scipy.signal.tukey(s_by.shape[0], 0.1) 
    s_bx = window * s_bx 
    s_by = window * s_by
    
    # get Frequencies
    freqB = np.fft.fftfreq(s_bx.shape[0], d = samp)
    for i in range(0,len(freqB)):
        if freqB[i] == 0:
            freqB[i] = 1e-99

    perB = freqB ** -1
    
    # Read tensors
    filename = str(tf_path) + "E"+ str(e_site) + "B"+ str(b_site) + "_s.j"
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    
    for index, line in enumerate(data):
        if line.startswith("ZXX"):
            index_zxx = index
        if line.startswith("ZXY"):
            index_zxy = index
        if line.startswith("ZYX"):
            index_zyx = index
        if line.startswith("ZYY"):
            index_zyy = index
        if line.startswith("TZX"):
            index_tzx = index
        if line.startswith("TZY"):
            index_tzy = index    
      
    data_zxx = data[index_zxx + 2 : index_zxy]
    zxx = np.loadtxt(data_zxx)    
    data_zxy = data[index_zxy + 2 : index_zyx]
    zxy = np.loadtxt(data_zxy)     
    data_zyx = data[index_zyx + 2 : index_zyy]
    zyx = np.loadtxt(data_zyx)    
    data_zyy = data[index_zyy + 2 : index_tzx]
    zyy = np.loadtxt(data_zyy)  
    per_z = zxx[:,0]
    
    # Select the periods we are interested on
    factor = np.ones(perB.shape[0])
    zxx_int=np.zeros([perB.shape[0],3])
    zxy_int=np.zeros([perB.shape[0],3])
    zyx_int=np.zeros([perB.shape[0],3])
    zyy_int=np.zeros([perB.shape[0],3])
    
    for i, v in enumerate (perB):
        if (v < low) or (v > hi):
           factor[i] = 0
    
    for i in range (0,3):
        zxx_int[:,i] = np.interp(perB, per_z, zxx[:,i+1])*factor
        zxy_int[:,i] = np.interp(perB, per_z, zxy[:,i+1])*factor
        zyx_int[:,i] = np.interp(perB, per_z, zyx[:,i+1])*factor
        zyy_int[:,i] = np.interp(perB, per_z, zyy[:,i+1])*factor
        
    # Calculate and propagate errors
    ex_calc=np.zeros([s_bx.shape[0],stat])
    ey_calc=np.zeros([s_bx.shape[0],stat])

    # Calculate error floor
    zzz_det = np.sqrt(np.abs(((zxy_int[:,0] + zxy_int[:,1]*1j) 
             * (zyx_int[:,0]+zyx_int[:,1]*1j)) 
             - ((zxx_int[:,0]+zxx_int[:,1]*1j) 
             * (zyy_int[:,0]+zyy_int[:,1]*1j))))

    for ik in range(0,len(zxx_int)):
        if zxx_int[ik,2] <= zzz_det[ik]*ef_tf:
            zxx_int[ik,2] = zzz_det[ik]*ef_tf 
        if zxy_int[ik,2] <= zzz_det[ik]*ef_tf:
            zxy_int[ik,2] = zzz_det[ik]*ef_tf 
        if zyx_int[ik,2] <= zzz_det[ik]*ef_tf:
            zyx_int[ik,2] = zzz_det[ik]*ef_tf 
        if zyy_int[ik,2] <= zzz_det[ik]*ef_tf:
            zyy_int[ik,2] = zzz_det[ik]*ef_tf 
    
    for i in range(0,stat):
        # Compute the electrics
        s_bx_fft = np.fft.fftpack.fft(s_bx 
                  + np.random.standard_normal(len(s_bx))*std_s_bx)
        s_by_fft = np.fft.fftpack.fft(s_by 
                  + np.random.standard_normal(len(s_by))*std_s_by)

        ex_1 = ((((zxx_int[:,0] + np.random.standard_normal() * zxx_int[:,2]) 
              +(zxx_int[:,1] + np.random.standard_normal() * zxx_int[:,2])*1j) 
              * (s_bx_fft )) 
              + (((zxy_int[:,0] + np.random.standard_normal()*zxy_int[:,2])
              + (zxy_int[:,1] + np.random.standard_normal()*zxy_int[:,2])*1j) 
              * (s_by_fft )))
        
        ey_1 = ((((zyx_int[:,0] + np.random.standard_normal()*zyx_int[:,2])
              + (zyx_int[:,1] + np.random.standard_normal()*zyx_int[:,2])*1j) 
              * (s_bx_fft )) 
              + (((zyy_int[:,0] + np.random.standard_normal()*zyy_int[:,2])
              + (zyy_int[:,1] + np.random.standard_normal()*zyy_int[:,2])*1j) 
              * (s_by_fft )))
        
        # Add error associated with the non-validity of the plane wave approx.
        ex_1 = ex_1 + np.random.standard_normal() * ex_1 * e_nvpwa
        ey_1 = ey_1 + np.random.standard_normal() * ey_1 * e_nvpwa

        # Compute ifft
        ex_calc[:,i] =  (ifft(ex_1 + np.conj(np.roll(ex_1[::-1],1))).real 
                       * (1000. / (4*np.pi*1e-7)))
        ey_calc[:,i] =  (ifft(ey_1 + np.conj(np.roll(ey_1[::-1],1))).real 
                       * (1000. / (4*np.pi*1e-7)))
    
    # Calculate the electric fields (mean value)
    mex=np.zeros(s_bx.shape[0])
    mey=np.zeros(s_bx.shape[0])
    
    for i in range(s_bx.shape[0]):
        mex[i]=ex_calc[i,:].mean()
        mey[i]=ey_calc[i,:].mean()
        
    # Calculate standard deviation / errorbars
    ex_s=np.zeros(s_bx.shape[0])
    ey_s=np.zeros(s_bx.shape[0])
    
    for i in range(s_bx.shape[0]):
        ex_s[i]=ex_calc[i,:].std()
        ey_s[i]=ey_calc[i,:].std()
    
    # Deffine outputs
    tf_ex = np.array(np.copy(mex))
    tf_ey = np.array(np.copy(mey))
    std_ex = np.array(np.copy(ex_s))
    std_ey = np.array(np.copy(ey_s))

    return (tf_ex, tf_ey, std_ex, std_ey)