#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:16:58 2017

Altgorithm related to Campanya et al. 2018 publication in Space Weather AGU Journal
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy.fft.fftpack import ifft as ifft
import scipy.signal
import seaborn as sns
import matplotlib.patches as mpatches
from scipy import signal
from scipy.stats import norm
import pandas as pd
sns.set()

#######################################################################
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

#######################################################################
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
    
    a = pd.read_csv(path, header = None, skiprows = None, sep='\s+')
    a = np.array(a)
    name = a[:,0]
    lat =  a[:,1]
    lon =  a[:,2]
    
    return(name, lat, lon)

#######################################################################
def read_electrics(e_path, e_site, samp, hi, low):

    """ Read measured electric time series
		
		Parameters
		-----------
		e_path = Folder with electric field time series
         e_site = Name of the site to compute electric fields
         samp = sampling rate (in seconds)
         hi = Maximum period to analyse (seconds)
         low = Minimum period to analyse (seconds)


		Returns
		-----------
		m_e = Measured electric time series (including x and y components)

		-----------------------------------------------------------------
    """

    
    # Read the electric field time series
    x_input = str(e_path) + str(e_site) + "Ex.txt"
    y_input = str(e_path) + str(e_site) + "Ey.txt"

    # Read electric time series
    f = open(x_input,'r')
    ex = np.loadtxt(f)
    f.close()

    f = open(y_input,'r')
    ey = np.loadtxt(f)
    f.close()   

    # Select periods of interest (frequency domain)
    perd, ex_fft, ey_fft = scr_fft(ex,ey,samp)

    factor = np.ones(perd.shape[0])        
    for i, v in enumerate (perd):
        if (v < low) or (v > hi):
            factor[i] = 0

    ex_cfft = (factor * ex_fft) 
    ey_cfft = (factor * ey_fft)   
    ex_m = np.array(ifft(ex_cfft + np.conj(np.roll(ex_cfft[::-1], 1))).real) 
    ey_m = np.array(ifft(ey_cfft + np.conj(np.roll(ey_cfft[::-1], 1))).real)
  
    # Create vector with electric time series
    m_e = np.array([ex_m, ey_m])

    return(m_e.T)

#######################################################################
def read_comp_elc(path, e_site, samp, storm, mode):
    """ Read modelled electric time series
		
		Parameters
		-----------
		path = Folder with modelled electric field time series
         e_site = Name of the site of interest
         samp = sampling rate (in seconds)
         storm = Name of the selected storm  
         mode = (1) for Approach #1 and (2) for Approach #2


		Returns
		-----------
		c_e = Modelled electric time series (including x and y components)
         stde = Standard deviation of the modelled electric time series
		-----------------------------------------------------------------
    """

    x_input = (str(path) + str(e_site) + "_Ex_" 
              + str(storm) + "_Approach_#" + str(mode) + ".dat")
            
    y_input = (str(path) + str(e_site) + "_Ey_" 
              + str(storm) + "_Approach_#" + str(mode) + ".dat")

    f = open(x_input,'r')
    ex = np.loadtxt(f)
    f.close()

    f = open(y_input,'r')
    ey = np.loadtxt(f)
    f.close()

    ce = np.array([ex[:,0], ey[:,0]])
    stde = np.array([ex[:,1], ey[:,1]])

    return(ce.T, stde.T)

#######################################################################
def read_comp_mag(path, e_site, samp, storm, mode):
    """ Read modelled electric time series
		
		Parameters
		-----------
		path = Folder with modelled electric field time series
         e_site = Name of the site of interest
         samp = sampling rate (in seconds)
         storm = Name of the selected storm  
         mode = (1) for Approach #1 and (2) for Approach #2


		Returns
		-----------
		c_e = Modelled electric time series (including x and y components)
         stde = Standard deviation of the modelled electric time series
		-----------------------------------------------------------------
    """

    x_input = (str(path)  + 'SECS_' + str(storm) + 'magBx' + str(e_site))
            
    y_input = (str(path) + 'SECS_' + str(storm) + 'magBy' + str(e_site))

    f = open(x_input,'r')
    ex = np.loadtxt(f)
    f.close()

    f = open(y_input,'r')
    ey = np.loadtxt(f)
    f.close()

    ce = np.array([ex, ey])

    return(ce.T)

#######################################################################
def read_reg_elc(path, e_site, samp, storm, mode):

    """ Read modelled regional electric time series
		
		Parameters
		-----------
		path = Folder with modelled regional electric field time series
         e_site = Name of the site of interest
         samp = sampling rate (in seconds)
         storm = Name of the selected storm  
         mode = (1) for Approach #1 and (2) for Approach #2


		Returns
		-----------
		c_e = Modelled electric time series (including x and y components)
         stde = Standard deviation of the modelled electric time series
		-----------------------------------------------------------------
    """
    
    # Read the time series magnetic fields
    x_input = (str(path) + str(e_site) + "_Ex_" 
              + str(storm) + "_Approach_#" + str(mode) + "regional.dat")

    y_input = (str(path) + str(e_site) + "_Ey_" 
              + str(storm) + "_Approach_#" + str(mode) + "regional.dat")

    f = open(x_input,'r')
    ex = np.loadtxt(f)
    f.close()

    f = open(y_input,'r')
    ey = np.loadtxt(f)
    f.close()
      
    ce = np.array([ex[:,0], ey[:,0]])
    stde = np.array([ex[:,1], ey[:,1]])

    return(ce.T, stde.T)

#######################################################################
def read_loc_elc(data_path, e_site, samp, storm, mode):
    """ Read modelled local electric time series
		
		Parameters
		-----------
		path = Folder with modelled local electric field time series
         e_site = Name of the site of interest
         samp = sampling rate (in seconds)
         storm = Name of the selected storm  
         mode = (1) for Approach #1 and (2) for Approach #2


		Returns
		-----------
		c_e = Modelled electric time series (including x and y components)
         stde = Standard deviation of the modelled electric time series
		-----------------------------------------------------------------
    """

    x_input = (str(data_path) + str(e_site) + "_Ex_" 
               + str(storm) + "_Approach_#" + str(mode) + "local.dat")
    y_input = (str(data_path) + str(e_site) + "_Ey_" 
               + str(storm) + "_Approach_#" + str(mode) + "local.dat")

    f = open(x_input,'r')
    ex = np.loadtxt(f)
    f.close()

    f = open(y_input,'r')
    ey = np.loadtxt(f)
    f.close()
       
    ce = np.array([ex[:,0], ey[:,0]])
    stde = np.array([ex[:,1], ey[:,1]])

    return(ce.T, stde.T)


#######################################################################
def polt_time_series(rmf, meas_f, mod_f, std_error, mint, maxt, val):

    """ Plot time series comparing measured and modelled electric fields
		
		Parameters
		-----------
		rmf = Name of the site of interest
         meas_f = Measured electric field at the site of interest
         mod_f = Modelled electric field at the site of interest
         std_error = etandard deviation for the modelled electric field
         mint = Time series properties. Starting point for the analysis
         maxt = Time series properties. End point for the analysis
         val = If 0 for total, 1 for regional, 2 for local 
		Returns
		-----------
		
		-----------------------------------------------------------------
    """

    lenlen = np.arange(0,len(meas_f),1)

    plt.style.use('ggplot')
    for ik in range(0,mod_f.shape[1]):
        n_rmsx = ((np.sqrt(np.average((meas_f[mint:maxt,ik,0]
                 - mod_f[mint:maxt,ik,0])**2)))
                 / (np.average((std_error[mint:maxt,ik,0]))))   
        
        rmsx = ((np.sqrt(np.average((meas_f[mint:maxt,ik,0]
                 - mod_f[mint:maxt,ik,0])**2))))   
        
        n_rmsy = ((np.sqrt(np.average((meas_f[mint:maxt,ik,1] 
                 - mod_f[mint:maxt,ik,1])**2)))
                 /(np.average(std_error[mint:maxt,ik,1])))   
        
        rmsy =  ((np.sqrt(np.average((meas_f[mint:maxt,ik,1]
                 - mod_f[mint:maxt,ik,1])**2))))   
                
        plt.figure()
        if val == 0:
            plt.suptitle('Electric time series at ' + str(rmf[ik]))
        if val == 1:
            plt.suptitle('Regional Electric time series at ' + str(rmf[ik]))
        if val == 2:
            plt.suptitle('Local Electric time series at ' + str(rmf[ik]))
       
        # Time series x direction
        plt.subplot(4,1,1)
        plt.plot(lenlen[mint:maxt],
                 meas_f[mint:maxt,ik,0],
                 'k',
                 lw = 1.0, 
                 label='Measured_x'
                 )        
        plt.plot(lenlen[mint:maxt],
                 mod_f[mint:maxt,ik,0], 
                 'r',
                 lw = 1.0, 
                 label='x_cal_' + 'nRMS: %1.1f ' %(n_rmsx)
                 + ' RMS: %1.1f ' %(rmsx)
                 )        
        plt.fill_between(lenlen[mint:maxt],
                         mod_f[mint:maxt,ik,0]-(std_error[mint:maxt,ik,0]),
                         mod_f[mint:maxt,ik,0]+(std_error[mint:maxt,ik,0]), 
                         facecolor = 'r',
                         alpha=0.5, 
                         linewidth=0.5
                         )
        plt.legend()
        plt.ylabel('Ex [mV/km]')
        plt.xlabel('Time [min]')

        
        # Time series y direction
        plt.subplot(4,1,2)
        plt.plot(lenlen[mint:maxt],
                 meas_f[mint:maxt,ik,1], 
                 'k',
                 lw = 1.0 ,
                 label='Measured_y')
        plt.plot(lenlen[mint:maxt],
                 mod_f[mint:maxt,ik,1], 
                 'c',
                 lw = 1.0, 
                 label='Y_cal_' +' nRMS: %1.1f ' %(n_rmsy)
                 +'RMS: %1.1f ' %(rmsy))
        plt.fill_between(lenlen[mint:maxt],
                         mod_f[mint:maxt,ik,1]-(std_error[mint:maxt,ik,1]),
                         mod_f[mint:maxt,ik,1]+(std_error[mint:maxt,ik,1]), 
                         facecolor = 'c',
                         alpha=0.5, 
                         linewidth=0.5)
        plt.legend()
        plt.ylabel('Ey [mV/km]')
        plt.xlabel('Time [min]')
        
        # Propagated errors (standard deviation)
        plt.subplot(4,1,3)
        plt.plot(lenlen[mint:maxt],
                 std_error[mint:maxt,ik,0], 
                 'r--',
                 lw = 1.0, 
                 label='ex error(std)')
        plt.plot(lenlen[mint:maxt],
                 std_error[mint:maxt,ik,1], 
                 'c--',
                 lw = 1.0, 
                 label='ey error(std)')
        plt.legend()
        er = np.max([std_error[mint:maxt,ik,1],
                      std_error[mint:maxt,ik,0],
                      mod_f[mint:maxt,ik,0] - meas_f[mint:maxt,ik,0],
                      mod_f[mint:maxt,ik,1] - meas_f[mint:maxt,ik,1]] 
                      )
        plt.ylim(0,er)
        plt.ylabel('[mV/km]')
        plt.xlabel('Time [min]')
       
        # Difference between measured and modelled e_fields        
        plt.subplot(4,1,4)
        plt.plot(lenlen[mint:maxt], 
                 np.abs(mod_f[mint:maxt,ik,0] - meas_f[mint:maxt,ik,0]), 
                 'r--',
                 lw = 1.0, 
                 label = '|Diff. ex|')
        plt.plot(lenlen[mint:maxt], 
                 np.abs(mod_f[mint:maxt,ik,1] - meas_f[mint:maxt,ik,1]), 
                 'c--',
                 lw = 1.0, 
                 label = '|Diff. ey|')
        plt.legend()
        plt.ylim(0,er)
        plt.ylabel('[mV/km]')
        plt.xlabel('Time [min]')
        
        if val == 0:
            # Signal-to-noise ratio    
            snr_x = 10*np.log10((np.sum(meas_f[mint:maxt,ik,0]**2))
                    / np.sum((mod_f[mint:maxt,ik,0] - meas_f[mint:maxt,ik,0])**2))
       
            snr_y = 10*np.log10((np.sum(meas_f[mint:maxt,ik,1]**2))
                    / np.sum((mod_f[mint:maxt,ik,1] - meas_f[mint:maxt,ik,1])**2))
         
            print('SNR_'+ str(rmf[ik]))
            print(snr_x)
            print(snr_y)
            
            # Coherence        
            coh_x = np.corrcoef(meas_f[mint:maxt,ik,0], 
                                mod_f[mint:maxt,ik,0])[0, 1]
    
            coh_y = np.corrcoef(meas_f[mint:maxt,ik,1], 
                                mod_f[mint:maxt,ik,1])[0, 1]
            
            print('Coh_'+ str(rmf[ik]))
            print(coh_x)
            print(coh_y)

    return()

#######################################################################
def polt_spectogram(rmf, meas_f, mod_f, samp, mint, maxt):

    """ Plot spectrograms comparing measured and modelled electric fields
		
		Parameters
		-----------
		rmf = Name of the site of interest
         meas_f = Measured electric field at the site of interest
         mod_f = Modelled electric field at the site of interest
         samp = Sampling rate (in seconds)
         mint = Time series properties. Starting point for the analysis
         maxt = Time series properties. End point for the analysis

		Returns
		-----------
		
		-----------------------------------------------------------------
    """
    # Deffine parameters of the time series for spectral analysis
    dp = 256      # data points for the spectral analysis
    fs = 1/samp   # Frequency
    
    # Compute the spectrogram
    for ik in range(0,mod_f.shape[1]):
        # Measured
        xf, xt, xts = signal.spectrogram(meas_f[:,ik,0], 
                                           fs,
                                           nperseg=dp, 
                                           noverlap=(dp/2), 
                                           nfft=None, 
                                           detrend=False, 
                                           return_onesided=False,
                                           scaling='density', 
                                           axis=-1, 
                                           mode = 'psd'
                                           )

        yf, yt, yts = signal.spectrogram(meas_f[:,ik,1], 
                                           fs, 
                                           nperseg=dp, 
                                           noverlap=(dp/2), 
                                           nfft=None, 
                                           detrend=False, 
                                           return_onesided=False, 
                                           scaling='density', 
                                           axis=-1, 
                                           mode='psd'
                                           )
        # Modelled
        xfm, xtm, xtsm = signal.spectrogram(mod_f[:,ik,0], 
                                              fs, 
                                              nperseg=dp, 
                                              noverlap=(dp/2), 
                                              nfft=None, 
                                              detrend=False, 
                                              return_onesided=False, 
                                              scaling='density', 
                                              axis=-1, 
                                              mode='psd'
                                              )
        
        yfm, ytm, ytsm = signal.spectrogram(mod_f[:,ik,1], 
                                              fs, 
                                              nperseg=dp, 
                                              noverlap=(dp/2), 
                                              nfft=None, 
                                              detrend=False, 
                                              return_onesided=False, 
                                              scaling='density', 
                                              axis=-1, 
                                              mode='psd'
                                              )
        
        difx=np.log10(np.abs(np.abs(xtsm)/np.abs(xts)))
        dify=np.log10(np.abs(np.abs(ytsm)/np.abs(yts)))
        
        # Plot Spectrograms
        plt.figure()
        plt.suptitle(str(rmf[ik]) , fontsize = 12)
        
        plt.subplot(3,2,1) 
        plt.title('Ex_measured')
        plt.pcolormesh(xt/samp,
                       np.log10(xf),
                       np.log10(xts), 
                       cmap = "magma", 
                       vmin=0, 
                       vmax=6)
        plt.xlim(mint,maxt)        
        plt.ylim(-4.2,-2)
        plt.tight_layout()
        plt.ylabel('log10 Freq. [Hz]')
        plt.xlabel('Time [min]')
        plt.colorbar()
        
        plt.subplot(3,2,2) 
        plt.title('Ey_measured')
        plt.pcolormesh(yt/samp,
                       np.log10(yf),
                       np.log10(yts), 
                       cmap = "magma", 
                       vmin=0, 
                       vmax=6)
        plt.xlim(mint,maxt)      
        plt.ylim(-4.2,-2)
        plt.tight_layout()
        plt.ylabel('log10 Freq. [Hz]')
        plt.xlabel('Time [min]')
        plt.colorbar()

        plt.subplot(3,2,3) 
        plt.title('Ex_modelled')
        plt.pcolormesh(xt/samp,
                       np.log10(xf),
                       np.log10(xtsm), 
                       cmap = "magma", 
                       vmin=0,
                       vmax=6)
        plt.xlim(mint,maxt)      
        plt.ylim(-4.2,-2)
        plt.tight_layout()
        plt.ylabel('log10 Freq. [Hz]')
        plt.xlabel('Time [min]')
        plt.colorbar()
        
        plt.subplot(3,2,4)
        plt.title('Ey_modelled')
        plt.pcolormesh(yt/samp,
                       np.log10(yf),
                       np.log10(ytsm), 
                       cmap = "magma", 
                       vmin=0, 
                       vmax=6)
        plt.xlim(mint,maxt)      
        plt.ylim(-4.2,-2)
        plt.tight_layout()
        plt.ylabel('log10 Freq. [Hz]')
        plt.xlabel('Time [min]')
        plt.colorbar()
        
        plt.subplot(3,2,5) 
        plt.title('log10(Ex_meas / Ex_mod)')
        plt.pcolormesh(xt/samp,
                       np.log10(xf),
                       difx, 
                       cmap = "bwr_r", 
                       vmin=-3, 
                       vmax=3)
        plt.xlim(mint,maxt)     
        plt.ylim(-4.2,-2)
        plt.tight_layout()
        plt.ylabel('log10 Freq. [Hz]')
        plt.xlabel('Time [min]')
        plt.colorbar()
               
        plt.subplot(3,2,6)
        plt.title('log10(Ey_meas / Ey_mod)')
        plt.pcolormesh(yt/samp,
                       np.log10(yf),
                       dify, 
                       cmap = "bwr_r", 
                       vmin=-3, 
                       vmax=3)
        plt.xlim(mint,maxt)         
        plt.ylim(-4.2,-2)
        plt.tight_layout()
        plt.ylabel('log10 Freq. [Hz]')
        plt.xlabel('Time [min]')
        plt.colorbar()
        
        plt.show()        
        
    return()        

#######################################################################        
def plot_modelled_vs_meaured(e_site, xx, yy):
    
    """ Statistical comparision between modelled and measured electric fields
		
		Parameters
		-----------
		e_site = Name of the site of interest
         xx = Measured electric field at the site of interest
         yy = Modelled electric field at the site of interest
         
		Returns
		-----------
		
		-----------------------------------------------------------------
    """

    plt.style.use('ggplot')
    for ik in range(0,xx.shape[1]):

        ptxx = np.append(xx[:,ik,0], xx[:,ik,1])
        ptyy = np.append(yy[:,ik,0], yy[:,ik,1])
        
        # Standard deviation of the differences
        plt.figure()
        plt.suptitle("Differences (Measured - Modelled) at " 
                     + str(e_site[ik]), 
                     fontsize = 12)
        dif = ptxx - ptyy
        sns.kdeplot(dif,shade = True, c='r')
        mudif, stddif = norm.fit(dif)
        values = mpatches.Patch(color = 'r', 
                                label = 'Av_diff: %1.1f' %(mudif)
                                + 'std: %1.1f' %(stddif))
        plt.legend(handles = [values])
        plt.tight_layout()
        
        # Measured vs Modelled
        g = sns.jointplot(ptxx, 
                          ptyy, 
                          kind = 'kde', 
                          color = 'steelblue', 
                          space = 0)
        g.plot_joint(plt.scatter, 
                     c = 'k', 
                     s = 10, 
                     linewidth = 0, 
                     marker = 'o', 
                     alpha = 0.2)
        
        g.ax_joint.collections[0].set_alpha(0)
        g.ax_joint.plot(ptxx, ptxx,'r')
        plt.suptitle(str(e_site[ik]), fontsize = 12)
        plt.tight_layout()