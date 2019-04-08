#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
 Algorithm to plot the modelled geoelectric fields and compare them 
 to measured data.

 Original author: Joan Campanya i Llovet (Trinity College Dublin) [joan.campanya@tcd.ie]
 08 - 2018, Dublin, Ireland.

 Inputs 
   1) Measured time series
   2) Modelled time series

   Input values are specified in the "inputs.py" file in the "scr" folder. 

 Outputs (in "out" folder)
   1) Comparison between measured and modelled time series in the time
      domain, including uncertanties.
   2) Comparison between measured and modelled time series in the frequency
      domain (spectrograms)
   3) Statistical analysis comparing modelled time series vs measured time series

Altgorithm related to Campanya et al. 2018 publication in Space Weather AGU Journal
"""
from scr import functions_EM_plotting as femp
from constants import *
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

#######################################################################
# 1) Read the input values

# exec(open('inputs.py').read())
#######################################################################
# 2) Read coordinates of the sites of interest
e_site, e_lat, e_lon = femp.read_co(in_path + str(sit_f))

#######################################################################
# 3) Deffine main variables
m_e_field = np.zeros([length,len(e_site),2])
comp_Efield = np.zeros([length,len(e_site), 2])
std_Efield = np.zeros([length,len(e_site), 2])

reg_e = np.zeros([length,len(e_site), 2])
std_reg_e = np.zeros([length,len(e_site), 2])

loc_e = np.zeros([length,len(e_site), 2])
std_loc_e = np.zeros([length,len(e_site), 2])

#######################################################################
# 4) Read measured and computed electric time series
for v, ip in enumerate(e_site):
    m_e_field[:,v,:] = femp.read_electrics(e_path, 
                                             ip, 
                                             samp, 
                                             hi, 
                                             low
                                             )
    
    
    comp_Efield[:,v,:], std_Efield[:,v,:] = femp.read_comp_elc(out_path, 
                                                                ip, 
                                                                samp, 
                                                                storm, 
                                                                mode
                                                                )

#######################################################################
# 5) Plot the time series and compute fit of the data


# 5.1) Comparision in time domain including errors
femp.polt_time_series(e_site, 
                      m_e_field, 
                      comp_Efield, 
                      std_Efield,
                      mint, 
                      maxt,
                      0
                      )

# Save measured E_fields with selected periods of range
for vi, ppi in enumerate(e_site):
    m_ex_path_id = (out_path 
                 + str(ppi) 
                 + '_Ex_' 
                 + str(storm) 
                 + '_measured'
                 + '.dat')
    
    m_ey_path_id = (out_path 
                 + str(ppi) 
                 + '_Ey_' 
                 + str(storm) 
                 + '_measured'
                 + '.dat')        
    
    ex = np.array([m_e_field[:,vi,0], std_Efield[:,vi,0]])
    ey = np.array([m_e_field[:,vi,1], std_Efield[:,vi,1]])
        
    np.savetxt(m_ex_path_id, ex.T, fmt=['%15.5f' , '%15.5f'])  
    np.savetxt(m_ey_path_id, ey.T, fmt=['%15.5f', '%15.5f'])   




## 5.2) Comparision using spectrograms
#femp.polt_spectogram(e_site, 
#                     m_e_field, 
#                     comp_Efield, 
#                     samp, 
#                     mint, 
#                     maxt
#                     )
##
## 5.3) Grafic measured vs modelled and differences evaluation 
#femp.plot_modelled_vs_meaured(e_site, 
#                              m_e_field, 
#                              comp_Efield)

# 5.4) Local and regional signals (Approach_#2)
if mode == 2:

    for v, ip in enumerate(e_site):
        reg_e[:,v,:], std_reg_e[:,v,:] = femp.read_reg_elc(out_path, 
                                                           ip, 
                                                           samp, 
                                                           storm, 
                                                           mode
                                                           )
        
        
        loc_e[:,v,:], std_loc_e[:,v,:] = femp.read_loc_elc(out_path, 
                                                           ip, 
                                                           samp, 
                                                           storm, 
                                                           mode
                                                           )
    
    ## Plot regional signal
    femp.polt_time_series(e_site, 
                          m_e_field, 
                          reg_e, 
                          std_reg_e, 
                          mint, 
                          maxt,
                          1
                          )

    ## Plot local signal
    femp.polt_time_series(e_site, 
                          m_e_field, 
                          loc_e, 
                          std_loc_e, 
                          mint, 
                          maxt,
                          2
                          )
