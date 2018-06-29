# Compute_GeoElectric_Fields

All the steps for computing GeoElectric fields during geomagnetic storms are based on magnetic measurements from magnetic observatories and electromagnetic tensor relationships.

This program was written by Joan Campanya (of TCD), and Sean Blake (of TCD). All python programs are written in Python 2.7.

## Dependencies needed
*numpy*

*pandas*

*matplotlib*

*scipy*

*seaborn*

*geopy*

*multiprocessing* (only to run it in parallel version)



## 1. Main Structure of the program:
The algorithm is divided in three folders:

a.  **in** folder. 

b.  **out** folder. contains sub-products created during the modelling, whithin SECS folder, and the computed electric time series specifying: 1) name of the site, 2) component of the electric field, 3) storm, and 4) approach used to compute the electric fields.

c.  **scr** folder. 

## 1.1 **in** folder

1.1.1. **data** folder with magnetic time series an tensor relationships.

1.1.2. **Observatories.dat** contains the name and coordinates, in degrees, of the magnetic observatories.

1.1.3. **Sites_interest.dat** contains the name and coordinates in degrees of the sites of interest where we want to compute the               electric fields.

The program requires of three main inputs:

a) Time series from the magnetic observatories.
a.1) Time series are located within **data** folder with an additional folder with the name of the geomagnetic storm. On the exaple: **22-23_06_2015**.

a.2) Electric and magnetic time series are located at different folders **B** and **E**. Note that electric time series are only used to compare the results from modelling with the measured electric fields.

a.3) Magnetic and electric time series are one file per type of data and per magnetic observatory.

b) Electromagnetic tensor relationships

b.1) The electromagnetic tensor relationships are located within **data** folder in the **TF** folder.

b.2) The tensor relationships relates time series measured at one site with time series measured at the same site or at a different site. The name of the files with the information of the tensor relationships is orgaised as: 1) Electric *E* or Magnetic *B* component of site 1 + 2) *Name of site 1* + 3) Electric *E* or Magnetic *B* component site 2 + 4) *Name of site 2* + *_s.j*. For example, the MT impedance tensor relating electrics and magnetics from the same site would be *ESITE1BSITE1_s.j*

b.3) The format of the files is the standard format from BIRRP processing code (j file, *http://www.complete-mt-solutions.com/mtnet/data/download_data.html*).

c) Data files with name and coordinates of magnetic observatories (*Observatories.dat*) and sites of interest (*Sites_interest.dat*)



## 1.2 out folder

This folder contains sub-products created during the modelling, whithin **SECS** folder, and the computed electric time series specifying: 1) *name of the site*, 2) *component of the electric field*, 3) *storm*, and 4) *approach used to compute the electric fields*.


## 1.3 **scr** folder

This filder contains the main programs and an input file with the main parameters:

1.3.1.1. **EM_modelling.py**: Model electric and magnetic time series based on approaches presented by Campanya et al., (2018).

1.3.1.2. **functions_EM_modelling.py** contains all the functions requested by the program EM_modelling.py.

1.3.2.1. **EM_plotting.py**: Compares measured and modelled time series in time and frequency domain.

1.3.2.2. **functions_EM_plotting.py** contains all the functions requested by the program EM_plotting.py.

1.3.4.1. **SECS_interpolation.py** interpolates the magnetic fields between observatories using spherical elementary current systems.

1.3.4.2. **secs_pre.py** Contains functions and vairables requested for SECS_interpolation.py.

1.3.5. **inputs.py** contains all the inputs that need to be specified for modelling and plotting IEFs. 

**inputs.py** is the only file that needs to be modified when all the necessary inputs (electromagnetic time series and electromagnetic tensor relationships) are at the selected folders (see details below). On the presented example, only the *main_path* in **inputs.py** will need to be modified and the example should run, first executing **EM_modelling.py** and then **EM_plotting.py**.
