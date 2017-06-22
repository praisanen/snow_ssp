# snow_ssp
Snow single-scattering parameterization

This repository contains the following files:

1) snow_ssp.f90
= A fortran90 program for computing the following single-scattering properties of
  snow in the wavelength range 0.199-2.7 µm:
   - extinction efficiency (assumed to be exactly 2)
   - single-scattering co-albedo = 1 - single-scattering albedo
   - asymmetry parameter
   - scattering phase function

The input data needed for the parameterization include 
   - real part of ice refractive index
   - imaginary part of ice refractive index
   - size parameter defined with respect to the snow grain volume-to-projected area equivalent radius
   
For more information, please see this publication: 
   
  P. Räisänen, A. Kokhanovsky, G. Guyot, O. Jourdan and T. Nousiainen, 2015:
  Parameterization of single-scattering properties of snow, The Cryosphere, 
  9, 1277-1301, doi:10.5194/tc-9-1277-2015.
 
2) files_for_NorESM_or_CESM 
 
= Directory containing files that may be used for testing the impact of snow grain nonsphericity in
  the NorESM and CESM earth system models (or any model using the SNICAR scheme for snow on land 
  (a part of the CLM land-surface model) and the CICE4 sea ice model)
 
3) gpl.txt 
= Gnu General Public License, version 3

4) README.md
= this file
