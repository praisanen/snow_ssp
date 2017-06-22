
This directory contains files that may be used for testing the effect of
non-spherical snow grains in the surface albedo calculation in either 
the Norwegian Earth System Model (NorESM) or in the Community Earth System 
Model (CESM), or in principle, any other model employing the SNICAR scheme
for snow on land and the CICE4 sea ice model.

The single-scattering properties (aka. optical properties, or inherent optical properties) 
of non-spherical snow grains are based on the Optimized habit combination (OHC) as described 
in Räisänen et al. (2015) and the method of computing band-mean values is described in the 
Appendix of Räisänen et al. (2017); see the References section below.

**************************************************
A) SINGLE-SCATTERING PROPERTIES OF SNOW OVER LAND
**************************************************

Three versions of the "snicar_optics_5bnd_*.nc" file are provided.

A.1) snicar_optics_5bnd_c090915.nc
----------------------------------
This is the default version in SNICAR / CLM4 containing the single-scattering
properties of spherical snow grains (for 1471 values of effective radius
ranging from 30 to 1500 µm) and several aerosol species, each defined
separately for the five shortwave spectral bands in SNICAR. This file has
been created by Mark Flanner, and it was obtained from the CESM public 
Subversion input data repository, directory
https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/snicardata/ .

A.2) snicar_optics_5bnd_spheres_PRaisanen.nc
--------------------------------------------
In this version, the single-sacttering properties of snow grains are still based
on spheres, but they have been recomputed as explained in the Appendix
of Räisänen et al. (2017). The aerosol single-scattering properties are unchanged.

A.3) snicar_optics_5bnd_OHC_PRaisanen.nc
----------------------------------------
In this version, the single-scattering properties of snow grains are based on
the "Optimized habit combination" (a mixture of three non-spherical shapes)
introduced in Räisänen et al. (2015), with spectral averaging handled
as explained in the Appendix of Räisänen et al. (2017). The aerosol single-scattering
properties are unchanged.

*****************************************************
B) SINGLE-SCATTERING PROPERTIES OF SNOW OVER SEA ICE
*****************************************************

Two versions of the CICE4 fortran file "ice_shortwave.F90" are provided.

B.1) ice_shortwave.F90 
--------------------------
This is the original version of "ice_shortwave.F90" in CICE4,
provided here just for reference. Spherical snow grains are assumed. This
version is from NorESM, but it is (virtually) identical to the CESM1.0 version
available (e.g.) here:
http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesmBbrowser/html_code/cice/ice_shortwave.F90.html

B.2) ice_shortwave_PRaisanen.F90
------------------------------
This is a modified version which supports the use of either spherical or 
non-spherical snow grains, the latter activated with the compiler directive 
"NONSPH_SNOW". Note that even in the spherical case, the snow single-scattering
properties (that is, arrays Qs_tab, ws_tab and gs_tab) differ slightly
from those in the original version, due to differences in assumed
ice refractive index and spectral averaging technique (see the Appendix
of Räisänen et al. 2017).

*************
C) REFERENCES
*************

Räisänen, P., Kokhanovsky, A., Guyot, G., Jourdan, O. and Nousiainen, T., 2015:
Parameterization of single-scattering properties of snow, The Cryosphere, 9, 
1277-1301, doi:10.5194/tc-9-1277-2015.

Räisänen, P., Makkonen. R., Kirkevåg, A. and Debernard, J. B., 2017:
Effect of snow grain shape on climate simulations: Sensitivity tests
with the Norwegian Earth System Model (The Cryosphere Discussions, submitted)
