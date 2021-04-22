# Quasi-Bessel-Beam
Calculation of quasi-Bessel beams after an ideal axicon

Python module for the calculation of quasi-Bessel beams (QBB) generated via an ideal Gaussian beam and an ideal axicon as per
P. Wu et al., Theoretical analysis of a quasi-Bessel beam for laser ablation, Photon. Res. / Vol. 2, No. 3 / June 2014
http://dx.doi.org/10.1364/PRJ.2.000082


## How to use
The main use is the calculation of the radial intensity profile of the QBB at a position rho away from the optical axis and an on-axis distance z behind the axicon via the function

`radial_intensity(rho, z, w = 1, I_0 = 1, n = 1.5, alpha = 5, wavelength = 1E-6)`

The module also allows to **approximate** the depth of field (DOF) and the size of the central peak via the functions

    dof(w, n = 1.5, alpha = 5)

    beam_diameter(wavelength = 1E-6, n = 1.5, alpha = 5)
  
