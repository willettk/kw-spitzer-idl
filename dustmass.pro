function dustmass, s60, dl, dtemp
;+
; NAME:
;       
; 	DUSTMASS
;
; PURPOSE:
;
;	Compute the mass of the warm dust sampled by the blackbody peaking near 60 um.
;
; INPUTS:
;
;	S60 - 60 um IRAS flux [Jy]
;
;	DL - luminosity distance of the target [Mpc]
;
;	DTEMP - measured dust temperature [K]
;
; OUTPUTS:
;
;	Returns log (M_dust [M_sun])
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> print, dustmass(0.66, 500, 50)
;
;	Gives the dust mass of an object with S_60 = 0.66 Jy, D_L = 500 Mpc, and T_dust = 50 K.
;
; NOTES:
;
;	Dust mass equation is based on Eq. 3, Yang et al (2007), ApJ, 660:1198
;
;	Dust absorption coefficient assumes a = 0.1 um, rho = 3.3 g cm^-2 from Draine & Lee (1984)
;
; REVISION HISTORY
;       Written by K. Willett                Apr 08
;-

; Physical constants (cgs)

c = 299792.458d * 1d5
kb = 1.381d-16
h = 6.626d-27
lambda = 60 * 1d-4
nu = c / lambda

; Conversions

mpc2cm = 3.086d24
jy2cgs = 1d-23
msun = 1.989d33

q60 = 3.5d-3
a = 0.1 * 1d-4
rho = 3.3

; Equation

s60 = s60 * jy2cgs
dl = dl * mpc2cm

k60 = 3d * q60 / (4d * a * rho)

b60 = ((2d * h * nu^3) / c^2) / (exp(h * nu / (kb * dtemp)) - 1d)

mdust = s60 * dl^2 / (k60 * b60)

return, alog10(mdust / msun)

end
