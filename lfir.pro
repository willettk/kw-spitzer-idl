function lfir, f60, f100, dl
;+
; NAME:
;       
;	LFIR
;
; PURPOSE:
;
;	Compute the far-IR total luminosity based on 60 and 100 um photometry	
;
; INPUTS:
;
;	f60 - 		60 um flux [Jy]
;
;	f100 - 		100 um flux [Jy]
;
;	dl - 		luminosity distance [Mpc]
;
; OUTPUTS:
;
;	LFIR - 		far-IR luminosity in log(L_sun)
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> print,lfir(2.56, 3.45, 104)
;
; NOTES:
;
; Formula comes from Sanders & Mirabel (1996), ARAA
;
; REVISION HISTORY
;       Written by K. Willett                
;-


mpc2cm = 3.086d24 ; cm/Mpc
C = 1.6	; scale factor correcting for extrapolated flux longward of 100 um IRAS filter. Typical range is from
	; 1.4 - 1.8
mks2cgs = 1d3	; 1 W m^-2 = 10^3 erg s^-1 cm^-2
lsun = 3.862d33 ; erg s^-1

fir_flux = 1.26d-14 * (2.58 * f60 + f100)	; W m^-2

fir_lum = (4d * !dpi * (dl*mpc2cm)^2 * C * fir_flux * mks2cgs) / lsun

; Note that Darling paper II does NOT use the fudge factor C (off from my values by a factor of 1.6, then)

return, alog10(fir_lum)

end
