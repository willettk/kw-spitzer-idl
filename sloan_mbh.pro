;+
; NAME:
;       
;	SLOAN_MBH
;
; PURPOSE:
;
;	Convert SDSS magnitudes of galaxies into estimates for bulge luminosity and black hole mass
;
; INPUTS:
;
;	G, R - 			apparent magnitudes in g, r, photometric bands from SDSS DR7
;	
;	DISTANCE - 		distance to object [Mpc]
;
; OUTPUTS:
;
;	LOG_MBH_MSUN - 		log of the black hole mass in units of M_sun
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> print, sloan_mbh
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
;	Added error computation		Dec 09
;-

function sloan_mbh, g=g, r=r, d=d, extinction = extinction, mv = mv, sigma_MV = sigma_MV

if not keyword_set(mv) then begin

	; Convert SDSS g, r magnitudes to V-band magnitude (Jester et al. 2005 - assume u-g > 0)
	
	apparent_v = g - 0.59 * (g-r) - 0.01

endif else apparent_v = mv

; Convert apparent to absolute magnitude using distance modulus (assume distance in Mpc)

	; Correct for Galactic extinction
	
	if not keyword_set(extinction) then A_V = 0 else A_V = extinction

absolute_v = apparent_v - 5 * alog10(d * 1d5) - A_V

; Convert V-band magnitude to luminosity (Bentz et al. 2009)

log_lv_lsun = 0.4 * (-1 * absolute_v + 4.83)

; Compute log (M_BH / M_sun) from relation in Bentz et al. (2009)

	; k, alpha use BCES method with AGN only
	
	k = -0.02
	alpha = 0.80

log_mbh_msun = k + alpha * log_lv_lsun + 8. - 10. * alpha

sigma_k = 0.06
sigma_alpha = 0.09
if not keyword_set(sigma_MV) then sigma_MV = 0.2

dMdalpha = alog10(-0.4 * absolute_v + 1.932)

dMdMV = -0.4 * alpha / (alog(10) * (-0.4 * absolute_v + 1.932))

sigma_mbh = sqrt(sigma_k^2 + sigma_alpha^2 * dMdalpha^2 + sigma_MV^2 * dMdMV^2)

;print, apparent_v
;print, absolute_v
;print, log_lv_lsun

mbh = [log_mbh_msun, sigma_mbh]

return, mbh

end
