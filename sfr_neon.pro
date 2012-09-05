function sfr_neon, fname
;+
; NAME:
;       
;	SFR_NEON
;
; PURPOSE:
;
;	Compute the star formation rate using the NeII and NeIII fluxes and the Ho/Keto relation
;
; INPUTS:
;
;	FNAME - 	string giving tag of object
;
; OUTPUTS:
;
;	SFR - 		star formation rate in solar masses/year
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> results = sfr_neon('cso001')
;
; NOTES:
;
;	Should be adapted to work on all objects, w/error catching if lines are not measured.
;
; REVISION HISTORY
;       Written by K. Willett                Jun 08
;-

restore,'~/Astronomy/Research/Spitzer/cso/linedata/neII.sav'
ne2ind = where(line.tag eq fname)
ne2flux = line.flux[ne2ind]

restore,'~/Astronomy/Research/Spitzer/cso/linedata/neIII.sav'
ne3ind = where(line.tag eq fname)
ne3flux = line.flux[ne3ind]

taglist = csodat('tag')
tagind = where(taglist eq fname)
dllist = csodat('dl')
dl = dllist(tagind)

mpc2cm = 3.086d24
ne2lum = 4d * !dpi * (dl * mpc2cm)^2 * ne2flux * 1d7
ne3lum = 4d * !dpi * (dl * mpc2cm)^2 * ne3flux * 1d7

nelum = ne2lum + ne3lum

sfr = 4.34d-41 * nelum / (0.75)

return,sfr

end
