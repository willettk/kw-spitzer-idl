function neV_lum, fname
;+
; NAME:
;       
;
;
; PURPOSE:
;
;	
;
; INPUTS:
;
;
;
; OUTPUTS:
;
;
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                
;-

flux = linelim(fname,'neV')

taglist = csodat('tag')
tagind = where(taglist eq fname)
dllist = csodat('dl')
dl = dllist(tagind)

mpc2cm = 3.086d24
lsun = 3.862d33
lum = 4d * !dpi * (dl * mpc2cm)^2 * flux * 1d7 * 1d-21


return, alog10(lum / lsun)



end
