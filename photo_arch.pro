;pro photo_arch
;+
; NAME:
;       PHOTO_ARCH
;
; PURPOSE:
;	Store the IR photometry of archived Spitzer OHMs and save to file
;
; INPUTS:
;
; OUTPUTS:
;
; KEYWORDS:
;
; EXAMPLE:
;	IDL> .r photo_arch
;
; REQUIRES:
;
; NOTES:
;
;
; REVISION HISTORY
;       Written by K. Willett                Dec 2007
;-

fnames = archdat('tag')
nf = n_elements(fnames)

iras = fltarr(nf,4)			; IRAS {12, 25, 60, 100} um [Jy]
iras_err = fltarr(nf,4)			; IRAS error [%] 

iso = fltarr(nf,4)			; ISO  {90, 150, 170, 200} um [Jy]
iso_err = fltarr(nf,4)			; ISO  error [%]

; arch001 
i = 0
iras(i,*)	= [5.118d-1,1.211,2.243,2.634]
iras_err(i,*) 	= [9,10,10,10]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch004 
i = 1
iras(i,*)	= [0,0.233,2.638,3.833]
iras_err(i,*) 	= [0,14,9,9]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch005 
i = 2
iras(i,*)	= [0,4.513d-1,5.685d,5.246]
iras_err(i,*) 	= [0,6,4,6]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch006 
i = 3
iras(i,*)	= [5.010d-1,8.129d-1,1.445d1,2.926d1]
iras_err(i,*) 	= [8,6,5,5]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch008 
i = 4
iras(i,*)	= [2.499d-1,1.034,1.154d1,2.023d1]
iras_err(i,*) 	= [15,6,7,5]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch010 
i = 5
iras(i,*)	= [2.945d-1,1.105,8.901,7.981]
iras_err(i,*) 	= [10,8,8,6]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch013 
i = 6
iras(i,*)	= [0,2.351d-1,2.281,1.816]
iras_err(i,*) 	= [0,25,6,9]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch016 
i = 7
iras(i,*)	= [0,3.735d-1,1.761,1.776]
iras_err(i,*) 	= [0,28,7,13]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch018 
i = 8
iras(i,*)	= [0,5.093d-1,8.503,9.976]
iras_err(i,*) 	= [0,13,6,8]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch019
i = 9
iras(i,*)	= [9.347d-1,9.320,4.068d1,3.280d1]
iras_err(i,*) 	= [6,7,9,11]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch020 
i = 10
iras(i,*)	= [1.872,8.662,3.199d1,3.029d1]
iras_err(i,*) 	= [5,5,5,4]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch021 
i = 11
iras(i,*)	= [3.431d-1,9.638d-1,1.113d1,2.086d1]
iras_err(i,*) 	= [10,11,6,6]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch022 
i = 12
iras(i,*)	= [0,1.272,1.793d1,1.813d1]
iras_err(i,*) 	= [0,8,8,6]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch023 
i = 13
iras(i,*)	= [2.621d-1,4.032d-1,1.174,7.134d-1]
iras_err(i,*) 	= [14,21,7,20]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch024 
i = 14
iras(i,*)	= [2.352d-1,2.282,2.174d1,2.138d1]
iras_err(i,*) 	= [10,4,4,4]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch025 
i = 15
iras(i,*)	= [0,6.695d-1,1.916,2.060]
iras_err(i,*) 	= [0,10,11,9]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch026 
i = 16
iras(i,*)	= [0,0,1.447,1.821]
iras_err(i,*) 	= [0,0,6,9]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch028 
i = 17
iras(i,*)	= [0.12,1.323,7.286,5.907]
iras_err(i,*) 	= [27,5,5,5]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch029 
i = 18
iras(i,*)	= [4.837d-1,7.907,1.038d2,1.124d2]
iras_err(i,*) 	= [5,5,4,3]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch031 
i = 19
iras(i,*)	= [2.361d-1,1.084,7.90,1.252d1]
iras_err(i,*) 	= [11,6,7,9]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch032 
i = 20
iras(i,*)	= [1.953d-1,1.658,3.114d1,3.490d1]
iras_err(i,*) 	= [24,6,6,6]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch033 
i = 21
iras(i,*)	= [0,3.433d-1,5.228,5.165]
iras_err(i,*) 	= [0,11,6,7]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch036
i = 22
iras(i,*)	= [0,5.495d-1,5.436,4.453]
iras_err(i,*) 	= [0,13,7,8]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; arch040 
i = 23
iras(i,*)	= [0,8.083d-1,7.087,8.363]
iras_err(i,*) 	= [0,7,10,6]
iso(i,*)	= [0,0,0,0]
iso_err(i,*)	= [0,0,0,0]

; Set errors to 20% if not given

iras_zeroes = where(iras_err eq 0)
iras_err(iras_zeroes) = 20.
iso_zeroes = where(iso_err eq 0)
iso_err(iso_zeroes) = 20.

iras_err = 1d-2 * iras_err * iras
iso_err  = 1d-2 * iso_err * iso

; Save to IDL file

save, filename = '~/Astronomy/Research/Spitzer/archived/data/idl_sav/photo_arch.sav', iras, iras_err, iso, iso_err

end
