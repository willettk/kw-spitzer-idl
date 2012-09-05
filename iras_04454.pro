
;+
; NAME:
;       
;	IRAS_04454
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

henrik_dir = '~/Astronomy/Research/Spitzer/henrik/'

henrik_hr = henrik_dir + '04454-4838.irshigh.xdr'
henrik_lr = henrik_dir + '04454-4838.irslow.xdr'

restore, henrik_hr

wave = sed.wave
flux = sed.flux.jy

plot, wave, flux, $
	xr = [30,40], $
	psym = 10









stop


end
