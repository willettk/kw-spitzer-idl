
;+
; NAME:
;       
;	IRAS12112
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
;       Written by K. Willett                Dec 08
;-

restore,'~/Astronomy/Research/Spitzer/12112+0305.irslow.xdr'

!p.multi=[0,1,1]
plot,sed.wave,sed.flux.jy, $
	/xlog, $
	xr=[5,35], /xstyle, $
	/ylog, $
	xtitle='Wavelength', $
	ytitle='Flux [Jy]', $
	title = 'IRAS 12112+0305', $
	charsize = 1.5, $
	psym = 10


restore,'~/Astronomy/Research/Spitzer/archived/data/structures/arch018.sav'

oplot,sed.wave_lr,sed.flux_lr, $
	psym = 10, $
	color=fsc_color("Red")

oplot,sed.wave_lr,sed.flux_lr * 1.2, $
	psym = 10, $
	color=fsc_color("Blue")


stop
end
