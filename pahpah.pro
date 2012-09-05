
;+
; NAME:
;       
;	PAHPAH
;
; PURPOSE:
;
;	Plot luminosities of 6.2 and 11.3 um PAH features against each other
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

ohm62 = ohmdat('pah62lum')
ohm11 = ohmdat('pah11lum')

arch62 = archdat('pah62lum')
arch11 = archdat('pah11lum')

con62 = condat('pah62lum')
con11 = condat('pah11lum')

cso62 = csodat('pah62lum')
cso11 = csodat('pah11lum')

; Kluge for limits on cso spline-fit PAH data

cso62[5 - 1] = 9.05
cso11[5 - 1] = 8.76

cso62[7 - 1] = 7.48
cso62[8 - 1] = 7.21

plot, ohm62, ohm11, $
	/nodata, $
	xtitle='6.2 um PAH luminosity', $
	ytitle='11.3 um PAH luminosity', $
	xrange=[5,10], /xstyle, $
	yrange=[5,10], /ystyle

oplot, ohm62, ohm11, psym = symcat(16)
oplot, arch62, arch11, psym = symcat(16)
oplot, con62, con11, psym = symcat(16)
oplot, cso62, cso11, psym = symcat(15), color=fsc_color("Blue")

; Limits

ind=[4,6,7]
arrow, cso62[ind], cso11[ind], cso62[ind] - 0.5, cso11[ind], color=fsc_color("Blue"),/data
arrow, cso62[4], cso11[4], cso62[4], cso11[4] - 0.5, color=fsc_color("Blue"),/data

oplot, findgen(20), linestyle=1

end
