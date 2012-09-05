;pro pah_comp
;+
; NAME:
;       
;	PAH_COMP
;
; PURPOSE:
;
;	Compare the results of measuring PAH features from PAHFIT and my own routines
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
;       Written by K. Willett                Sep 08
;-


; Read in PAHFIT data

tags = ohmdat('tag')
ntags = n_elements(tags)
pahdir = '~/Astronomy/Research/Spitzer/ohm/pahfit/'

pahfit62 = dblarr(ntags)
pahfit11 = dblarr(ntags)

for i = 0, ntags - 1 do begin

	restore, pahdir+tags[i]+'_fit.sav'

	power_62 = fit.dust_features[2].int_strength		; 10^-26 W/m^2 with /NO_MEGAJANSKY_SR passed into it
	power_112 = fit.dust_features[10].int_strength
	power_113 = fit.dust_features[11].int_strength

	power_11 = power_112 + power_113

	pahfit62[i] = power_62 * 1d-30
	pahfit11[i] = power_11 * 1d-30

endfor

; Read in my data

my62 = double(ohmdat('pah62flux'))					; W/cm^2
my11 = double(ohmdat('pah11flux'))

!p.multi = [0,1,2]

plot, pahfit62, my62, /nodata, $
	xtitle = 'PAHFIT data', $
	ytitle = 'My data'

oplot, pahfit62, my62, psym = symcat(16)
oplot, pahfit11, my11, psym = symcat(9)

x = fillarr(1d-22,0,1d-19)
oplot,x,x

legend, /top, /left, ['6.2 PAH','11.3 PAH'], psym=[16,9]

plot, pahfit62/pahfit11, my62 / my11, psym = symcat(16), yr = [0,2]
hor, 1;w, /data

end
