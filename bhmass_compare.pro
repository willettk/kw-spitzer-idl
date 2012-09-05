;+
; NAME:
;       
;	BHMASS_COMPARE
;
; PURPOSE:
;
;	Plot black-hole mass tracers from different sources for CSO sample
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
;       Written by K. Willett                Dec 09
;-

; BH masses from bulge luminosity - does not include N5793 or VIIZw485

mbh_bulge = [8.78, 7.61, 8.39, 8.81, 8.81, 8.22, 8.62, 8.73]

; Masses and limits from OIV linewidth

mbh_OIV = [8.16, 7.96, 9.23, 0, 7.79, 0, 8.62, 0]

; Mass from optical line in Rodriguez+06

mbh_optical = [0, 8.74, 0, 0, 0, 0, 0, 0]

obj = ['4C31.04','4C37.11','1146+59','4C12.50','OQ 208','PKS 1413+135','PKS 1718-649','1946+70']

plot, mbh_bulge, mbh_OIV, $
	xr = [7.5, 9.5], $
	yr = [7.5, 9.5], $
	xtitle='MBH_bulge [log M_sun]', $
	ytitle='MBH_line [log M_sun]', $
	psym = symcat(16)

oplot, mbh_bulge, mbh_optical, psym=symcat(14), color=fsc_color("Red")

oplot, indgen(20), linestyle=2

arrow, mbh_bulge[[0,1,4]], mbh_OIV[[0,1,4]], mbh_bulge[[0,1,4]], mbh_OIV[[0,1,4]] - 0.2, /data

xyouts, mbh_bulge + 0.1, mbh_oIV, obj, /data
xyouts, mbh_bulge + 0.1, mbh_optical, obj, /data, color=fsc_color("Red")

end
