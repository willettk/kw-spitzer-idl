
;+
; NAME:
;       
;	PAH62ICELIST
;
; PURPOSE:
;
;	Save list of targets to use the water ice-subtracted continuum for the 6.2 um PAH measurement
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
;       Written by K. Willett                Aug 08
;-

icelist62 = [$
	'mega004','mega008','mega017','mega023','mega027','mega028','mega033', $		; mega017 seems to show H2O ice, but it does
	'arch004','arch005','arch007','arch008','arch010','arch013','arch014', $		; not affect the EW measurement (PAH spline = sil spline)
	'arch017','arch018','arch024','arch026','arch029','arch030','arch031', $
	'arch032','arch033','arch035','arch036','arch039','arch040','arch048', $
	'control004','control034','control035', $
	'cso002']

save, filename='~/Astronomy/Research/Spitzer/icelist62.sav', $
	icelist62

end
