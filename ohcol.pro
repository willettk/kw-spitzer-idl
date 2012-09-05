function ohcol, ew
;+
; NAME:
;       
;	OHCOL
;
; PURPOSE:
;
;	Compute OH column density for 34.6 um absorption feature
;
; INPUTS:
;
;	EW - 		equivalent width (in microns)
;
; OUTPUTS:
;
;	Returns the column density in cm^-2
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> print,ohcol(0.01)
;
; NOTES:
;
;	General eqn. for optically thin absorption from Lequeux, pp. 56
;
; REVISION HISTORY
;       Written by K. Willett                Apr 09
;-

c = 2.99792458d10

; 35 um data

lambda = 34.616d-3
gl = 4.
gu = 6.
a_ul = 1.74d-2

nl_35 = (ew * 8d * !dpi * c * gl) / (a_ul * lambda^4 * gu)

; 53 um data

lambda = 53.29d-3
gl = 4.
gu = 4.
a_ul = 4.56d-2

nl_53 = (ew * 8d * !dpi * c * gl) / (a_ul * lambda^4 * gu)

return, nl_35

end
