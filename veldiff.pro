function veldiff, varz, refz
;+
; NAME:
;       
;	VELDIFF
;
; PURPOSE:
;
;	Compute the velocity difference between two redshifts
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
;       Written by K. Willett                Apr 08
;-

; Transform the desired redshift into a velocity

c = 299792.458d	; km/s

var_vel = c * (varz^2 + 2d*varz) / (varz^2 + 2d*varz + 2d)

ref_vel = c * (refz^2 + 2d*refz) / (refz^2 + 2d*refz + 2d)



return, var_vel - ref_vel



end
