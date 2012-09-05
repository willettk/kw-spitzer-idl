; pro kmpl
;+
; NAME:
;       
;	Compute the Kaplan-Meier product-limit estimator for survival analysis
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
;       Written by K. Willett                Nov 08
;-

; Data

x = [30,24,11,19,27,11,24,28]
lim = [1,0,1,0,1,1,1,0]			; 0 indicates that the value is an UPPER limit

x_prime = max(x) - x

sortind = sort(x_prime)

x_prime = x_prime[sortind]
lim_prime = lim[sortind]

x_prime_lim = x_prime[where(lim_prime eq 1)]

dj = intarr(n_elements(x_prime_lim))
nj = intarr(n_elements(x_prime_lim))
fhat_xl = dblarr(n_elements(x_prime_lim))

fhat_xl_temp = 1.

for i = 0, n_elements(x_prime_lim) - 1 do begin

	dj_temp = where(x_prime_lim eq x_prime_lim[i], dcount)
	dj[i] = dcount

	nj_temp = where(x_prime ge x_prime_lim[i], ncount)
	nj[i] = ncount

	if x_prime_lim[i] eq 0 then fhat_xl_temp = 1. else fhat_xl_temp = fhat_xl_temp * (1d - float(dj[i-1]) / float(nj[i-1]))
	fhat_xl[i] = fhat_xl_temp

	print, nj[i], dj[i], 1d - float(dj[i])/float(nj[i]), fhat_xl[i]

endfor


stop
end
