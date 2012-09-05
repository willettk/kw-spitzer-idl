;+
; NAME:
;       
;	ERRORGRID_PARSE
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
;	SPHERE - 	for clumpy, spherical distribution models run by CLUMPY
;
;	DUSTY - 	for smooth spherical models run in DUSTY
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
;       Written by K. Willett                Nov 09
;	Added SPHERE keyword - KW, Nov 09
;	Returns parameter values, not indices - Jan 10
; 	Added full range of Y for new models - Mar 10
;-

function errorgrid_parse, index, sphere = sphere, dusty = dusty, alltau = alltau

if keyword_set(sphere) then begin

	array_Y = [5, 10, 30]
	array_N0 = [1, 2, 4, 6, 8, 10, 12, 15, 20]
	array_q = [0.0, 1.0, 2.0, 3.0]
	array_tauv = [10, 30, 60, 80, 100, 200, 300, 500]

	Y_ind = index mod 3
	N0_ind = (index / 3) mod 9
	q_ind = (index / 3 / 9) mod 4
	tauv_ind = (index / 3 / 9 / 4) mod 8

	retarr = [array_Y[Y_ind], array_N0[N0_ind], array_q[q_ind], array_tauv[tauv_ind]]

endif else if keyword_set(dusty) then begin

	;Y = [2, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 400, 500, 750, 1000]
	Y = [2, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 250,300,350,400,450,500,600,650,700,750,800,850,900,950,1000]
	q = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

	if keyword_set(alltau) then taumod = 180 else taumod = 90

	ny = n_elements(y)
	nq = n_elements(q)

	Y_ind = index mod ny
	q_ind = (index / ny) mod nq
	s_ind = (index / ny / nq) mod taumod

	retarr = [Y[Y_ind], q[q_ind], s_ind]

endif else begin

	array_sig = [15, 30, 45, 60, 75]
	array_Y = [5, 10, 30, 100, 200]
	array_N0 = [1, 5, 10, 15, 20]
	array_q = [0.0, 1.0, 2.0, 3.0]
	array_tauv = [10, 30, 60, 80, 100, 200, 300, 500]
	array_inc = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

	sig_ind = index mod 5
	Y_ind = (index / 5) mod 5
	N0_ind = (index / 5 / 5) mod 5
	q_ind = (index / 5 / 5 / 5) mod 4
	tauv_ind = (index / 5 / 5 / 5 / 4) mod 8
	inc_ind = (index / 5 / 5 / 5 / 4 / 8) mod 10

	retarr = [array_sig[sig_ind], array_Y[Y_ind], array_N0[N0_ind], array_q[q_ind], array_tauv[tauv_ind], array_inc[inc_ind]]

endelse

return, retarr

end
