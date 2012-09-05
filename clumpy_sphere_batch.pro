;+
; NAME:
;       
;	CLUMPY_SPHERE_BATCH
;
; PURPOSE:
;
;	Run CLUMPY on spherical shell models
;
; INPUTS:
;
;	OBJ - 		tag of galaxy to run models on (eg, 'arch004')
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
;	IDL> clumpy_sphere_batch, 'arch004'
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
;	Generalized to all objects - Jan 10
;-

pro clumpy_sphere_batch, obj, verbose = verbose

;;;;;;;;;;;;;;;;;;;;

tag, obj, dirtag

; Parameter space (coarsely sampled - only models downloaded to local machine)

Y = [5, 10, 30]
N0 = [1, 2, 4, 6, 8, 10, 12, 15, 20]
q = [0.0, 1.0, 2.0, 3.0]
tauv = [10, 30, 60, 80, 100, 200, 300, 500]

errorgrid_sphere = fltarr(n_elements(Y), n_elements(N0), n_elements(q), n_elements(tauv))

timestart = systime(1)

for j = 0, n_elements(Y) - 1 do begin
	for k = 0, n_elements(N0) - 1 do begin
		for l = 0, n_elements(q) - 1 do begin
			for m = 0, n_elements(tauv) - 1 do begin

				clumpy_sphere, Y[j], N0[k], q[l], tauv[m], $
					errormerit, $
					obj=obj, $
					/noplot, $
					/quiet

				errorgrid_sphere[j,k,l,m] = errormerit

				if keyword_set(verbose) then begin
					print, 'Y = ', Y[j]
					print, 'N0 = ', N0[k]
					print, 'q = ', q[l]
					print, 'tau_V = ', tauv[m]
				endif

			endfor
		endfor
	endfor
endfor

timeend = systime(1)

save, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/sphere_sav/'+obj+'_sphere_grid.sav', minvalues, errorgrid_sphere

print,''
print, 'Time elapsed: ',(timeend - timestart),' seconds'

; Find best-fit model

minerror = where(errorgrid_sphere eq min(errorgrid_sphere))
minvalues = errorgrid_parse(minerror,/sphere)

print,''
print,'Sphere best fit parameters: '
print,''
print, 'Y = ', minvalues[0]
print, 'N0 = ', minvalues[1]
print, 'q = ', minvalues[2]
print, 'tau_V = ', minvalues[3]

clumpy_sphere, minvalues[0], minvalues[1], minvalues[2], minvalues[3], $
	obj=obj, $
	/quiet

end
