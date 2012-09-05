;+
; NAME:
;       
;	CLUMPY_BATCH
;
; PURPOSE:
;
;	Run all downloaded CLUMPY models and find best fit
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
;	IDL> clumpy_batch, 'cso005'
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
;-

pro clumpy_batch, obj

;;;;;;;;;;;;;;;;;;;;

; Parameter space

sig  = [15, 30, 45, 60, 75]
Y    = [5, 10, 30, 100, 200]
N0   = [1, 5, 10, 15, 20]
q    = [0.0, 1.0, 2.0, 3.0]
tauv = [10, 30, 60, 80, 100, 200, 300, 500]
inc  = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

;;;;;;;;;;;;;;;;;;;;

errorgrid = fltarr(n_elements(sig), n_elements(Y), n_elements(N0), n_elements(q), n_elements(tauv), n_elements(inc))

timestart = systime(1)

for i = 0, n_elements(sig) - 1 do begin
	for j = 0, n_elements(Y) - 1 do begin
		for k = 0, n_elements(N0) - 1 do begin
			for l = 0, n_elements(q) - 1 do begin
				for m = 0, n_elements(tauv) - 1 do begin
					for n = 0, n_elements(inc) - 1 do begin

					clumpy, sig[i], Y[j], N0[k], q[l], tauv[m], inc[n], $
						errormerit, obj = obj, $
						/noplot, /quiet
					
					errorgrid[i,j,k,l,m,n] = errormerit
					
					print, 'sigma = ', sig[i], ' deg'
					print, 'Y     = ', Y[j]
					print, 'N0    = ', N0[k]
					print, 'q     = ', string(q[l],format='(f9.2)')
					print, 'tau_V = ', tauv[m]
					print, 'inc   = ', inc[n], ' deg'
					
					endfor
				endfor
			endfor
		endfor
	endfor
endfor

timeend = systime(1)

tag, obj, dirtag
save, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/torus_sav/'+obj+'_torus_grid.sav', errorgrid

print,''
print, 'Time elapsed: ',fix((timeend - timestart)/6.)/10.,' minutes'

; Find best-fit model

minerror = where(errorgrid eq min(errorgrid))
minvalues = errorgrid_parse(minerror)

print,''
print,'Torus best fit parameters: '
print,''
print, 'sigma = ', minvalues[0],' deg'
print, 'Y     = ', minvalues[1]
print, 'N0    = ', minvalues[2]
print, 'q     = ', minvalues[3]
print, 'tau_V = ', minvalues[4]
print, 'inc   = ', minvalues[5],' deg'

clumpy, minvalues[0], minvalues[1], minvalues[2], $
	minvalues[3], minvalues[4], minvalues[5], $
	obj=obj, $
	/quiet

end
