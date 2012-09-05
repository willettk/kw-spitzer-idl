;+
; NAME:
;       
;	CLUMPY_ANALYSIS
;
; PURPOSE:
;
;	Analyze best-fit dust models to PKS 1413+135
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
;       Written by K. Willett                Nov 09
;-

pro clumpy_analysis, obj, sphere=sphere, error = error

tag, obj, dirtag
restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/torus_sav/'+obj+'_torus_grid.sav'
restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'

; Parameter space 

sig = [15, 30, 45, 60, 75]
Y = [5, 10, 30, 100, 200]
N0 = [1, 5, 10, 15, 20]
q = [0.0, 1.0, 2.0, 3.0]
tauv = [10, 30, 60, 80, 100, 200, 300, 500]
inc = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

; Find best-fit models in parameter grid

minind = where(errorgrid eq min(errorgrid))
minvalues = errorgrid_parse(minind)

er = 100 * (errorgrid - min(errorgrid)) / min(errorgrid)

if not keyword_set(error) then error=20
error_ind = where(er le error)

nerr = n_elements(error_ind)

errorarr = fltarr(6, nerr)

for i = 0, nerr - 1 do begin
	
	ind = errorgrid_parse(error_ind[i])
	errorarr[*, i] = minvalues

endfor

medianfits = median(errorarr, dim=2)
print,''
print,'Median best fit in '+strtrim(error,2)+'% range (TORUS, '+strtrim(nerr,2)+' models): ', medianfits
print,'Best fit (TORUS):                     ', $
	minvalues[0], minvalues[1], minvalues[2], $
	minvalues[3], minvalues[4], minvalues[5]

!p.multi=[0,2,3]

cs = 3.0

histoplot, errorarr[0,*],  charsize = cs, xtitle='!7r!3', xr=[-5,95], min_value=1, bin=10., title=sed.obj
histoplot, errorarr[1,*],  charsize = cs, xtitle='Y', min_value=1, bin=2., title=sed.tag, xr=[0, 210]
histoplot, errorarr[2,*],  charsize = cs, xtitle='N!I0!N', xr=[0,25], min_value=1, bin=1.
histoplot, errorarr[3,*],  charsize = cs, xtitle='q', xr=[-0.5,3], /xstyle, bin=0.5, /half, min_value=1
histoplot, errorarr[4,*],  charsize = cs, xtitle='!7s!3!IV!N', min_value=1, bin=2., xr=[0, 510]
histoplot, errorarr[5,*],  charsize = cs, xtitle='i', xr=[-5, 95], min_value=1, bin=10., /half

!p.multi=[0,1,1]


print,''

end
