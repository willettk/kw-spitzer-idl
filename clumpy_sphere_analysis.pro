;+
; NAME:
;       
;	CLUMPY_SPHERE_ANALYSIS
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
;       Written by K. Willett                Jan 10
;-

pro clumpy_sphere_analysis, obj, sphere=sphere, error = error

tag, obj, dirtag
restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/sphere_sav/'+obj+'_sphere_grid.sav'
restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'

; Parameter space 

; Y = [5, 10, 30]
; N0 = [1, 2, 4, 6, 8, 10, 12, 15, 20]
; q = [0.0, 1.0, 2.0, 3.0]
; tauv = [10, 30, 60, 80, 100, 200, 300, 500]

; Find best-fit models in parameter grid

minind = where(errorgrid_sphere eq min(errorgrid_sphere))
minvalues = errorgrid_parse(minind, /sphere)

er = 100 * (errorgrid_sphere - min(errorgrid_sphere)) / min(errorgrid_sphere)

if not keyword_set(error) then error=20
error_ind = where(er le error)

nerr = n_elements(error_ind)

errorarr = fltarr(4, nerr)

for i = 0, nerr - 1 do begin
	
	ind = errorgrid_parse(error_ind[i], /sphere)
	errorarr[*, i] = minvalues

endfor

if nerr gt 1 then medianfits = median(errorarr, dim=2) else medianfits = errorarr

print,''
print,'Median best fit in '+strtrim(error,2)+'% range (SPHERE, '+strtrim(nerr,2)+' models): ', medianfits
print,'Best fit (SPHERE):                     ', $
	minvalues[0], minvalues[1], minvalues[2], minvalues[3]

!p.multi=[0,2,2]

cs = 1.5

histoplot, errorarr[0,*],  charsize = cs, xtitle='Y', min_value=1, bin=2., title=sed.tag, xr=[0, 210]
histoplot, errorarr[1,*],  charsize = cs, xtitle='N!I0!N', xr=[0,25], min_value=1, bin=1.
histoplot, errorarr[2,*],  charsize = cs, xtitle='q', xr=[-0.5,3], /xstyle, bin=0.5, /half, min_value=1
histoplot, errorarr[3,*],  charsize = cs, xtitle='!7s!3!IV!N', min_value=1, bin=2., xr=[0, 510]

!p.multi=[0,1,1]

print,''

end
