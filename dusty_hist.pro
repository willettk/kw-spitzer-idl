;+
; NAME:
;       
;	DUSTY_HIST
;
; PURPOSE:
;
;	Analyze best-fit dust models with histograms of parameter space
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
;	PAHFIT		- runs on models fit to PAHFIT-subtracted continuum
;
;	ERRPERCENT	- error percentage over which to plot models (eg, ERRPERCENT = 20 plots all models
;				with a chi^2 value within 20% of the best fit)
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	ALLTAU is now default - would need to run simulations for CSOs again to fully populate grid
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
;	Modified for new tau_V grid, added PAHFIT keyword - Feb 10
;-

pro dusty_hist, $
	obj, $
	y_min, q_min, tauv_min, tdust_min, $
	pahfit=pahfit, errpercent = errpercent, $
	noplot = noplot, silent = silent

; Restore data

tag, obj, dirtag
targets, obj, r, iras, d

if keyword_set(pahfit) then pahfit_string = 'pahfit' else pahfit_string = ''

restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+obj+'_sphere_'+pahfit_string+'grid.sav'

; Parameter space 

array_Y = [2, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 400, 500, 750, 1000]
array_q = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

; Find best-fit models in parameter grid

minind = where(errorgrid eq min(errorgrid[where(errorgrid gt 0.)]))
minvalues = errorgrid_parse(minind,/dusty,/alltau)

	Y_min      = minvalues[0]
	q_min      = minvalues[1]
	s_min_temp = minvalues[2] + 1

	if s_min_temp gt 90 then begin
		s_min = s_min_temp - 90
		taudir = 'hightau'
		taustring = '_hightau'
		hightau = 1
	endif else begin
		s_min = s_min_temp
		taudir = 'regtau'
		taustring = ''
		hightau = 0
	endelse

outfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/out/'+taudir+'/kw_sphere_Y'$
	+string(Y_min,format='(i04)')+'_q'+string(q_min,format='(f3.1)')+taustring+'.out'
readcol, outfile, model, tau0, $
	f1, r1, r1_rc, theta1, tdust_arr, err, $
	/silent
tdust_min = tdust_arr[s_min - 1]

sppfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/spp/'+taudir+'/kw_sphere_Y'+$
	string(Y_min,format='(i04)')+'_q'+string(q_min,format='(f3.1)')+taustring+'.spp'
readcol, sppfile, model, tau0, format='i,f', /silent
tauv_min = tau0[s_min - 1]

; Create error grid around best fit

er = 100 * (errorgrid - min(errorgrid[where(errorgrid gt 0.)])) / min(errorgrid[where(errorgrid gt 0.)])

if n_elements(errpercent) eq 0 then errpercent = 20.
errarr = where(er le errpercent and er gt 0.)

n20 = n_elements(errarr)

arr20 = fltarr(3, n20)

for i = 0, n20 - 1 do begin
	
	ind = errorgrid_parse(errarr[i],/dusty,/alltau)
	if ind[2] gt 89 then s_ind = ind[2] - 90 else s_ind = ind[2]
	arr20[*, i] = [ind[0], ind[1], tau0[s_ind]]

endfor

if n20 gt 1 then medianfits = median(arr20, dim=2) else medianfits = arr20

if ~keyword_set(silent) then begin
	print,''
	print,'Median best fit in 20% range (SPHERE, '+strtrim(n20,2)+' models): ', medianfits
	print,'Best fit (TORUS):                                   ', Y_min, q_min, tauv_min
endif

; Plot histograms for distribution of Y, q, and tau_V

!p.multi=[0,1,3]

if ~keyword_set(noplot) then begin
	if n20 gt 1 and ~keyword_set(noplot) then begin
	
		plothist, arr20[0,*], /half, charsize = 2.0, xtitle='Y', xr=[0,1100], title=iras+' ('+obj+')'
		if min(arr20[1,*]) ne max(arr20[1,*]) then $
			plothist, arr20[1,*], /half, charsize = 2.0, xtitle='q', xr=[0,3], bin=0.5 $
			else $
				plot, [min(arr20[1,*])], [n20], xr=[0,3], yr=[0,n20+5], $
					psym=symcat(16), charsize=2.0, xtitle='q'
		plothist, arr20[2,*], /half, charsize = 2.0, xtitle='!7s!3!IV!N', xr=[0,500]
	
	endif else print,'No entries found within'+string(errpercent,format='(i3)')+'% of best fit for '+iras+' ('+obj+')'
endif

!p.multi=[0,1,1]

if ~keyword_set(silent) then print,''

end
