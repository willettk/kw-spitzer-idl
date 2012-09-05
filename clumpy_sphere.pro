;+
; NAME:
;       
;	CLUMPY_SPHERE
;
; PURPOSE:
;
;	Compare IRS spectra to a library of models to constrain the dust geometry in a sphere
;
; INPUTS:
;
;	Y 	- ratio of outer and inner radii of the dust torus
;
;	N0 	- average number of dust clumps seen along the line of sight
;
;	q 	- index of power law describing radial distribution of dust clumps (goes as r^(-q))
;
;	tau_V	- optical depth of a dust cloud
;
; OUTPUTS:
;
;	ERRORMERIT - 		parameter which assesses the goodness-of-fit of a particular model (try to minimize)
;
; KEYWORDS:
;
;	NOPLOT - 	set to suppress plots of models and data
;
;	QUIET - 	set to suppress printing results to screen
;
;	STOP - 		set to stop at end of program
;
; EXAMPLE:
;
;	Y = 10
;	N0 = 15
;	q = 0.0
;	tau_V = 10
;
;	IDL> clumpy_sphere, 75, 10, 15, 0, 10, 20
;
; NOTES:
;
;	Based on models available at http://newton.pa.uky.edu/~clumpyweb
;
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
; 	Generalized for all Spitzer galaxies - Jan 10
;-

;;;;;;;;;;;;;;;;;;;;

; Parameter space

; Y - 5, 10, 30
; N0 - 1, 2, 4, 6, 8, 10, 12, 15, 20
; q - 0.0, 1.0, 2.0, 3.0
; tv - 10, 30, 60, 80, 100, 200, 300, 500

;;;;;;;;;;;;;;;;;;;;

pro clumpy_sphere, Y, N, q, tauv, errormerit, obj=obj, noplot = noplot, quiet=quiet, stop=stop

common func_lfl, fmodel_lfl, fdata_lfl, ferror_lfl, function_lambda
common func_fl,  fmodel_fl,  fdata_fl,  ferror_fl

tag, obj, dirtag

; Best-fit keyword (if batch mode has already been run)

if keyword_set(bestfit) then begin

	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/sphere_sav/'+obj+'_sphere_grid.sav'

	array_Y = [5, 10, 30]
	array_N0 = [1, 2, 4, 6, 8, 10, 12, 15, 20]
	array_q = [0.0, 1.0, 2.0, 3.0]
	array_tauv = [10, 30, 60, 80, 100, 200, 300, 500]

	minind = where(errorgrid eq min(errorgrid))
	minvalues = errorgrid_parse(minind,/sphere)
	
	Y    = minvalues[0]
	N    = minvalues[1]
	q    = minvalues[2]
	tauv = minvalues[3]

endif

; Read in data from downloaded CLUMPY models

readcol, '~/Astronomy/Research/Spitzer/CSO/CLUMPY/sphere/SED+AGN-AA00-SPHERE-Y'+$
	string(Y, format='(i03)')+'-N'+$
	string(n,format='(i02)')+'-q'+$
	string(q,format='(f3.1)')+'-tv'+$
	string(tauv,format='(f05.1)')+'.txt', $
	lambda, sphereflux, $
	/silent

; Read in data for object

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'

flux = sed.flux_lr
wave = sed.wave_lr
err = sed.err_lr

c = 299794.258

; Rebin the DATA to the same wavelengths as the model

minwave = closeto(lambda, wave[0])
maxwave = closeto(lambda, wave[n_elements(wave) - 1])

truncated_model_lambda = lambda[minwave:maxwave]
rebinned_flux = fltarr(n_elements(truncated_model_lambda))
rebinned_err  = fltarr(n_elements(truncated_model_lambda))

for i=0, n_elements(truncated_model_lambda) - 1 do begin

	modelwaveind = closeto(wave, truncated_model_lambda[i])			
	rebinned_flux[i] = flux[modelwaveind]
	rebinned_err[i] = err[modelwaveind]

endfor

; Convert rebinned flux, errors from Janskys to W/m^2

flambda     = rebinned_flux * c / (truncated_model_lambda * 1d-6)^2 * 1d-26
flambda_err = rebinned_err  * c / (truncated_model_lambda * 1d-6)^2 * 1d-26

fdata_lfl = 1d-6 * flambda / 1d-15
fdata_fl = flambda / 1d-10

function_lambda = truncated_model_lambda
ferror_lfl = function_lambda * 1d-6 * flambda_err / 1d-15
ferror_fl = flambda_err / 1d-10

fmodel_lfl = sphereflux[minwave:maxwave]
fmodel_fl = sphereflux[minwave:maxwave] / lambda[minwave:maxwave]

result_lfl = amoeba(1.0d-5, scale=1d1, p0=[10.], function_value=fval_lfl, function_name='errorfit_lfl', ncalls=ncalls_lfl)
result_fl  = amoeba(1.0d-5, scale=1d1, p0=[10.], function_value=fval_fl , function_name='errorfit_fl' , ncalls=ncalls_fl)

if not keyword_set(quiet) then begin

	print,''
	print, 'Y = ', Y
	print, 'N0 = ',n
	print, 'q = ',q
	print, 'tau_V = ', tauv
	print,''


	print, 'lambda F_lambda'
	print,''
	print, 'Result: ',result_lfl
	print, 'Function value: ', fval_lfl[0]
	print, 'ncalls: ', ncalls_lfl
	print,''
;	print, 'F_lambda'
;	print,''
;	print, 'Result: ',result_fl
;	print, 'Function value: ', fval_fl[0]
;	print, 'ncalls: ', ncalls_fl
endif

errormerit = fval_lfl[0]

; Plotting

if not keyword_set(noplot) then begin

	!p.multi=[0,2,1]
	
	; LAMBDA * F_LAMBDA
	
	plot, lambda, sphereflux, $
		/xlog, $
		/ylog, $
		title='N!I0!N='+string(n,format='(i2)')+$
			', q='+string(q,format='(f3.1)')+$
			', !7s!3!IV!N='+string(tauv,format='(f5.1)')+ $
			', Y='+string(y,format='(i3)'), $
		xtitle='Rest wavelength [!7l!3m]', $
		ytitle='!7k!3F!I!7k!3!N', $
		charsize=1.5, $
		thick=2, $
		yr=[1d-3,1d0], /ystyle, $
		xr=[3,40], /xstyle
	
	;oplot, truncated_model_lambda, truncated_model_lambda * 1d-6 * flambda / 1d-15, psym=10
	oplot, wave, wave * 1d-6 * flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-15, psym=10

	; Overplot scaled model
	
	oplot, truncated_model_lambda, result_lfl[0] * fmodel_lfl, thick=2, color=fsc_color("Red")
	
	; F_LAMBDA
	
	plot, lambda, sphereflux/lambda, $
		/xlog, $
		/ylog, $
		xtitle='Rest wavelength [!7l!3m]', $
		ytitle='F!I!7k!3!N', $
		charsize=1.5, $
		thick=2, $
		yr=[1d-3,1d0], /ystyle, $
		xr=[3,40], /xstyle
	

	;oplot, truncated_model_lambda, flambda / 1d-10, psym=10
	oplot, wave, flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-10, psym=10

	; Overplot scaled model
	
	oplot, truncated_model_lambda, result_fl[0] * fmodel_fl, thick=2, color=fsc_color("Red")

	!p.multi=[0,1,1]

endif

if keyword_set(stop) then stop

end
