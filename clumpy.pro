;+
; NAME:
;       
;	CLUMPY
;
; PURPOSE:
;
;	Compare IRS spectra to a library of models to constrain the dust geometry
;
; INPUTS:
;
;	sigma 	- angular thickness of the torus as viewed from equatorial plane [degrees]
;
;	Y 	- ratio of outer and inner radii of the dust torus
;
;	N0 	- average number of dust clumps seen along the line of sight
;
;	q 	- index of power law describing radial distribution of dust clumps (goes as r^(-q))
;
;	tau_V	- optical depth of a dust cloud
;
;	inc 	- inclination angle of the dust torus wrt the observer [degrees]
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
;	sigma = 75 degrees
;	Y = 10
;	N0 = 15
;	q = 0.0
;	tau_V = 10
;	inc = 20 degrees
;
;	IDL> clumpy, 75, 10, 15, 0, 10, 20
;
; NOTES:
;
;	Based on models available at http://newton.pa.uky.edu/~clumpyweb
;
;	Note: models must be downloaded to local machine for CLUMPY.pro to run properly. Might consider adding a 
;		spawn command for wget to generate models if they don't exist?
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
; 	Modified for proper rebinning of data to model spectral resolution - KW, Nov 09
;	ERRORFIT_LFL, ERRORFIT_FL stored as separate functions - Nov 09
;-

;;;;;;;;;;;;;;;;;;;;

; Parameter space

; Y - 5, 10, 30, 100, 300
; N0 - 1, 5, 10, 15, 20
; sig - 15, 30, 45, 60, 75
; q - 0.0, 1.0, 2.0, 3.0
; tv - 10, 30, 60, 80, 100, 200, 500, 500
; i - 0, 10, 20, 30, 40, 50, 60, 70, 80, 90

;;;;;;;;;;;;;;;;;;;;

pro clumpy, sig, Y, N, q, tauv, inc, $
	errormerit, $
	wave, flux, modelwave, modelflux, $
	xr = xr, yr=yr, $
	noplot = noplot, quiet=quiet, stop=stop, obj=obj, bestfit = bestfit

common func_lfl, fmodel_lfl, fdata_lfl, ferror_lfl, function_lambda
common func_fl,  fmodel_fl,  fdata_fl,  ferror_fl

tag, obj, dirtag

; Best-fit keyword (if batch mode has already been run)

if keyword_set(bestfit) then begin

	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/torus_sav/'+obj+'_torus_grid.sav'

	array_sig = [15, 30, 45, 60, 75]
	array_Y = [5, 10, 30, 100, 200]
	array_N0 = [1, 5, 10, 15, 20]
	array_q = [0.0, 1.0, 2.0, 3.0]
	array_tauv = [10, 30, 60, 80, 100, 200, 300, 500]
	array_inc = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

	minind = where(errorgrid eq min(errorgrid))
	minvalues = errorgrid_parse(minind)
	
	sig = minvalues[0]
	Y=minvalues[1]
	N=minvalues[2]
	q=minvalues[3]
	tauv=minvalues[4]
	inc=minvalues[5]

endif

; Read in data from downloaded CLUMPY models

readcol, '~/Astronomy/Research/Spitzer/CSO/CLUMPY/sig'+$
	string(sig,format='(i2)')+'/SED+AGN-AA00-TORUSG-sig'+$
	string(sig,format='(i2)')+'-Y'+$
	string(Y, format='(i03)')+'-N'+$
	string(n,format='(i02)')+'-q'+$
	string(q,format='(f3.1)')+'-tv'+$
	string(tauv,format='(f05.1)')+'.txt', $
	lambda, i0, i10, i20, i30, i40, i50, i60, i70, i80, i90, $
	/silent

case inc of
	0:  begin
		modelinc = i0	
		fitcolor=fsc_color("White")
	    end
	10: begin
		modelinc = i10	
		fitcolor=fsc_color("Red")
	    end
	20: begin
		modelinc = i20	
		fitcolor=fsc_color("Orange")
	    end
	30: begin
		modelinc = i30	
		fitcolor=fsc_color("Yellow")
	    end
	40: begin
		modelinc = i40	
		fitcolor=fsc_color("Green")
	    end
	50: begin
		modelinc = i50	
		fitcolor=fsc_color("Cyan")
	    end
	60: begin
		modelinc = i60	
		fitcolor=fsc_color("Blue")
	    end
	70: begin
		modelinc = i70	
		fitcolor=fsc_color("Violet")
	    end
	80: begin
		modelinc = i80	
		fitcolor=fsc_color("Pink")
	    end
	90: begin
		modelinc = i90	
		fitcolor=fsc_color("Brown")
	    end
endcase

; Read in data

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'

flux = sed.flux_lr
wave = sed.wave_lr
err = sed.err_lr

noneg = where(flux gt 0)

flux = flux[noneg]
wave = wave[noneg]
err = err[noneg]

c = 299794.258

; Rebin the DATA to the same wavelengths as the model

minwave = closeto(lambda, wave[0])
maxwave = closeto(lambda, wave[n_elements(wave) - 1])

truncated_model_lambda = lambda[minwave:maxwave]
rebinned_flux = fltarr(n_elements(truncated_model_lambda))
rebinned_err  = fltarr(n_elements(truncated_model_lambda))

;for i=0, n_elements(wave) - 1 do begin
for i=0, n_elements(truncated_model_lambda) - 1 do begin

	modelwaveind = closeto(wave, truncated_model_lambda[i])			
	rebinned_flux[i] = flux[modelwaveind]
	rebinned_err[i] = err[modelwaveind]

	;modelwaveind = closeto(lambda, wave[i])
	;fmodel_lfl[i] = modelinc[modelwaveind]
	;fmodel_fl[i]  = modelinc[modelwaveind] / wave[i]

endfor

; Convert rebinned flux, errors in data from Janskys to W/m^2

flambda     = rebinned_flux * c / (truncated_model_lambda * 1d-6)^2 * 1d-26
flambda_err = rebinned_err  * c / (truncated_model_lambda * 1d-6)^2 * 1d-26

fdata_lfl = 1d-6 * flambda / 1d-15
fdata_fl = flambda / 1d-10

function_lambda = truncated_model_lambda
ferror_lfl = function_lambda * 1d-6 * flambda_err / 1d-15
ferror_fl = flambda_err / 1d-10

fmodel_lfl = modelinc[minwave:maxwave]
fmodel_fl = modelinc[minwave:maxwave] / lambda[minwave:maxwave]

; max(1./n * sqrt( total( ( (p[0] * fmodel_lfl - function_lambda * fdata_lfl)/ferror_lfl)^2) ) )

; Need:

; p[0] - 		fitted scale factor corresponding to AGN bolometric luminosity
; fmodel - 		model flux
; function_lambda - 	wavelength grid
; fdata - 		data
; ferror - 		measured error

result_lfl = amoeba(1.0d-5, scale=1d1, p0=[10.], function_value=fval_lfl, function_name='errorfit_lfl', ncalls=ncalls_lfl)
result_fl  = amoeba(1.0d-5, scale=1d1, p0=[10.], function_value=fval_fl , function_name='errorfit_fl' , ncalls=ncalls_fl)

if not keyword_set(quiet) then begin

	print,''
	print, 'sigma = ', sig,' deg'
	print, 'Y     = ', Y
	print, 'N0    = ',n
	print, 'q     = ',q
	print, 'tau_V = ', tauv
	print, 'inc   = ', inc,' deg'
	print,''


	print, 'lambda F_lambda'
	print,''
	print, 'Result:         ',result_lfl
	print, 'Function value: ', fval_lfl[0]
	print, 'ncalls:         ', ncalls_lfl
	print,''
;	print, 'F_lambda'
;	print,''
;	print, 'Result:         ',result_fl
;	print, 'Function value: ', fval_fl[0]
;	print, 'ncalls:         ', ncalls_fl
endif

errormerit = fval_lfl[0]

; Plotting

if not keyword_set(noplot) then begin

	!p.multi=[0,2,1]
	
	; LAMBDA * F_LAMBDA
	
	if not keyword_set(xr) then xr = [3,40]
	if not keyword_set(yr) then yr = [min(wave * 1d-6 * flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-15),1d0]

	plot, lambda, i0, $
		/xlog, $
		/ylog, $
		title='N!I0!N='+string(n,format='(i2)')+$
			', q='+string(q,format='(f3.1)')+$
			', !7s!3!IV!N='+string(tauv,format='(f5.1)')+ $
			', Y='+string(y,format='(i3)')+ $
			', !7r!3='+string(sig,format='(i2)')+'!Eo!N'+ $
			', i='+string(inc,format='(i2)')+'!Eo!N', $
		xtitle='Rest wavelength [!7l!3m]', $
		ytitle='!7k!3F!I!7k!3!N', $
		charsize=1.5, $
		thick=2, $
		xr = xr, $
		yr = yr, $
		/xstyle, $
		/ystyle
	
	oplot, lambda, i10, thick=2, color=fsc_color("Red")
	oplot, lambda, i20, thick=2, color=fsc_color("Orange")
	oplot, lambda, i30, thick=2, color=fsc_color("Yellow")
	oplot, lambda, i40, thick=2, color=fsc_color("Green")
	oplot, lambda, i50, thick=2, color=fsc_color("Cyan")
	oplot, lambda, i60, thick=2, color=fsc_color("Blue")
	oplot, lambda, i70, thick=2, color=fsc_color("Violet")
	oplot, lambda, i80, thick=2, color=fsc_color("Pink")
	oplot, lambda, i90, thick=2, color=fsc_color("Brown")
	
	;oplot, truncated_model_lambda, truncated_model_lambda * 1d-6 * flambda / 1d-15, psym=10
	oplot, wave, wave * 1d-6 * flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-15, psym=10

	; Overplot scaled model
	
	oplot, truncated_model_lambda, result_lfl[0] * fmodel_lfl, thick=2, color=fitcolor

	modelwave = truncated_model_lambda
	modelflux = result_lfl[0] * fmodel_lfl
	
	; F_LAMBDA
	
	plot, lambda, i0/lambda, $
		/xlog, $
		/ylog, $
		title=sed.obj, $
		xtitle='Rest wavelength [!7l!3m]', $
		ytitle='F!I!7k!3!N', $
		charsize=1.5, $
		thick=2, $
		yr=[min(flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-10),1d0], /ystyle, $
		xr=[3,40], /xstyle
	
	oplot, lambda, i10/lambda, thick=2, color=fsc_color("Red")
	oplot, lambda, i20/lambda, thick=2, color=fsc_color("Orange")
	oplot, lambda, i30/lambda, thick=2, color=fsc_color("Yellow")
	oplot, lambda, i40/lambda, thick=2, color=fsc_color("Green")
	oplot, lambda, i50/lambda, thick=2, color=fsc_color("Cyan")
	oplot, lambda, i60/lambda, thick=2, color=fsc_color("Blue")
	oplot, lambda, i70/lambda, thick=2, color=fsc_color("Violet")
	oplot, lambda, i80/lambda, thick=2, color=fsc_color("Pink")
	oplot, lambda, i90/lambda, thick=2, color=fsc_color("Brown")
	
	;oplot, truncated_model_lambda, flambda / 1d-10, psym=10
	oplot, wave, flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-10, psym=10

	; Overplot scaled model
	
	oplot, truncated_model_lambda, result_fl[0] * fmodel_fl, thick=2, color=fitcolor

	!p.multi=[0,1,1]

endif

; Save best-fit data to IDL file

if keyword_set(bestfit) then save, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/torus_sav/'+obj+'_torus_bestfit.sav', $
	wave, flux, modelwave, modelflux

if keyword_set(stop) then stop

end
