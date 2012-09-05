;+
; NAME:
;       
;	DUSTY_BESTFIT
;
; PURPOSE:
;
;	Compare IRS spectra to a library of models to constrain the dust geometry
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
;	QUIET - 	set to suppress printing results to screen
;
;	STOP - 		set to stop at end of program
;
; EXAMPLE:
;
;	IDL> dusty_bestfit
;
; NOTES:
;
;
; REVISION HISTORY
;       Written by K. Willett                Feb 10
;-

pro dusty_bestfit, obj, $
	quiet=quiet, stop=stop,  $
	pahfit=pahfit

common func_lfl, fmodel_lfl, fdata_lfl, ferror_lfl, function_lambda

tag, obj, dirtag

if keyword_set(pahfit) then begin
	pahfit_string = 'pahfit' 
	dusty_ohm_pahfit, obj, wave, flux, err, /noplot
endif else begin
	pahfit_string = ''
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'
	flux = sed.flux_lr
	wave = sed.wave_lr
	err = sed.err_lr
endelse

dustfile='~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+obj+'_sphere_'+pahfit_string+'grid.sav'

restore, filename=dustfile

y = y_min
q = q_min
s = s_min

; Read in data from previously-run DUSTY model

readcol, '~/Astronomy/Research/Spitzer/CSO/DUSTY/spectra/'+taudir+'/kw_sphere_Y'+$
	string(Y,format='(i04)')+'_q'+$
	string(q,format='(f3.1)')+taustring+'.s'+ $
	string(s,format='(i03)'),$
	lambda, fTot, xAtt, xDs, xDe, fInp, tauT, albedo, $   
	/silent

; Read in data

noneg = where(flux gt 0)

flux = flux[noneg]
wave = wave[noneg]
err = err[noneg]

c = 299794.258

; Rebin the data to the same wavelengths as the model

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

; Convert rebinned flux, errors in data from Janskys to W/m^2

flambda     = rebinned_flux * c / (truncated_model_lambda * 1d-6)^2 * 1d-26
flambda_err = rebinned_err  * c / (truncated_model_lambda * 1d-6)^2 * 1d-26

fdata_lfl = 1d-6 * flambda / 1d-15

function_lambda = truncated_model_lambda
ferror_lfl = function_lambda * 1d-6 * flambda_err / 1d-15

fmodel_lfl = fTot[minwave:maxwave]

result_lfl = amoeba(1.0d-5, scale=1d1, p0=[10.], function_value=fval_lfl, function_name='errorfit_lfl', ncalls=ncalls_lfl)

if not keyword_set(quiet) then begin

	print,''
	print, 'Y     = ', Y
	print, 'q     = ',q
	print, 'model = ',s
	print, 'tauv = ',tauv_min
	print, 'T_dust = ',tdust
	print,''


	print, 'lambda F_lambda'
	print,''
	print, 'Result:         ',result_lfl
	print, 'Function value: ', fval_lfl[0]
	print, 'ncalls:         ', ncalls_lfl
	print,''
endif

errormerit = fval_lfl[0]

; Plotting

if not keyword_set(noplot) then begin

	; LAMBDA * F_LAMBDA
	
	plot, lambda, ftot, $
		/nodata, $
		/xlog, $
		/ylog, $
		title='q='+string(q,format='(f3.1)')+$
			', !7s!3!IV!N='+string(tauv_min,format='(f5.1)')+ $
			', Y='+string(y,format='(i4)'), $
		xtitle='Rest wavelength [!7l!3m]', $
		ytitle='!7k!3F!I!7k!3!N', $
		charsize=1.5, $
		thick=2, $
		;yr=[1d-4,1d-1], /ystyle, $
		yr=[min(wave * 1d-6 * flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-15)<min(result_lfl[0] * fmodel_lfl),1d0], /ystyle, $
		xr=[3,40], /xstyle
	
	;oplot, truncated_model_lambda, truncated_model_lambda * 1d-6 * flambda / 1d-15, psym=10
	oplot, wave, $
		wave * 1d-6 * flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-15, $
		psym=10

	; Overplot scaled model
	
	oplot, truncated_model_lambda, result_lfl[0] * fmodel_lfl, thick=2, color=fsc_color("Red")
	
endif

if keyword_set(stop) then stop

end
