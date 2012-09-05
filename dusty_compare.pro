;+
; NAME:
;       
;	DUSTY_COMPARE
;
; PURPOSE:
;
;	Compare IRS spectra to a library of models to constrain the dust geometry
;
; INPUTS:
;
;	Y 	- ratio of outer and inner radii of the dust torus
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
;	IDL> dusty_compare, 75, 10, 15, 0, 10, 20
;
; NOTES:
;
;	Based on models available at http://newton.pa.uky.edu/
;
;	Note: models must be downloaded to local machine for routine to run properly. Might consider adding a 
;		spawn command for wget to generate models if they don't exist?
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
; 	Modified for proper rebinning of data to model spectral resolution - KW, Nov 09
;	ERRORFIT_LFL, ERRORFIT_FL stored as separate functions - Nov 09
; 	Updated directory definitions - Feb 10
;	Removed BESTFIT keyword; created separate program as DUSTY_BESTFIT - Feb 10
;	Only run lambda f_lambda for faster batch mode
;-

;;;;;;;;;;;;;;;;;;;;

; Parameter space

;;;;;;;;;;;;;;;;;;;;

pro dusty_compare, Y, q, s, $
	errormerit, $
	noplot = noplot, quiet=quiet, stop=stop, obj=obj, $
	hightau = hightau, $
	pahfit=pahfit

common func_lfl, fmodel_lfl, fdata_lfl, ferror_lfl, function_lambda
;common func_fl,  fmodel_fl,  fdata_fl,  ferror_fl

tag, obj, dirtag

if keyword_set(pahfit) then begin
	pahfit_string = 'pahfit' 
endif else begin
	pahfit_string = ''
endelse

if keyword_set(hightau) then begin
	taudir = 'hightau' 
	taustring = '_hightau'
	alltau = 1
endif else begin 
	taudir = 'regtau'
	taustring = ''
	alltau = 0
endelse

; Translate model number to optical depth

sppfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/spp/'+taudir+'/kw_sphere_Y'+string(Y,format='(i04)')+'_q'+string(q,format='(f3.1)')+taustring+'.spp'

readcol, sppfile, model, tau0, $
	format='i,f', $
;	tau0, psi, fv, fk, f12, c21, c31, c43, b8, b14, b9, b11, r9, $
	/silent

tauv = tau0[s - 1]

outfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/out/'+taudir+'/kw_sphere_Y'+string(Y,format='(i04)')+'_q'+string(q,format='(f3.1)')+taustring+'.out'

readcol, outfile, model, tau0, $
	f1, r1, r1_rc, theta1, tdust_arr, err, $
	/silent

tdust = tdust_arr[s - 1]

; Read in data from previously-run DUSTY model

readcol, '~/Astronomy/Research/Spitzer/CSO/DUSTY/spectra/'+taudir+'/kw_sphere_Y'+$
	string(Y,format='(i04)')+'_q'+$
	string(q,format='(f3.1)')+taustring+'.s'+ $
	string(s,format='(i03)'),$
	lambda, fTot, xAtt, xDs, xDe, fInp, tauT, albedo, $   
	/silent

; Read in data

if keyword_set(pahfit) then begin
	
	dusty_ohm_pahfit, obj, wave, flux, err, /noplot

endif else begin
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'
	
	flux = sed.flux_lr
	wave = sed.wave_lr
	err = sed.err_lr

endelse

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

;fdata_fl = flambda / 1d-10
;ferror_fl = flambda_err / 1d-10
;fmodel_fl = fTot[minwave:maxwave] / lambda[minwave:maxwave]
;result_fl  = amoeba(1.0d-5, scale=1d1, p0=[10.], function_value=fval_fl , function_name='errorfit_fl' , ncalls=ncalls_fl)

if not keyword_set(quiet) then begin

	print,''
	print, 'Y     = ', Y
	print, 'q     = ',q
	print, 'model = ',s
	print, 'tauv = ',tauv
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

	!p.multi=[0,1,1]
;	!p.multi=[0,2,1]
	
	; LAMBDA * F_LAMBDA
	
	plot, lambda, ftot, $
		/nodata, $
		/xlog, $
		/ylog, $
		title='q='+string(q,format='(f3.1)')+$
			', !7s!3!IV!N='+string(tauv,format='(f5.1)')+ $
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
	
;	; F_LAMBDA
;	
;	plot, lambda, ftot/lambda, $
;		/nodata, $
;		/xlog, $
;		/ylog, $
;;		title=sed.obj, $
;		xtitle='Rest wavelength [!7l!3m]', $
;		ytitle='F!I!7k!3!N', $
;		charsize=1.5, $
;		thick=2, $
;		yr=[min(flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-10),1d0], /ystyle, $
;		xr=[3,40], /xstyle
;	
;	;oplot, truncated_model_lambda, flambda / 1d-10, psym=10
;	oplot, wave, flux * c / (wave * 1d-6)^2 * 1d-26 / 1d-10, psym=10
;
;	; Overplot scaled model
;	
;	oplot, truncated_model_lambda, result_fl[0] * fmodel_fl, thick=2, color=fsc_color("Red")

	!p.multi=[0,1,1]

endif

if keyword_set(stop) then stop

end
