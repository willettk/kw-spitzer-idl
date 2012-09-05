;+
; NAME:
;       
;	DUSTY_COMPARE_BATCH
;
; PURPOSE:
;
;	Run DUSTY on models to compare to data
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
;	IDL> dusty_compare_batch, 'cso005'
;
; NOTES:
;
;	PAHFIT must be compiled before running the program with the PAHFIT keyword
;
; REVISION HISTORY
;       Written by K. Willett                Nov 09
;-

pro dusty_compare_batch, obj, pahfit=pahfit, stop = stop

tag, obj, dirtag

;;;;;;;;;;;;;;;;;;;;

; Parameter space

Y = [2, 5, 10, 15, 20, 30, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 650, 700, 750, 800, 850, 900, 950, 1000]
q = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

;;;;;;;;;;;;;;;;;;;;

nregtau = 90
nhightau = 90

errorgrid = fltarr(n_elements(Y), n_elements(q), nregtau + nhightau)

;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; REGTAU ;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

timestart = systime(1)

for j = 0, n_elements(Y) - 1 do begin
	for l = 0, n_elements(q) - 1 do begin

		if (Y[j] eq 50 and q[l] eq 1.0) then nregtau = 87		; DUSTY didn't generate models 88-90 for these parameters
		if (Y[j] eq 75 and q[l] eq 1.0) then nregtau = 88		; DUSTY didn't generate models 89-90 for these parameters

		for m = 1, nregtau do begin

			dusty_compare, Y[j], q[l], m, $
				errormerit, obj = obj, $
				/quiet, $
				/noplot, $
				pahfit = pahfit
			
			errorgrid[j,l,m-1] = errormerit
			
		endfor
	endfor

	print, 'Y = '+string(Y[j],format='(i4)')+' completed for regtau'

endfor

timeend = systime(1)

print,''
print, 'Time elapsed for regtau: ',fix((timeend - timestart)/6.)/10.,' minutes'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; HIGHTAU ;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

timestart = systime(1)

for j = 0, n_elements(Y) - 1 do begin
	for l = 0, n_elements(q) - 1 do begin

		if (Y[j] eq 50 and q[l] eq 1.0) then nhightau = 71		; DUSTY didn't generate models for these parameters
		if (Y[j] eq 75 and q[l] eq 1.0) then nhightau = 58		; DUSTY didn't generate models for these parameters

		for m = 1, nhightau do begin

			dusty_compare, Y[j], q[l], m, $
				errormerit, obj = obj, $
				/quiet, $
				/hightau, $
				/noplot, $
				pahfit = pahfit
			
			errorgrid[j,l,m-1 + nregtau] = errormerit
			
		endfor
	endfor

	print, 'Y = '+string(Y[j],format='(i4)')+' completed for hightau'

endfor

timeend = systime(1)

print,''
print, 'Time elapsed for hightau: ',fix((timeend - timestart)/6.)/10.,' minutes'

; Find best-fit model

minerror = where(errorgrid eq min(errorgrid[where(errorgrid gt 0.)]))
if n_elements(minerror) gt 1 then print, 'Number of minima in errorgrid: ',n_elements(minerror)
minerrorval = errorgrid_parse(minerror[0],/dusty, /alltau)

	Y_min      = minerrorval[0]
	q_min      = minerrorval[1]
	s_min_temp = minerrorval[2] + 1

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

; Load .spp file to associate the model number with an optical depth

sppfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/spp/'+taudir+'/kw_sphere_Y'+string(Y_min,format='(i04)')+'_q'+string(q_min,format='(f3.1)')+taustring+'.spp'

readcol, sppfile, model, tau0, $
	format='i,f', $
;	tau0, psi, fv, fk, f12, c21, c31, c43, b8, b14, b9, b11, r9, $
	/silent

tauv = tau0[s_min - 1]
tauv_min=tauv

; Load .out file to associate model with dust temperature at the outer edge

outfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/out/'+taudir+'/kw_sphere_Y'+string(Y_min,format='(i04)')+'_q'+string(q_min,format='(f3.1)')+taustring+'.out'

readcol, outfile, model_arr, tauoutarr, $
	f1, r1, r1_rc, theta1, tdust_arr, $
	/silent

tdust = tdust_arr[s_min - 1]

; Print best fit to screen

print,''
print,'Torus best fit parameters: '
print,''
print, 'Y     = ', Y_min
print, 'q     = ', q_min
print, 'tauv  = ', tauv
print, 'Tdust  = ', tdust
print,''
print,taudir
print,''

; Save results to IDL files

if keyword_set(pahfit) then pahfit_string = 'pahfit' else pahfit_string = ''
savefile = '~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+obj+'_sphere_'+pahfit_string+'grid.sav'

save, filename=savefile, $
	errorgrid, Y_min, q_min, tauv_min, tdust, s_min, taudir, taustring

; Plot best-fit on screen

dusty_bestfit, obj, pahfit = pahfit

if keyword_set(stop) then stop

end
