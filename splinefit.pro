pro splinefit, wave, flux, splresult, pivotpoints, spl=spl
;+
; NAME:
;       
;	SPLINEFIT
;
; PURPOSE:
;
;	Fit an interpolated spline to a spectrum from a selected range of continuum fitting choices
;
; INPUTS:
;
;	WAVE - array of wavelengths for spectrum
;
;	FLUX - array of measured fluxes for spectrum
;
; OUTPUTS:
;
;	SPLRESULT - array with fitted spline over same range as WAVE
;
;	PIVOTPOINTS - indices of locations where spline is fit
;
; KEYWORDS:
;
;	SPL - 		integer giving type of spline fit for used for the continuum
;
;			1: continuum-dominated - power-law from 5.0-7.5 and 27.0-31.5 um, intermediate spline anchor at 14.0 um
;			2: PAH-dominated - power-law from 5.5-14.5 um, spline anchors at 14.5 and 27.0 um
;			3: absorption-dominated (DEFAULT) - spline fit anchored at 5.2, 5.6, 14.0, and 27.0 um
;			4: absorption-dominated w/no PAH flux at 8 um (adds anchor at 7.8 um)
;			5+: custom spline fits
;
;
; EXAMPLE:
;
;	IDL> splinefit, wave, flux, spl=1
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Jan 10
;-

case spl of
	1: begin
		if wave(0) gt 5.0 then begind = 5.5 else begind = 5.0
		m1 = alog10(flux(closetomed(wave,flux,7.5)) / flux(closetomed(wave,flux,begind))) $
			/ alog10(wave(closetomed(wave,flux,7.5)) / wave(closetomed(wave,flux,begind)))
	   	b1 = alog10(flux(closetomed(wave,flux,7.5))) - m1 * alog10(wave(closetomed(wave,flux,7.5)))
	   	powerfit1 = 10^(m1 * alog10(wave(0:closetomed(wave,flux,7.5))) + b1)
		m2 = alog10(flux(closetomed(wave,flux,31.5)) / flux(closetomed(wave,flux,27.0))) $
			/ alog10(wave(closetomed(wave,flux,31.5)) / wave(closetomed(wave,flux,27.0)))
	   	b2 = alog10(flux(closetomed(wave,flux,31.5))) - m2 * alog10(wave(closetomed(wave,flux,31.5)))
	   	powerfit2 = 10^(m2 * alog10(wave(closetomed(wave,flux,27.0):n_elements(wave)-1)) + b2)
		wcont = [wave(closetomed(wave,flux,7.5)),wave(closetomed(wave,flux,14.0)),wave(closetomed(wave,flux,27.0))]
		fcont = [flux(closetomed(wave,flux,7.5)),flux(closetomed(wave,flux,14.0)),flux(closetomed(wave,flux,27.0))]
		yspl = spl_init(wcont,fcont)
	   	spltemp = spl_interp(wcont, fcont, yspl, wave(closetomed(wave,flux,7.5):closetomed(wave,flux,27.0)))
		splresult = [powerfit1(0:n_elements(powerfit1)-2), spltemp, powerfit2(1:n_elements(powerfit2)-1)]
		pivotpoints = [closetomed(wave,flux,7.5),closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]
	   end
	2: begin 
		if wave(n_elements(wave)-1) lt 30.0 then endind = 28.0 else endind = 30.0
		wcont = wave([closetomed(wave,flux,14.5),closetomed(wave,flux,27.0),closetomed(wave,flux,endind)]) 
	   	fcont = flux([closetomed(wave,flux,14.5),closetomed(wave,flux,27.0),closetomed(wave,flux,endind)]) 
	   	yspl = spl_init(wcont, fcont)
	   	spltemp = spl_interp(wcont, fcont, yspl, wave(closetomed(wave,flux,14.5):n_elements(wave)-1))
	   	m = alog10(flux(closetomed(wave,flux,14.5)) / flux(closetomed(wave,flux,5.5))) $
			/ alog10(wave(closetomed(wave,flux,14.5)) / wave(closetomed(wave,flux,5.5)))
	   	b = alog10(flux(closetomed(wave,flux,14.5))) - m * alog10(wave(closetomed(wave,flux,14.5)))
	   	powerfit = 10^(m * alog10(wave(0:closetomed(wave,flux,14.5)-1)) + b)
	   	splresult = [powerfit,spltemp]
		pivotpoints = [closetomed(wave,flux,5.5),closetomed(wave,flux,14.5),closetomed(wave,flux,27.0),closetomed(wave,flux,endind)]
	   end
	3: begin
		wcont = wave([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]) 
		fcont = flux([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]) 
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]
	   end
	4: begin
		wcont = wave([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,7.8),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]) 
		fcont = flux([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,7.8),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]) 
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,7.8),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,27.0)]
	   end
	5: begin
		pp = [5,6,11.7,14,27]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
	6: begin
		pp = [5,6,11.7,14,27,29]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
	7: begin
		pp = [6,7,14,27,28,29]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
	8: begin
		pp = [5.5,7,14,27,29]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
	9: begin
		pp = [5,7,14,27,28,29]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
	10: begin
		pp = [5,6.7,13.5,26.5,31]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
	11: begin
		pp = [5.4,6.7,13.9,26.5,31]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
endcase


end
