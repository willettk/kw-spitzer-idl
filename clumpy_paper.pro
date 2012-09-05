;+
; NAME:
;       
;	CLUMPY_PAPER
;
; PURPOSE:
;
;	Create figure w/spectra and CLUMPY model fit for the CSO paper
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
;	IDL> clumpy, 75, 10, 15, 0, 10, 20
;
; NOTES:
;
;
; REVISION HISTORY
;       Written by K. Willett                Dec 09
;-

;;;;;;;;;;;;;;;;;;;;

pro clumpy_paper, ps=ps, stop=stop

; Read in data

c = 299794.258

restore, '~/Astronomy/Research/Spitzer/CSO/CLUMPY/torus_sav/cso007_torus_bestfit.sav'

wave_1946 = wave
flux_1946 = flux
modelwave_1946 = modelwave
modelflux_1946 = modelflux

restore, '~/Astronomy/Research/Spitzer/CSO/CLUMPY/torus_sav/cso004_torus_bestfit.sav'

wave_oq208 = wave
flux_oq208 = flux
modelwave_oq208 = modelwave
modelflux_oq208 = modelflux

restore, '~/Astronomy/Research/Spitzer/CSO/CLUMPY/torus_sav/cso005_torus_bestfit.sav'

wave_pks1413 = wave
flux_pks1413 = flux
modelwave_pks1413 = modelwave
modelflux_pks1413 = modelflux

; Plotting

if not keyword_set(noplot) then begin

	!p.multi=[0,1,1]

	if keyword_set(ps) then begin
	;	!p.font=0
		set_plot,'ps'
		device, filename = '~/Astronomy/Research/Spitzer/cso/papers/clumpy_paper.ps', /color, /portrait
		cthick = 5
		lthick = 5
		cs = 1.7
	endif else begin
		cthick = 1
		lthick = 1
		cs = 1
	endelse

	; LAMBDA * F_LAMBDA
	
	xr = [4,40]
	yr = [2d-4,2d-1]

	plot, indgen(10), $
		/nodata, $
		/xlog, $
		/ylog, $
		xtitle='Rest wavelength [!7l!3m]', $
		ytitle='!7k!3F!I!7k!3!N [10!E-15!N Wm!E-2!N]', $
		charsize=cs, $
		charthick=cthick, $
		yminor=5, $
		yticks=2, $
		ytickname=['10!E-1!N','10!E-2!N','10!E-3!N'], $
		ytickv=[1d-1,1d-2,1d-3], $
		xticks=3, $
		xtickname=['5','10','20','30'], $
		xtickv=[5,10,20,30], $
		thick=lthick, $
		xthick=lthick, $
		ythick=lthick, $
		xr = xr, $
		yr = yr, $
		/xstyle, $
		/ystyle
	
	oplot, wave_1946, wave_1946 * 1d-6 * flux_1946 * c / (wave_1946 * 1d-6)^2 * 1d-26 / 1d-15, $
		psym=10, thick = lthick
	oplot, wave_oq208, wave_oq208 * 1d-6 * flux_oq208 * c / (wave_oq208 * 1d-6)^2 * 1d-26 / 1d-15, $
		psym=10, thick = lthick
	oplot, wave_pks1413, wave_pks1413 * 1d-6 * flux_pks1413 * c / (wave_pks1413 * 1d-6)^2 * 1d-26 / 1d-15, $
		psym=10, thick = lthick

	; Overplot scaled model
	
	oplot, modelwave_1946, modelflux_1946, $
		color=fsc_color("Red"), thick = lthick
	oplot, modelwave_oq208, modelflux_oq208, $
		color=fsc_color("Red"), thick = lthick
	oplot, modelwave_pks1413, modelflux_pks1413, $
		color=fsc_color("Red"), thick = lthick

	if keyword_set(ps) then begin
		device,/close
		set_plot,'x'
	endif


	!p.multi=[0,1,1]

endif

if keyword_set(stop) then stop

end
