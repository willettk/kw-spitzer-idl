
;+
; NAME:
;       
;	Create spectra with examples of strong absorption features for OHM Paper I
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
;       Written by K. Willett                Mar 2010
;-


restore,'~/Astronomy/Research/Spitzer/archived/data/structures/arch029.sav'

wave_sh = sed.wave_sh
flux_sh = sed.flux_sh

wave_lh = sed.wave_lh
flux_lh = sed.flux_lh

wave_lr = sed.wave_lr
flux_lr = sed.flux_lr

!p.font=-1

set_plot,'ps'
device,filename = '~/Astronomy/Research/Spitzer/papers/abs_example.ps'
thick=3
cthick=3

multiplot, [2,2], $
	mxtitle = 'Wavelength [!7l!3m]', $
	mytitle = 'Flux density [Jy]', $
	mxtitsize = 1.2, $
	mytitsize = 1.2, $
	mcharthick = cthick, $
	mxcharthick = cthick, $
	mycharthick = cthick

; H2O water ice

plot, wave_lr, flux_lr, $
	xstyle=9, $
	ystyle=9, $
	psym=10, $
	thick=thick, $
	xthick=thick, $
	ythick=thick, $
	xr=[5.1,7.4], $
	yr=[0.1,0.5]

axis, xaxis=1, xstyle=1, xtickformat='(f4.1)', charthick=3, xthick=thick
axis, yaxis=0, ystyle=1, ytickformat='(f4.1)', charthick=3

plots, [6.0,6.0]-0.05, [0.35,0.2], linestyle=2, thick=thick
plots, [6.8,6.8], [0.4,0.3], linestyle=2, thick=thick
xyouts, 5.65, 0.38, 'H!I2!NO ice', /data, charthick=3
xyouts, 6.7, 0.42, 'HAC', /data, charthick=3

multiplot

; Gas-phase C2H2, HCN, CO2

plot, wave_sh, flux_sh, $
	xstyle=9, $
	ystyle=9, $
	psym=10, $
	thick=thick, $
	xthick=thick, $
	ythick=thick, $
	xr=[13.45,14.25], $
	yr=[0.61,1.1]

axis, xaxis=1, xstyle=1, xtickformat='(f4.1)', charthick=3, xthick=thick
axis, yaxis=1, ystyle=1, ytickformat='(f4.1)', charthick=3, ythick=thick

plots, [13.7,13.7], [1.0,0.85], linestyle=2, thick=thick
plots, [14.04,14.04], [1.0,0.9], linestyle=2, thick=thick
xyouts, 13.67, 1.03, 'C!I2!NH!I2!N', /data, charthick=3
xyouts, 14.0, 1.03, 'HCN', /data, charthick=3

multiplot

; 34.6 um OH

plot, wave_lh, flux_lh, $
	psym = 10, $
	xtickformat='(f4.1)', $
	xstyle=1, $
	ystyle=9, $
	thick=thick, $
	xthick=thick, $
	ythick=thick, $
	charthick=3, $
	xr=[33.5,35.4], $
	yr=[30,58]
multiplot

axis, yaxis=0, ystyle=1, ytickformat='(f4.1)', charthick=3
plots, [34.6,34.6], [50,42], linestyle=2, thick=thick
xyouts, 34.52, 52, 'OH', /data, charthick=3

; Crystalline silicate

fname='arch029'
tag, fname, dirtag
saveset_file = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/silicate/'+fname+'_cryssil.sav'

saveset_test = file_test(saveset_file)
restore, saveset_file

posflux = where(flux gt 0)
wave = sed.wave_lr[posflux]
flux = sed.flux_lr[posflux]

; Spline fits

if wave[n_elements(wave)-1] lt 30. then redend = 28. else redend = 30.

wcont = wave([closetomed(wave,flux,5.6),closetomed(wave,flux,7.1),closetomed(wave,flux,redend)]) 
fcont = flux([closetomed(wave,flux,5.6),closetomed(wave,flux,7.1),closetomed(wave,flux,redend)]) 
yspl = spl_init(wcont, fcont)
splresult = spl_interp(wcont, fcont, yspl, wave)

tau = alog(flux / splresult)
	
wcont_amorph = fltarr(n_elements(amorph_sil))
tcont_amorph = fltarr(n_elements(amorph_sil))
pivotpoints_amorph = fltarr(n_elements(amorph_sil))

for i = 0, n_elements(amorph_sil) - 1 do begin
	wcont_amorph[i] = wave[closetomed(wave,tau,amorph_sil[i])]
	tcont_amorph[i] = tau[closetomed(wave,tau,amorph_sil[i])]
endfor

yspl_amorph = spl_init(wcont_amorph, tcont_amorph)
splresult_amorph = spl_interp(wcont_amorph, tcont_amorph, yspl_amorph, wave)
cryst_resid = tau - splresult_amorph


; Plot the results

plotind = where(cryst_resid lt 0.1)
plot, wave[plotind], cryst_resid[plotind], $
	xtickformat='(f4.1)', $
	xr = [12.1,20], $
	yr=[-0.3,0.3], $
	xstyle=1, $
	charthick=3, $
	thick=thick, $
	xthick=thick, $
	ythick=thick, $
	psym = 10

axis, yaxis=1, ystyle=1, ytickformat='(f4.1)', charthick=3
plots, [16.0,16.0], [0.2,0.0], linestyle=2, thick=thick
xyouts, 14.8,0.25, 'Crystalline sil.', /data, charthick=3

multiplot

multiplot,/reset

device,/close
set_plot,'x'

end
