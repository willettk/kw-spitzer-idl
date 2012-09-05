;+
; NAME:
;       
;	CRYSTALLINE
;
; PURPOSE:
;
;	Measure crystalline silicate components from calibrated lo-res IRS spectra
;
; INPUTS:
;
;	FNAME - 	tag of object to reduce (eg, 'mega005')
;
; OUTPUTS:
;
;	TAU16, TAU16ERR		- maximum optical depth and error of crystalline silicate at 16 um
;
;	TAU23, TAU23ERR		- maximum optical depth and error of crystalline silicate at 23 um
;
; KEYWORDS:
;
;	SILMIN16, SILMIN23	- two-element array giving the rest wavelength [um] in which to measure the depth of the
;					silicate features. Ex: [15.5,16.5]
;
;	AMORPH_SIL		- array with pivot points for secondary spline fit to continuum. 
;					Default is [10, 10.5, 12, 12.5, 13, 13.5, 14.5, 15, 16.8, 17.5, 20, 21, 25, 30]
;
;	MANUAL			- set to override the latest saved values (if they exist) 
;					for AMORPH_SIL, SILMIN16, SILMIN23
;
;	SMOOTHED		- set to smooth spectra with a Gaussian kernel at FWHM = 0.02 um
;
;	NO_TRUNCATION		- final spectrum (optical depth) is plotted with all residuals greater than 0.1 
;					cut off (easier to see features on a single plot)
;
;	SAVE_SETTINGS		- for objects in which AMORPH_SIL, SILMIN16, or SILMIN23 must be adjusted from
;					the default values, the SAVE keyword saves the values and uses them
;					as the new default the next time CRYSTALLINE.pro is run. Overridden
;					by the MANUAL keyword. 
;
; EXAMPLE:
;
;	IDL> crystalline, 'mega001', /ps
;
; REQUIRES:
;
;
;
; NOTES:
;
;	Based on results and continuum fits in Spoon et al. (2006), ApJ, 638, 759
;
; REVISION HISTORY
;       Written by K. Willett                Sep 09
;	Added outputs TAU16, TAU23 - Feb 10
;-

pro crystalline, fname, $
	tau16, tau16err, tau23, tau23err, $
	silmin16 = silmin16, silmin23 = silmin23, $
	amorph_sil = amorph_sil, $
	savesettings = savesettings, manual = manual, $
	ps = ps, noplot = noplot, stop = stop, $
	smoothed = smoothed, $
	no_truncation = no_truncation, $
	quiet = quiet, yr=yr, xr=xr

;device, window_state = state

; Read in the full LR spectrum

tag, fname, dirtag
saveset_file = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/silicate/'+fname+'_cryssil.sav'
restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr
err = sed.err_lr

; If best fits are already made, restore from IDL sav file

saveset_test = file_test(saveset_file)

if keyword_set(manual) or saveset_test eq 0 then begin
	if n_elements(silmin16) eq 0 then silmin16 = [15.0,17.0]
	if n_elements(silmin23) eq 0 then silmin23 = [22.0,24.0]
	
	; Remove amorphous silicate features with a second spline fit
	
	if n_elements(amorph_sil) le 0 then begin
		if wave[n_elements(wave)-1] gt 30. then $
			amorph_sil = [10, 10.5, 12, 12.5, 13, 13.5, 14.5, 15, 16.8, 17.5, 20, 21, 25, 30, 33] else $
			amorph_sil = [10, 10.5, 12, 12.5, 13, 13.5, 14.5, 15, 16.8, 17.5, 20, 21, 25, 30]
	endif
	
endif else restore, saveset_file

if keyword_set(smoothed) then flux = gauss_smooth(wave, flux, 0.2)

; Remove negative pixels, since they represent an unphysical measurement and produce
; wacky errors in a log-log plot

posflux = where(flux gt 0)

wave = wave[posflux]
flux = flux[posflux]
err  = err[posflux]

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

tau16 = min(cryst_resid(closeto(wave,silmin16[0]):closeto(wave,silmin16[1])))
tau23 = min(cryst_resid(closeto(wave,silmin23[0]):closeto(wave,silmin23[1])))

tau16err = stddev(cryst_resid(closeto(wave,15.5):closeto(wave,16.5)))
tau23err = stddev(cryst_resid(closeto(wave,22.5):closeto(wave,23.5)))
	
; Plot the results

if not keyword_set(noplot) then begin

	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/silicate/'+fname+'_cryst_sil.ps', /color
		cs = 1
		ls = 2
	endif else begin
		cs = 1
		ls = 1
	endelse

       cs = 0.5
       th = 1
       
	!p.multi = [0,1,3]
       
	if not keyword_set(xr) then xr = [8,40]
	if not keyword_set(yr) then yr = [1d-3,1d-1]

	; Raw spectra and local continuum fit

	plot, wave, flux, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		/xlog, $
		/ylog, $
		/xstyle, $
		charsize = cs, $
		thick = ls, $
		xthick = ls, $
		ythick = ls, $
		charthick = ls, $
		title = sed.obj
	
	oplot, wcont, fcont, psym = 4, symsize = 2, color = fsc_color("Yellow"), thick = ls
	oplot, wave, splresult, color = fsc_color("Yellow"), thick = ls
	
	; Optical depth spectrum and amorphous silicate fit
	
	plot, wave, tau, $
;		/xlog, $
		xr = [8,40], $
;		xr = xr, $
		/xstyle, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Optical depth !7s!3', $
		psym = 10, $
		charsize = cs, $
		thick = ls, $
		xthick = ls, $
		ythick = ls, $
		charthick = ls, $
		title = sed.obj
	
	hor, 0, linestyle = 1, thick = ls
	
	oplot, wcont_amorph, tcont_amorph, psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	oplot, wave, splresult_amorph, color = fsc_color("Red"), thick = ls

	; Residual optical depth after amorphous silicate subtraction

	if keyword_set(no_truncation) then plotind = indgen(n_elements(wave)) else plotind = where(cryst_resid lt 0.1)

	plot, wave[plotind], cryst_resid[plotind], $
;		/xlog, $
		xr = [8,40], $
		yr=[-0.8,0.5], $
		/xstyle, $
;		xtitle = 'Wavelength [!7l!3m]', $
;		ytitle = 'Residual optical depth', $
;		title = sed.obj, $
		psym = 10, $
		charsize = cs, $
		thick = th, $
		charthick = th
	

	; Label poss. crystalline silicate features

	ver, 11, linestyle=2 & xyouts, 11.1, 0.4, '11', /data
	ver, 16, linestyle=2 & xyouts, 16.1, 0.4, '16', /data
	ver, 19, linestyle=2 & xyouts, 19.1, 0.4, '19', /data
	ver, 23, linestyle=2 & xyouts, 23.1, 0.4, '23', /data
	ver, 28, linestyle=2 & xyouts, 28.1, 0.4, '28', /data

	hor, 0, linestyle=1

	ind10 = closeto(wave,10)
	ind30 = closeto(wave,30)
;	hor, min(cryst_resid[ind10:ind30]), color=fsc_color("Dodger Blue"), linestyle=2, thick = th

	ver, silmin16(0), color=fsc_color("Dodger Blue")
	ver, silmin16(1), color=fsc_color("Dodger Blue")
	
	ver, silmin23(0), color=fsc_color("Dodger Blue"), linestyle=2
	ver, silmin23(1), color=fsc_color("Dodger Blue"), linestyle=2
	
	hor, tau16, color=fsc_color("Dodger Blue")
	hor, tau23, color=fsc_color("Dodger Blue"), linestyle=2
	
	xyouts, 33, 0.2, /data, sed.obj, charsize=1.0, charthick=th

;	!p.multi=[0,1,1]

;	if keyword_set(ps) then begin
;		device, /close
;		set_plot,'x'
;	endif

endif

if not keyword_set(quiet) then begin
	print,''
	print, '16 um optical depth: ', string(tau16,format='(f5.2)'),' +- ',string(tau16err,format='(f5.3)'),' for ',sed.obj
	print,''
	print, '23 um optical depth: ', string(tau23,format='(f5.2)'),' +- ',string(tau23err,format='(f5.3)')
	print,''
endif
	
;savefile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/silicate/'+fname+'_sil.sav'
;sil_wave = wave & sil_spline = splresult
;save, filename=savefile, sil_wave, sil_spline

if keyword_set(savesettings) then begin
	save, filename=saveset_file, amorph_sil, silmin16, silmin23
endif

if keyword_set(stop) then stop

end
