;+
; NAME:
;       SILICATE
;
; PURPOSE:
;	Measure the silicate 9.7 and 18 um optical depths from calibrated lo-res IRS spectra
;
; INPUTS:
;
;	FNAME - 	tag of object to reduce (eg, 'mega005')
;
;	SPL - 		integer giving type of spline fit for used for the continuum
;
;			1: continuum-dominated - power-law from 5.0-7.5 and 27.0-31.5 um, intermediate spline anchor at 14.0 um
;			2: PAH-dominated - power-law from 5.5-14.5 um, spline anchors at 14.5 and 27.0 um
;			3: absorption-dominated (DEFAULT) - spline fit anchored at 5.2, 5.6, 14.0, and 27.0 um
;			4: absorption-dominated w/no PAH flux at 8 um (adds anchor at 7.8 um)
;			5+: custom spline at user-defined pivot points 
;
;	SILMIN10 - 	two-element array giving region in which to find the local minimum for the 9.7 um silicate
;				strength. The default is to look between 8.0 and 12.0 um - this option is
;				useful for avoiding noise spikes that give a false, deeper depth than actually
;				exists (especially for objects with very weak continuum). 
;
;	SILMIN18 - 	two-element array giving region in which to find the local minimum for the 18 um silicate
;				strength. The default is to look between 17.0 and 19.0 um.
;
;	- calibrated lo-res spectra from all four modules via SPEC2PUFLUX
;
; OUTPUTS:
;
;	TAU10 - 	optical depth of 9.7 um silicate feature
;
;	TAU10ERR - 	error in optical depth of 9.7 um silicate feature
;
;	TAU18 - 	optical depth of 18 um silicate feature
;
;	TAU18ERR - 	error in optical depth of 18 um silicate feature
;
;	WAVE - 		array w/wavelength values of spectrum (saved to IDL file)
;
;	SPLRESULT - 	spline fit to continuum (saved to IDL file)
;
; KEYWORDS:
;
;	PS - hard copies of spectra with spline fits and continuum removed
;
;	SILMAX10 - if set, then measures the 9.7 um silicate emission (common in PG quasars) rather than absorption
;
;	SILMAX18 - if set, then measures the 18 um silicate emission 
;
; EXAMPLE:
;
;	IDL> silicate, /ps
;
; REQUIRES:
;
;	TAG.pro
;	READCOL.pro
;	SPLINEFIT.pro
;
; NOTES:
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 2007
;	Added bonus order - KW, Aug 07
;	Only does silicate measurement now; PAH fits done in PAHSPL.pro - KW, Sep 07
;	Added SILMIN input - KW, Nov 07
;	Output the spline fit for use in PAH EW ice correction (not really necessary; recomputed in PAHSPL_NU) - Aug 08
;	Mid-IR spline point at 26.1 um moved to 27.0 um to avoid OIV/FeII emission - Jul 09
;	Added measurement of 18 um silicate feature - Jul 09
;	Made SPLINEFIT standalone routine - Jan 10
;-

pro silicate, fname, $
	tau10, tau10err, tau18, tau18err, $
	ps = ps, spl = spl, noplot = noplot, stop = stop, $
	silmin10 = silmin10, silmax10 = silmax10, silmin18 = silmin18, silmax18 = silmax18, $
	quiet = quiet, yr=yr, xr=xr, nosave = nosave

device, window_state = state

; Read in the full LR spectrum

tag, fname, dirtag
restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr

; Remove negative pixels, since they represent an unphysical measurement and produce
; wacky errors in a log-log plot

negflux = where(flux le 0)
posflux = setdifference(indgen(n_elements(flux)),negflux)

wave = wave(posflux)
flux = flux(posflux)

; Spline fits for the continuum

if not keyword_set(spl) then spl = 2

splinefit, spl=spl, wave, flux, splresult, pivotpoints

tau = alog(flux / splresult)

if n_elements(silmin10) eq 0 then silmin10 = [8.0,12.0]
if n_elements(silmin18) eq 0 then silmin18 = [17.0,19.0]

if keyword_set(silmax10) then tau10 = max(tau(closeto(wave,silmin10(0)):closeto(wave,silmin10(1)))) else $
	tau10 = min(tau(closeto(wave,silmin10(0)):closeto(wave,silmin10(1))))
if keyword_set(silmax18) then tau18 = max(tau(closeto(wave,silmin18(0)):closeto(wave,silmin18(1)))) else $
	tau18 = min(tau(closeto(wave,silmin18(0)):closeto(wave,silmin18(1))))

tau10err = stddev(tau(closeto(wave,9.5):closeto(wave,9.9)))
tau18err = stddev(tau(closeto(wave,17.0):closeto(wave,19.0)))
	
;if fname eq 'cso008' then tau10err = stddev([tau(closeto(wave,8.5):closeto(wave,9.5)),tau(closeto(wave,9.9):closeto(wave,10.5))])

; Plot the results

if not keyword_set(noplot) then begin

	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/silicate/'+fname+'_sil.ps', /color
		cs = 1
		ls = 2
	endif else begin
		cs = 2
		ls = 1
	endelse
	
	!p.multi = [0,1,2]
	
	if not keyword_set(xr) then xr = [4,40]
	if not keyword_set(yr) then yr = [1d-3,1d-1]

	plot, wave, flux, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		title = sed.obj, $
		/xlog, $
		/ylog, $
		xrange = xr, $
		yrange = yr, $
		/xstyle, $
		charsize = cs, $
		thick = ls, $
		charthick = ls
	
	oplot, wave(pivotpoints), flux(pivotpoints), psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	oplot, wave, splresult, color = fsc_color("Yellow"), thick = ls
	
	xyouts,5,1,charsize=1.5,fname
	
	plot, wave, tau, $
		/xlog, $
		xr = xr, $
		/xstyle, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Optical depth !7s!3', $
		title = sed.obj, $
		charsize = cs, $
		thick = ls, $
		charthick = ls
	
	oplot, wave(pivotpoints), tau(pivotpoints), psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	hor, 0, linestyle = 1, thick = ls
	xyouts, 0.8, 0.95, fname, /normal, charsize = cs, charthick = ct, color=axiscolor
	
	ver, silmin10(0), color=fsc_color("Dodger Blue")
	ver, silmin10(1), color=fsc_color("Dodger Blue")
	
	ver, silmin18(0), color=fsc_color("Dodger Blue"), linestyle=2
	ver, silmin18(1), color=fsc_color("Dodger Blue"), linestyle=2
	
	hor, tau10, color=fsc_color("Red")
	hor, tau18, color=fsc_color("Red"), linestyle=2
	
	!p.multi=[0,1,1]

	if keyword_set(ps) then begin
		device, /close
		set_plot,'x'
	endif

endif

if not keyword_set(quiet) then begin
	print,''
	print, '9.7 um silicate: ', tau10,' +- ',string(tau10err,format='(f5.3)'),' for ',fname
	print,''
	print, '18  um silicate: ', tau18,' +- ',string(tau18err,format='(f5.3)'),' for ',fname
endif
	
if not keyword_set(nosave) then begin
	savefile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/silicate/'+fname+'_sil.sav'
	sil_wave = wave & sil_spline = splresult
	save, filename=savefile, sil_wave, sil_spline
endif

if keyword_set(stop) then stop

end
