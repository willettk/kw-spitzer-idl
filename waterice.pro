;+
; NAME:
;
;       WATERICE
;
; PURPOSE:
;	
;	Measure the water ice depth using the silicate-based spline continuum
;
; INPUTS:
;
;	FNAME - 	tag of object to reduce (eg, 'mega005')
;
;	SPL - 		integer giving type of spline fit for used for the 9.7 um continuum
;
;			1: continuum-dominated
;			2: PAH-dominated
;			3: absorption-dominated (DEFAULT)
;			4: absorption-dominated w/no PAH flux at 8 um
;			5+: custom spline at user-defined pivot points 
;
; OUTPUTS:
;
; KEYWORDS:
;
;	PS - hard copies of spectra with spline fits and continuum removed
;
;	QUIET - do not print results to screen (default shows fluxes and EW)
;
; EXAMPLE:
;
;	IDL> waterice, 'mega001', depth, err
;
; REQUIRES:
;
;	TAG.pro
;	READCOL,pro
;
; NOTES:
;
;	Adapted from PAHSPL_NU.pro
;
; REVISION HISTORY
;       
;	Written by KW, Dec 08
;	Added 6th SPL option - Jan 10
;-

pro waterice, fname, icedepth, icedepth_err, $
	ps = ps, spl = spl, noplot = noplot, stop = stop, icemin = icemin, icemax = icemax, wait = wait, quiet = quiet

device, window_state = state

; Read in the full LR spectrum

tag, fname,dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr
err = sed.err_lr

; Remove negative pixels, since they represent an unphysical measurement and produce
; wacky errors in a log-log plot

negflux = where(flux le 0)
posflux = setdifference(indgen(n_elements(flux)),negflux)

wave = wave(posflux)
flux = flux(posflux)
err = err(posflux)

; Spline fits for the 9.7 um continuum

if not keyword_set(spl) then spl = 2

splinefit, spl=spl, wave, flux, splresult, pivotpoints

; Plot the results

;if state(0) eq 1 then wset, 0 else window, 0

tau = alog(flux / splresult)

if n_elements(icemin) eq 0 then icemin = [5.5,6.3]

if keyword_set(icemax) then tauice = max(tau(closeto(wave,icemin(0)):closeto(wave,icemin(1)))) else $
	tauice = min(tau(closeto(wave,icemin(0)):closeto(wave,icemin(1))))
tauiceerr = stddev(tau(closeto(wave,5.8):closeto(wave,6.1)))
	
if not keyword_set(noplot) then begin

	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/waterice/'+fname+'_ice.ps', /color
		cs = 1
		ls = 2
	endif else begin
		cs = 2
		ls = 1
	endelse
	
	!p.multi = [0,1,2]
	
	plot, wave, flux, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
	;	/xlog, $
		/ylog, $
		xrange = [5,8], $
		yrange = [1d-4,1d1], $
		/xstyle, $
		charsize = cs, $
		thick = ls, $
		charthick = ls, $
		title = sed.obj
	
	oplot, wave(pivotpoints), flux(pivotpoints), psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	oplot, wave, splresult, color = fsc_color("Yellow"), thick = ls
	
	plot, wave, tau, $
	;	/xlog, $
		xr = [5,8], $
		/xstyle, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Optical depth !7s!3', $
		charsize = cs, $
		thick = ls, $
		charthick = ls, $
		title = sed.obj
	
	oplot, wave(pivotpoints), tau(pivotpoints), psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	hor, 0, linestyle = 1, thick = ls
	xyouts, 0.8, 0.95, fname, /normal, charsize = cs, charthick = ct, color=axiscolor
	
	ver, icemin(0), color=fsc_color("Dodger Blue")
	ver, icemin(1), color=fsc_color("Dodger Blue")
	
	hor, tauice, color=fsc_color("Red")
	
	if not keyword_set(quiet) then begin
		print,''
		print, 'Water ice tau: ', string(tauice, format='(f5.2)'),' +- ',string(tauiceerr,format='(f5.2)'),' for ',fname
	endif
	
	!p.multi=[0,1,1]
	
endif

;savefile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/waterice/'+fname+'_waterice.sav'
;ice_wave = wave & ice_spline = splresult
;save, filename=savefile, ice_wave, ice_spline

icedepth = tauice
icedepth_err = tauiceerr

if keyword_set(stop) then stop
if keyword_set(wait) then wait, 3
end
