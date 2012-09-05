pro sirocky, fname, tau10, tau10err, tau18, tau18err, $
	ps = ps, spl = spl, noplot = noplot, stop = stop, quiet = quiet, $
	silmin10 = silmin10, silmax10 = silmax10, silmin18 = silmin18, silmax18 = silmax18
;+
; NAME:
;       SIROCKY
;
; PURPOSE:
;	Measure the silicate 10 and 18 um optical depths from calibrated lo-res IRS spectra
;
;
; INPUTS:
;
;	FNAME - 	tag of object to reduce (eg, 'mega005')
;
;	SPL - 		integer giving type of spline fit for used for the um continuum
;
;			1: continuum-dominated - spline fit to 5-7, 13.8-14.2, and 26.5-31.5 um
;			2: absorption-dominated - spline fit to 5.2-5.6, 7.8, 13.2-14.5, 26.0-31.5
;			3: PAH-dominated - power-law fit from mean(5.3-5.7) to mean(14.0-15.0), spline from 26.0-31.5 (DEFAULT)
;
;	SILMIN10 - 	two-element array giving region in which to find the local minimum for the 10 um silicate strength
;
; OUTPUTS:
;
;	TAU10 - 	optical depth of 10 um silicate feature
;
;	TAU10ERR - 	error in optical depth of 10 um silicate feature
;
;	TAU18 - 	optical depth of 18 um silicate feature
;
;	TAU18ERR - 	error in optical depth of 18 um silicate feature
;
; KEYWORDS:
;
;	PS - hard copies of spectra with spline fits and continuum removed
;
;	SILMAX - if set, then measures the silicate emission (common in PG quasars) rather than absorption
;
; EXAMPLE:
;
;	IDL> sirocky, 'mega001', spl=3, /quiet
;
; REQUIRES:
;
;	TAG.pro
;	READCOL,pro
;
; NOTES:
;
;	Spline-fit method is from Sirocky et al. 2008
;
; REVISION HISTORY
; 	Adapted from SILICATE.pro 			Apr 09
;-

; Read in the full LR spectrum

tag, fname,dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr
err = sed.err_lr

; Remove negative pixels, since they represent an unphysical measurement and produce
; wacky errors in a log-log plot

wave = wave[where(flux ge 0)]
flux = flux[where(flux ge 0)]
err = err[where(flux ge 0)]

; Spline fits for the continuum

if not keyword_set(spl) then spl = 3
case spl of
	;
	; Continuum-dominated
	;
	1: begin
		if wave[0] gt 5.0 then begwave = 5.5 else begwave = 5.0
		if wave[n_elements(wave)-1] lt 31.5 then endwave = 28.0 else endwave = 31.5
		
		; Resample spectrum

		seg1_wave = wave[closeto(wave,begwave):closeto(wave,7.0)]
		seg1_flux = flux[closeto(wave,begwave):closeto(wave,7.0)]
		size1 = size(seg1_wave)
		binsize1 = 10
;		seg1_wave_r = seg1_wave[indgen(size1[1] / binsize1 + 1) * binsize1]
;		seg1_flux_r = seg1_flux[indgen(size1[1] / binsize1 + 1) * binsize1]
		seg1_wave_r = [seg1_wave[0],seg1_wave[n_elements(seg1_wave)-1]]
		seg1_flux_r = [seg1_flux[0],seg1_flux[n_elements(seg1_flux)-1]]

		seg2_wave = wave[closeto(wave,13.8):closeto(wave,14.2)]
		seg2_flux = flux[closeto(wave,13.8):closeto(wave,14.2)]
		size2 = size(seg2_wave)
		binsize2 = 2
;		seg2_wave_r = seg2_wave[indgen(size2[1] / binsize2 + 1) * binsize2]
;		seg2_flux_r = seg2_flux[indgen(size2[1] / binsize2 + 1) * binsize2]
;		seg2_wave_r = [seg2_wave[0],seg2_wave[n_elements(seg2_wave)-1]]
;		seg2_flux_r = [seg2_flux[0],seg2_flux[n_elements(seg2_flux)-1]]
		seg2_wave_r = wave[closetomed(wave,flux,14.0)]
		seg2_flux_r = flux[closetomed(wave,flux,14.0)]

		seg3_wave = wave[closeto(wave,26.5):closeto(wave,endwave)]
		seg3_flux = flux[closeto(wave,26.5):closeto(wave,endwave)]
		size3 = size(seg3_wave)
		binsize3 = 4
;		seg3_wave_r = seg3_wave[indgen(size3[1] / binsize3 + 1) * binsize3]
;		seg3_flux_r = seg3_flux[indgen(size3[1] / binsize3 + 1) * binsize3]
		seg3_wave_r = [seg3_wave[0],seg3_wave[n_elements(seg3_wave)-1]]
		seg3_flux_r = [seg3_flux[0],seg3_flux[n_elements(seg3_flux)-1]]

		wcont = [seg1_wave_r, seg2_wave_r, seg3_wave_r]
		fcont = [seg1_flux_r, seg2_flux_r, seg3_flux_r]
		yspl = spl_init(wcont,fcont)
	   	splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,begwave),closetomed(wave,flux,7.0),$
;			closetomed(wave,flux,13.8),closetomed(wave,flux,14.2),$
			closetomed(wave,flux,14.0),$
			closetomed(wave,flux,26.5),closetomed(wave,flux,endwave)]
	   end
	;
	; Absorption-dominated
	;
	2: begin 
		if wave[0] gt 5.0 then begwave = 5.4 else begwave = 5.2
		if wave[n_elements(wave)-1] lt 31.5 then endwave = 28.0 else endwave = 31.5

		; Resample spectrum

		seg1_wave = wave[closeto(wave,begwave):closeto(wave,5.6)]
		seg1_flux = flux[closeto(wave,begwave):closeto(wave,5.6)]
		size1 = size(seg1_wave)
		binsize1 = 3
;		seg1_wave_r = seg1_wave[indgen(size1[1] / binsize1 + 1) * binsize1]
;		seg1_flux_r = seg1_flux[indgen(size1[1] / binsize1 + 1) * binsize1]
		seg1_wave_r = [seg1_wave[0],seg1_wave[n_elements(seg1_wave)-1]]
		seg1_flux_r = [seg1_flux[0],seg1_flux[n_elements(seg1_flux)-1]]

		seg2_wave = wave[closetomed(wave,flux,7.8)]
		seg2_flux = flux[closetomed(wave,flux,7.8)]

		seg3_wave = wave[closeto(wave,13.2):closeto(wave,14.5)]
		seg3_flux = flux[closeto(wave,13.2):closeto(wave,14.5)]
		size3 = size(seg3_wave)
		binsize3 = 5
;		seg3_wave_r = seg3_wave[indgen(size3[1] / binsize3 + 1) * binsize3]
;		seg3_flux_r = seg3_flux[indgen(size3[1] / binsize3 + 1) * binsize3]
		seg3_wave_r = [seg3_wave[0],seg3_wave[n_elements(seg3_wave)-1]]
		seg3_flux_r = [seg3_flux[0],seg3_flux[n_elements(seg3_flux)-1]]

		seg4_wave = wave[closeto(wave,26.5):closeto(wave,endwave)]
		seg4_flux = flux[closeto(wave,26.5):closeto(wave,endwave)]
		size4 = size(seg4_wave)
		binsize4 = 4
;		seg4_wave_r = seg4_wave[indgen(size4[1] / binsize4 + 1) * binsize4]
;		seg4_flux_r = seg4_flux[indgen(size4[1] / binsize4 + 1) * binsize4]
		seg4_wave_r = [seg4_wave[0],seg4_wave[n_elements(seg4_wave)-1]]
		seg4_flux_r = [seg4_flux[0],seg4_flux[n_elements(seg4_flux)-1]]

;		wcont = [seg1_wave_r, seg2_wave, seg3_wave_r, seg4_wave_r]
;		fcont = [seg1_flux_r, seg2_flux, seg3_flux_r, seg4_flux_r]
		wcont = [seg1_wave_r, seg3_wave_r, seg4_wave_r]
		fcont = [seg1_flux_r, seg3_flux_r, seg4_flux_r]
		yspl = spl_init(wcont,fcont)
	   	splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,begwave),closetomed(wave,flux,5.6),$
			closetomed(wave,flux,7.8), $
			closetomed(wave,flux,13.2),closetomed(wave,flux,14.5),$
			closetomed(wave,flux,26.5),closetomed(wave,flux,endwave)]
	   end
	;
	; PAH-dominated
	;
	3: begin
		if wave[0] gt 5.0 then begwave = 5.5 else begwave = 5.0
		if wave[n_elements(wave)-1] lt 31.5 then endwave = 28.0 else endwave = 31.5

		flux_powerlaw_lo = mean(flux[closeto(wave,5.3):closeto(wave,5.7)])
		flux_powerlaw_hi = mean(flux[closeto(wave,14.0):closeto(wave,15.0)])

		wave_powerlaw_lo = closeto(flux[closeto(wave,5.3):closeto(wave,5.7)],  flux_powerlaw_lo) + closeto(wave,5.3)
		wave_powerlaw_hi = closeto(flux[closeto(wave,14.0):closeto(wave,15.0)],flux_powerlaw_hi) + closeto(wave,14.0)

		m = alog10(flux_powerlaw_hi / flux_powerlaw_lo) $
			/ alog10(wave[wave_powerlaw_hi] / wave[wave_powerlaw_lo])
	   	b = alog10(flux_powerlaw_hi) - m * alog10(wave[wave_powerlaw_hi])
	   	powerfit = 10^(m * alog10(wave[0:wave_powerlaw_hi - 1]) + b)

		; Kluge for objects where spline is thrown off by features near 26 um

		arr26 = ['mega002','arch012','arch018']
		junk = where(arr26 eq fname, jcount)
		if jcount gt 0 then oIV_wave = 26.5 else oIV_wave = 27.5

		seg1_wave = wave[closeto(wave,oIV_wave):closeto(wave,endwave)]
		seg1_flux = flux[closeto(wave,oIV_wave):closeto(wave,endwave)]
;		seg1_wave = wave[closeto(wave,26.5):closeto(wave,endwave)]
;		seg1_flux = flux[closeto(wave,26.5):closeto(wave,endwave)]
		size1 = size(seg1_wave)
		binsize1 = 6
;		seg1_wave_r = seg1_wave[indgen(size1[1] / binsize1 + 1) * binsize1]
;		seg1_flux_r = seg1_flux[indgen(size1[1] / binsize1 + 1) * binsize1]
		seg1_wave_r = [seg1_wave[0],seg1_wave[n_elements(seg1_wave)-1]]
		seg1_flux_r = [seg1_flux[0],seg1_flux[n_elements(seg1_flux)-1]]
		seg1_wave_r = [seg1_wave[0],seg1_wave[fix(n_elements(seg1_wave)/2)], seg1_wave[n_elements(seg1_wave)-1]]
		seg1_flux_r = [seg1_flux[0],seg1_flux[fix(n_elements(seg1_wave)/2)], seg1_flux[n_elements(seg1_flux)-1]]

		wcont = [wave[wave_powerlaw_hi], seg1_wave_r]
		fcont = [flux_powerlaw_hi,       seg1_flux_r]
		yspl = spl_init(wcont, fcont)
		spltemp = spl_interp(wcont, fcont, yspl, wave[wave_powerlaw_hi:n_elements(wave)-1])
		splresult = [powerfit[0:n_elements(powerfit)-1], spltemp]
		pivotpoints = [closetomed(wave,flux,5.5),closetomed(wave,flux,14.5),$
			closeto(wave,26.5), closeto(wave,seg1_wave[fix(n_elements(seg1_wave)/2)]), closeto(wave,endwave)]
	   end
endcase

;plot, wave, flux, psym=10, /xlog,/ylog
;oplot, wave, splresult, color=fsc_color("Yellow")
;stop

; Measure depth of silicate features

tau = alog(flux / splresult)

	; 10 um

	if n_elements(silmin10) eq 0 then silmin10 = [8.0,12.0]
	
	if min(tau[closeto(wave,9.5):closeto(wave,10.0)]) gt 0. then $
		tau10 = max(tau[closeto(wave,silmin10[0]):closeto(wave,silmin10[1])]) else $
		tau10 = min(tau[closeto(wave,silmin10[0]):closeto(wave,silmin10[1])])
	
	tau10err = stddev(tau(closeto(wave,9.5):closeto(wave,9.9)))
	
	; 18 um

	if n_elements(silmin18) eq 0 then silmin18 = [16.0,20.0]
	
	if min(tau[closeto(wave,17.5):closeto(wave,18.5)]) gt 0. then $
		tau18 = max(tau[closeto(wave,silmin18[0]):closeto(wave,silmin18[1])]) else $
		tau18 = min(tau[closeto(wave,silmin18[0]):closeto(wave,silmin18[1])])
	
	tau18err = stddev(tau(closeto(wave,17.0):closeto(wave,19.0)))
	
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
	
	plot, wave, flux, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		/xlog, $
		/ylog, $
		xrange = [4,40], $
		yrange = [1d-3,1d1], $
		/xstyle, $
		charsize = cs, $
		thick = ls, $
		charthick = ls, $
		title = sed.obj+' - '+fname
	
	oplot, wave(pivotpoints), flux(pivotpoints), psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	oplot, wave, splresult, color = fsc_color("Yellow"), thick = ls
	
	plot, wave, tau, $
		/xlog, $
		xr = [4,40], $
		/xstyle, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Optical depth !7s!3', $
		charsize = cs, $
		thick = ls, $
		charthick = ls, $
		title = sed.obj+' - '+fname
	
	oplot, wave(pivotpoints), tau(pivotpoints), psym = 4, symsize = 2, color = fsc_color("Red"), thick = ls
	hor, 0, linestyle = 1, thick = ls
	
	ver, silmin10(0), color=fsc_color("Blue")
	ver, silmin10(1), color=fsc_color("Blue")
	hor, tau10, color=fsc_color("Blue")
	
	ver, silmin18(0), color=fsc_color("Dodger Blue")
	ver, silmin18(1), color=fsc_color("Dodger Blue")
	hor, tau18, color=fsc_color("Dodger Blue")
	
	!p.multi=[0,1,1]

endif

if not keyword_set(quiet) then begin
	print,''
	print, '10 um silicate: ', string(tau10,format='(f5.2)'),' +- ',string(tau10err,format='(f5.2)'),' for ',fname
	print,''
	print, '18 um silicate: ', string(tau18,format='(f5.2)'),' +- ',string(tau18err,format='(f5.2)'),' for ',fname
	print,''
endif
	
if keyword_set(stop) then stop

end
