;+
; NAME:
;       PAHSPL_NU
;
; PURPOSE:
;	Measure the 6.2 and 11.3 um PAH features, giving their fluxes, EWs, and errors. This program
;		measures the PAH flux using a spline fit. 
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
;			5: custom spline at user-defined pivot points 
;
;	PIVOTS_62 - 	pivot points to define continuum for spline fit of 6.2 um PAH
;
;	PIVOTS_11 - 	pivot points to define continuum for spline fit of 11.3 um PAH
;
;	- calibrated lo-res spectra from all four modules via SPEC2PUFLUX
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
;	IDL> 
;
; REQUIRES:
;
;	TAG.pro
;	READCOL,pro
;
; NOTES:
;
;	PAHSPL.pro, the previous incarnation of this program, did not properly take into account the 
;	difference between f_nu and f_lambda wrt Janskys, and should not be used. - KW, Oct 28, 2007.
;
; REVISION HISTORY
;       Written by K. Willett                Sep 2007
;	Edited to get units of flux/EW right		- KW, Oct 07
;	Added inputs for setting continuum pivot points at 6.2 and 11.3 um - KW, Nov 07
;	Added proper measurements for getting the water-ice corrected 6.2 EW - Aug 08
;	Added upper limits - Aug 08
;	Added 6th SPL option - Jan 10
;-

pro pahspl_nu, fname, $
	pahew62, pahew62err, $
	pahew62_ice, pahew62_iceerr, $
	pahew11, pahew11err, $
	pah62flux, pah11flux, $
	pah62flux_lim, pahew62_lim, $
	pah11flux_lim, pahew11_lim, $
	pah62fwhm, pah11fwhm, $
	spl = spl, ps = ps, pivots_62 = pivots_62, pivots_11 = pivots_11, $
	noplot = noplot, quiet = quiet, stop = stop, wait = wait, $
	xr62 = xr62, yr62 = yr62, xr11 = xr11, yr11 = yr11, $
	limit = limit, cont62=cont62

if not keyword_set(ps) then device, window_state = state

limfac = 3.			; Factor by which to multiply positive PAH flux (for a conservative upper limit)
c = 299792.458d * 1d9		; um / s

; Read in the full LR spectrum

tag,fname,dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr
err = sed.err_lr

; Spline fits for the 9.7 um continuum

if not keyword_set(spl) then spl = 2

splinefit, spl=spl, wave, flux, splresult, pivotpoints

; Add feature to measure the PAH 6.2 um EW for placement in Spoon's classification scheme

; Designate continuum regions on either side

if n_elements(pivots_62) eq 0 then begin
	if closeto(wave,5.15) ne 0 then pahpiv = [closetomed(wave,flux,5.15),closetomed(wave,flux,5.55),$
		closetomed(wave,flux,5.95),closetomed(wave,flux,6.55),closetomed(wave,flux,7.1)] else $
		pahpiv = [closetomed(wave,flux,5.55),closetomed(wave,flux,5.95),closetomed(wave,flux,6.55),closetomed(wave,flux,7.1)]
endif else begin
	npiv = n_elements(pivots_62)
	pahpiv = fltarr(npiv)
	for j = 0,npiv-1 do pahpiv(j) = closetomed(wave,flux,pivots_62(j))
endelse


bpah = closetomed(wave,flux,5.95)
epah = closetomed(wave,flux,6.55)

; Spline fit to PAH continuum

wp = wave(pahpiv)
fp = flux(pahpiv)
pspl1 = spl_init(wp,fp)
pahspl = spl_interp(wp,fp,pspl1,wave)

; Add the flux over designated area

newpah = flux - pahspl								; Continuum subtracted flux array
addpah = fltarr(2,epah - bpah + 1)
for i = bpah, epah do begin
	addpah(0,i-bpah) = newpah(i)
	addpah(1,i-bpah) = (wave(i+1) - wave(i)) * c / (wave(i))^2
endfor

; 6.2 PAH flux and error

pah62fluxarr = addpah(0,*) * addpah(1,*) 					; Flux density in Jy, d_nu in Hz
pah62flux = total(pah62fluxarr) * 1d-30						; Flux in W/cm^2

cont62 = pahspl(closetomed(wave,flux,6.2))
pahew62 = pah62flux / cont62 * (6.2)^2 / c * 1d30				; EW in um
baseerr62 = stddev(flux(closetomed(wave,flux,6.1):closetomed(wave,flux,6.3)))	; Estimate error in baseline flux density [Jy]
pahew62err = baseerr62 * pah62flux / cont62^2 * (6.2)^2/ c*1d30			; Error in equivalent width

; 6.2 um FWHM

peak62 = max(newpah[bpah:epah])
peakwave62 = where(newpah[bpah:epah] eq peak62)
peakind62 = peakwave62 + bpah

i = 0 & nplevel = peak62
while nplevel gt 0.5 * peak62 do begin
	i = i+1
	nplevel = newpah[peakind62 - i]
endwhile

if i gt 0 then pah62fwhm = wave[peakind62 + i] - wave[peakind62 - i] else pah62fwhm = 0.

; Also define a PAH EW using the spline-based silicate continuum (correcting for water ice)

ice_62cont = splresult(closetomed(wave,splresult,6.2)) 
pahew62_ice = pah62flux / ice_62cont * (6.2)^2 / c * 1d30
baseerr_spl = stddev(splresult(closetomed(wave,splresult,6.1):closetomed(wave,splresult,6.3)))
pahew62_iceerr = baseerr_spl * pah62flux / ice_62cont^2 * (6.2)^2 / c*1d30

; Upper limits for PAH 6.2

if keyword_set(limit) then begin

	limrms = stddev(newpah[bpah:epah])
	fwhm62 = 0.14				; um - taken from average of all detected 6.2 um features in arch sample

	flux_cgs = limfac * limrms / (6.2)^2 * c * 1d-30
	sigwidth = fwhm62 / (2d * sqrt(2d * alog(2d)))
	pah62flux_lim = flux_cgs * sigwidth * sqrt(2d * !dpi) 			; Units of W/cm^2

	pahew62_lim = pah62flux_lim / cont62 * (6.2)^2 / c * 1d30
endif

; Measure the PAH 11 um EW

; Fit continuum w/pivot points

if n_elements(pivots_11) eq 0 then pivots_11 = [10.1,10.9,11.8,12.4]

npiv_11 = n_elements(pivots_11)
pahpiv11 = fltarr(npiv_11)
for j = 0, npiv_11-1 do pahpiv11(j) = closetomed(wave,flux,pivots_11(j))

bpah11 = closetomed(wave,flux,10.9)
epah11 = closetomed(wave,flux,11.8)

; Spline fit to PAH continuum

wp11 = wave(pahpiv11)
fp11 = flux(pahpiv11)
pspl11 = spl_init(wp11,fp11)
pahspl11 = spl_interp(wp11,fp11,pspl11,wave)

newpah11 = flux - pahspl11
addpah11 = fltarr(2,epah11 - bpah11 + 1)
for i = bpah11, epah11 do begin
	addpah11(0,i-bpah11) = newpah11(i)
	addpah11(1,i-bpah11) = (wave(i+1) - wave(i)) * c / (wave(i))^2
endfor

pah11fluxarr = addpah11(0,*) * addpah11(1,*) 
pah11flux = total(pah11fluxarr) * 1d-30
cont11 = pahspl11(closetomed(wave,flux,11.3))
pahew11 = pah11flux / cont11 * (11.3)^2 / c * 1d30
baseerr11 = stddev(pahspl11(closetomed(wave,flux,11.2):closetomed(wave,flux,11.4)))
pahew11err = baseerr11 * pah11flux / (pahspl11(closetomed(wave,flux,11.3)))^2 * (11.3)^2/ c*1d30

; 11.3 um FWHM

peak11 = max(newpah11[bpah11:epah11])
peakwave11 = where(newpah11[bpah11:epah11] eq peak11)
peakind11 = peakwave11 + bpah11

i = 0 & nplevel11 = peak11
while nplevel11 gt 0.5 * peak11 do begin
	i = i+1
	nplevel11 = newpah11[peakind11 - i]
endwhile

if i gt 0 then pah11fwhm = wave[peakind11 + i] - wave[peakind11 - i] else pah11fwhm = 0.

; Upper limits for PAH 11.3

if keyword_set(limit) then begin
	limrms = stddev(newpah[bpah:epah])
	fwhm11 = 0.24				; um - taken from average of all detected 6.2 um features in arch sample

	c = 299792.458d * 1d9
	flux_cgs = limfac * limrms / (11.3)^2 * c * 1d-30
	sigwidth = fwhm11 / (2d * sqrt(2d * alog(2d)))				; Assumes a Gaussian, not Voigt (not _quite_ correct)
	pah11flux_lim = flux_cgs * sigwidth * sqrt(2d * !dpi) 			; Units of W/cm^2

	pahew11_lim = pah11flux_lim / cont11 * (11.3)^2 / c * 1d30
endif

; Plot the results

if not keyword_set(noplot) then begin
	
	!p.multi = [0,2,1]
	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/pah/'+fname+'_pahsplnu.ps', /color
		cs = 1
		ls = 2
	endif else begin
		cs = 2
		ls = 1
	endelse
	
	if n_elements(xr62) eq 0 then xr62 = [5,8]
	if n_elements(yr62) eq 0 then yr62 = [1d-4,1d0]
	
	; 6.2 um PAH
	
	plot, wave, flux, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		psym=10, $
	;	/xlog, $
		/ylog, $
		xrange = xr62, $
		yrange = yr62, $
		/xstyle, $
		/ystyle, $
		charsize = cs, $
		thick = ls, $
		charthick = ls, $
		title = sed.obj
	
	oplot, wp, fp, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
	oplot, wave, pahspl, color=fsc_color("Green"), thick = ls
	oplot, wave, splresult, color=fsc_color("Yellow"), thick = ls
	ver, wave(bpah), linestyle = 1
	ver, wave(epah), linestyle = 1
	xyouts, 0.5,0.95, fname, /normal
	
	;xyouts, 5.2, 0.07, 'PAH 6.2 um EW  : '+string(pahew62,format='(f6.3)')+' um', /data, charsize = cs
	;xyouts, 5.2, 0.05, 'PAH 6.2 um flux : '+string(pah62flux,format='(e9.3)')+' W/cm^2', /data, charsize = cs
	
	; 11.3 um PAH
	
	if n_elements(xr11) eq 0 then xr11 = [8,15]
	if n_elements(yr11) eq 0 then yr11 = [1d-4,1d1]
	
	plot, wave, flux, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		psym=10, $
	;	/xlog, $
		/ylog, $
		xrange = xr11, $
		yrange = yr11, $
		/xstyle, $
		charsize = cs, $
		thick = ls, $
		charthick = ls, $
		title = sed.obj
	
	oplot, wp11, fp11, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
	oplot, wave, pahspl11, color=fsc_color("Green"), thick = ls
	ver, wave(bpah11), linestyle = 1
	ver, wave(epah11), linestyle = 1
	
	;xyouts, 5.2, 0.7, 'PAH 11um EW  : '+string(pahew11,format='(f6.3)')+' um', /data, charsize = cs
	;xyouts, 5.2, 0.5, 'PAH 11um flux : '+string(pah11flux,format='(e9.3)')+' W/cm^2', /data, charsize = cs
	
	if keyword_set(ps) then begin
		device,/close
		set_plot,'x'
	endif
	
	!p.multi=[0,1,1]

endif	; noplot


; Print results to screen

if not keyword_set(quiet) then begin

	print,''
	print,'PAH 6.2 EW [um]      :  ',string(pahew62,format='(f6.3)'),' +- ',string(pahew62err,format='(f6.3)')
	print,'PAH 6.2 EW [um] (ice):  ',string(pahew62_ice,format='(f6.3)'),' +- ',string(pahew62_iceerr,format='(f6.3)')
	print,'PAH 11.3 EW [um]     :  ',string(pahew11,format='(f6.3)'),' +- ',string(pahew11err,format='(f6.3)')
	print,''
	print,'PAH 6.2 flux [W/cm^2] :  ',string(pah62flux,format='(e8.2)')
	if keyword_set(limit) then print,'PAH 6.2 flux limit [W/cm^2] :  ',string(pah62flux_lim,format='(e8.2)')
	print,'PAH 11.3 flux [W/cm^2]:  ',string(pah11flux,format='(e8.2)')
	if keyword_set(limit) then print,'PAH 11.3 flux limit [W/cm^2] :  ',string(pah11flux_lim,format='(e8.2)')

endif

if n_elements(wait) gt 0 then wait,wait
if keyword_set(stop) then stop

end
