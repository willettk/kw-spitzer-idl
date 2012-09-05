pro stitch_orders, fname, ps = ps, stop = stop, write = write, nobonus = nobonus, $
	nocal = nocal, nospicehdr = nospicehdr, quiet=quiet
;+
; NAME:
;       STITCH_ORDERS
;
; PURPOSE:
; 	Stitch together the different IRS modules and output them as files that
;	can be read by SMART. 
;
; INPUTS:
;
;
; KEYWORDS:
;
;	NOCAL - runs routine without calibration to add suffixes and read artificial fluxes in SMART
;
; EXAMPLE:
;	IDL> stitch_orders, 'mega001', /ps
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;	READCOL.pro
;	SPICEHDR.pro
;
; NOTES:
;
;	This program adapts previous versions of SPEC2PUFLUX_ARCH and SPICEHDR. Formerly known as SMARTMAKE.pro
;
; REVISION HISTORY
;       Written by K. Willett                Feb 2008
;	Default is now to do no scaling of the HR modules and stitch the LR together through
;		multiplicative scaling - KW, May 08
;	Renamed to STITCH_ORDERS.pro - May 08
;-

if not keyword_set(quiet) then begin
	print,''
	print, 'Running ',fname
endif

; Set device to read in colors

device, decomposed = 0

	yellow = fsc_color('Yellow')
	green = fsc_color('Green')
	red = fsc_color('Red')
	blue = fsc_color('Blue')
	orange = fsc_color('Orange')

; Find directory and object name

tag, fname, dirtag 
targets, fname, redshift, obj

; Read in the spectrum

	specpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/coadd/'

	; Lores

	readcol, specpath+fname+'_sl1_coadd.tbl', $
		det_sl1, wave_sl1, flux_sl1, err_sl1, bit_sl1, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_sl2_coadd.tbl', $
		det_sl2, wave_sl2, flux_sl2, err_sl2, bit_sl2, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_sl3_coadd.tbl', $
		det_sl3, wave_sl3, flux_sl3, err_sl3, bit_sl3, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll1_coadd.tbl', $
		det_ll1, wave_ll1, flux_ll1, err_ll1, bit_ll1, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll2_coadd.tbl', $
		det_ll2, wave_ll2, flux_ll2, err_ll2, bit_ll2, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll3_coadd.tbl', $
		det_ll3, wave_ll3, flux_ll3, err_ll3, bit_ll3, format = 'i,f,f,f,i', skipline = 1, /silent

	; Hires

	readcol, specpath+fname+'_sh_coadd.tbl', $
		det_sh, wave_sh, flux_sh, err_sh, bit_sh, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_lh_coadd.tbl', $
		det_lh, wave_lh, flux_lh, err_lh, bit_lh, format = 'i,f,f,f,i', skipline = 1, /silent

; Read in the filters

pupath = '~/Astronomy/Research/Spitzer/spitzer/'

readcol, pupath+'bluePUtrans.txt', bwave, btrans, /silent
readcol, pupath+'redPUtrans.txt', rwave, rtrans, /silent


; LL1 module (18-35 um) is fixed

calflux_ll1 = flux_ll1

; Anchor LL2 to LL1

	first_pixel = wave_ll1(0)
	last_pixel= wave_ll2(n_elements(wave_ll2)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_ll2(n_elements(flux_ll2)-npix:n_elements(flux_ll2)-1)
		first_sec = calflux_ll1(0:npix-1)

		ll2_ll1_frac = mean(first_sec) / mean(last_sec)	
		
		calflux_ll2 = flux_ll2 * ll2_ll1_frac
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_ll2 ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin
			junk1 = where(abs(wave_ll2(overlap_pixels(0)+i) - wave_ll1) eq $
				min(abs(wave_ll2(overlap_pixels(0)+i) - wave_ll1)))
			overlap_bin(i) = fix(junk1(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_ll2(overlap_pixels) * (2. * i / nsteps)
			caltot(i) = abs(total(tempcal - calflux_ll1(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		ll2_ll1_frac = (2. * junkindex(0) / nsteps)

		calflux_ll2 = flux_ll2 * ll2_ll1_frac
	endelse

; Anchor SL1 to LL2

	first_pixel = wave_ll2(0)
	last_pixel= wave_sl1(n_elements(wave_sl1)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_sl1(n_elements(flux_sl1)-npix:n_elements(flux_sl1)-1)
		first_sec = calflux_ll2(0:npix-1)

		sl1_ll2_frac = mean(first_sec) / mean(last_sec)	
		
		calflux_sl1 = flux_sl1 * sl1_ll2_frac
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_sl1 ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin
			junk2 = where(abs(wave_sl1(overlap_pixels(0)+i) - wave_ll2) eq $
				min(abs(wave_sl1(overlap_pixels(0)+i) - wave_ll2)))
			overlap_bin(i) = fix(junk2(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_sl1(overlap_pixels) * (4. * i / nsteps)
			caltot(i) = abs(total(tempcal - calflux_ll2(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		sl1_ll2_frac = (4. * junkindex(0) / nsteps)

		calflux_sl1 = flux_sl1 * sl1_ll2_frac
	endelse

; Anchor SL2 to SL1

	first_pixel = wave_sl1(0)
	last_pixel= wave_sl2(n_elements(wave_sl2)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_sl2(n_elements(flux_sl2)-npix:n_elements(flux_sl2)-1)
		first_sec = calflux_sl1(0:npix-1)

		sl2_sl1_frac = mean(first_sec) / mean(last_sec)	
		
		calflux_sl2 = flux_sl2 * sl2_sl1_frac
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_sl2 ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin 
			junk3 = where(abs(wave_sl2(overlap_pixels(0)+i) - wave_sl1) eq $
				min(abs(wave_sl2(overlap_pixels(0)+i) - wave_sl1)))
			overlap_bin(i) = fix(junk3(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_sl2(overlap_pixels) * (2. * i / nsteps)
			caltot(i) = abs(total(tempcal - calflux_sl1(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		sl2_sl1_frac = (2. * junkindex(0) / nsteps)

		calflux_sl2 = flux_sl2 * sl2_sl1_frac
	endelse

; Calibrate bonus orders to first-order modules

calflux_sl3 = flux_sl3 * sl1_ll2_frac
calflux_ll3 = flux_ll3


; Open the individual nods and calibrate them using the same scale factors as the coadded spectra

nodpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/nods/'

	; Lores

	readcol, nodpath+fname+'_sl1_1p.tbl', $
		det_sl1_1p, wave_sl1_1p, flux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl1_2p.tbl', $
		det_sl1_2p, wave_sl1_2p, flux_sl1_2p, err_sl1_2p, bit_sl1_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl2_1p.tbl', $
		det_sl2_1p, wave_sl2_1p, flux_sl2_1p, err_sl2_1p, bit_sl2_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl2_2p.tbl', $
		det_sl2_2p, wave_sl2_2p, flux_sl2_2p, err_sl2_2p, bit_sl2_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl3_1p.tbl', $
		det_sl3_1p, wave_sl3_1p, flux_sl3_1p, err_sl3_1p, bit_sl3_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl3_2p.tbl', $
		det_sl3_2p, wave_sl3_2p, flux_sl3_2p, err_sl3_2p, bit_sl3_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll1_1p.tbl', $
		det_ll1_1p, wave_ll1_1p, flux_ll1_1p, err_ll1_1p, bit_ll1_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll1_2p.tbl', $
		det_ll1_2p, wave_ll1_2p, flux_ll1_2p, err_ll1_2p, bit_ll1_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll2_1p.tbl', $
		det_ll2_1p, wave_ll2_1p, flux_ll2_1p, err_ll2_1p, bit_ll2_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll2_2p.tbl', $
		det_ll2_2p, wave_ll2_2p, flux_ll2_2p, err_ll2_2p, bit_ll2_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll3_1p.tbl', $
		det_ll3_1p, wave_ll3_1p, flux_ll3_1p, err_ll3_1p, bit_ll3_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll3_2p.tbl', $
		det_ll3_2p, wave_ll3_2p, flux_ll3_2p, err_ll3_2p, bit_ll3_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	; Hires

	readcol, nodpath+fname+'_sh_1p.tbl', $
		det_sh_1p, wave_sh_1p, flux_sh_1p, err_sh_1p, bit_sh_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sh_2p.tbl', $
		det_sh_2p, wave_sh_2p, flux_sh_2p, err_sh_2p, bit_sh_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_lh_1p.tbl', $
		det_lh_1p, wave_lh_1p, flux_lh_1p, err_lh_1p, bit_lh_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_lh_2p.tbl', $
		det_lh_2p, wave_lh_2p, flux_lh_2p, err_lh_2p, bit_lh_2p, format = 'i,f,f,f,i', skipline = 1, /silent

if not keyword_set(nocal) then begin

; Calibrate the nods to peakups

; Lo-res

	calflux_sl1_1p = flux_sl1_1p * sl1_ll2_frac
	calflux_sl2_1p = flux_sl2_1p * sl2_sl1_frac
	calflux_sl3_1p = flux_sl3_1p * sl1_ll2_frac
	calflux_ll1_1p = flux_ll1_1p
	calflux_ll3_1p = flux_ll3_1p 
	calflux_ll2_1p = flux_ll2_1p * ll2_ll1_frac
	
	calflux_sl1_2p = flux_sl1_2p * sl1_ll2_frac
	calflux_sl2_2p = flux_sl2_2p * sl2_sl1_frac
	calflux_sl3_2p = flux_ll1_2p * sl1_ll2_frac
	calflux_ll1_2p = flux_ll1_2p 
	calflux_ll3_2p = flux_ll3_2p 
	calflux_ll2_2p = flux_ll2_2p * ll2_ll1_frac
	
; Option for non-calibrated spectra

endif

if keyword_set(nocal) then begin

	calflux_sl1 = flux_sl1
	calflux_sl2 = flux_sl2
	calflux_sl3 = flux_sl3
	calflux_ll1 = flux_ll1
	calflux_ll2 = flux_ll2
	calflux_ll3 = flux_ll3
	
	calflux_sl1_1p = flux_sl1_1p
	calflux_sl2_1p = flux_sl2_1p
	calflux_sl3_1p = flux_sl3_1p
	calflux_ll1_1p = flux_ll1_1p
	calflux_ll2_1p = flux_ll2_1p
	calflux_ll3_1p = flux_ll3_1p
	
	calflux_sl1_2p = flux_sl1_2p
	calflux_sl2_2p = flux_sl2_2p
	calflux_sl3_2p = flux_sl3_2p
	calflux_ll1_2p = flux_ll1_2p
	calflux_ll2_2p = flux_ll2_2p
	calflux_ll3_2p = flux_ll3_2p

endif

; Default is to perform no calibration on high-res modules

calflux_sh = flux_sh
calflux_lh = flux_lh

calflux_sh_1p = flux_sh_1p
calflux_sh_2p = flux_sh_2p
calflux_lh_1p = flux_lh_1p
calflux_lh_2p = flux_lh_2p

; Plot results

!p.multi = [0,1,2]

; LORES spectra

pspath = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibration/'

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = pspath+fname+'_calflux.ps', /landscape, /color, bits_per_pixel = 8
	axiscolor = fsc_color('Black')
	bgcolor = fsc_color('White')
	cs = 1
	ct = 2
endif else begin
	axiscolor = fsc_color('White')
	bgcolor = fsc_color('Black')
	cs = 2
	ct = 1
endelse

plot, wave_sl1, flux_sl1, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
	xr = [4,40], /xstyle, $
	yr = [1d-4,1d1], /ystyle, $
	/xlog,/ylog, $
	title = obj, $
	charsize = cs, $
	color = axiscolor, $
	background = bgcolor, $
	linestyle = 2, $
	thick = ct, $
	charthick = ct, $
	/nodata


oplot, wave_sl1, flux_sl1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_sl2, flux_sl2, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll1, flux_ll1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll2, flux_ll2, color = axiscolor, linestyle = 2, thick = ct

oplot, wave_ll2, calflux_ll2, color = blue, thick = ct
oplot, wave_ll1, calflux_ll1, color = red, thick = ct
oplot, wave_sl2, calflux_sl2, color = yellow, thick = ct
oplot, wave_sl1, calflux_sl1, color = green, thick = ct

if not keyword_set(nobonus) then begin
	oplot, wave_sl3, flux_sl3, color = axiscolor, linestyle = 2, thick = ct
	oplot, wave_ll3, flux_ll3, color = axiscolor, linestyle = 2, thick = ct
	oplot, wave_ll3, calflux_ll3, color = orange, thick = ct
	oplot, wave_sl3, calflux_sl3, color = orange, thick = ct
endif

oplot, rwave, rtrans / max(rtrans) * max(flux_ll1), linestyle = 1, color=red, thick = ct
oplot, bwave, btrans / max(btrans) * max(flux_ll1), linestyle = 1, color=blue, thick = ct


xyouts, 0.2, 0.85, 'Lo-res', /normal, charsize = cs, charthick = ct, color = axiscolor

; HIRES spectra

plot, wave_sh, flux_sh, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
	xr = [5,40], $
	yr = [0,max(flux_lh)], $
	title = obj, $
	charsize = cs, $
	linestyle = 2, $
	color = axiscolor, $
	background = bgcolor, $
	thick = ct, $
	charthick = ct, $
	/nodata

;oplot, wave_sh, flux_sh, color = axiscolor, linestyle = 1, thick = ct, psym = 10
;oplot, wave_lh, flux_lh, color = axiscolor, linestyle = 1, thick = ct, psym = 10
oplot, wave_sh, flux_sh, color = blue, linestyle = 0, thick = ct, psym = 10
oplot, wave_lh, flux_lh, color = red, linestyle = 0, thick = ct, psym = 10

oplot, rwave, rtrans / max(rtrans) * max(flux_lh), linestyle = 1, color=red, thick = ct
oplot, bwave, btrans / max(btrans) * max(flux_lh), linestyle = 1, color=blue, thick = ct

xyouts, 0.2, 0.35, 'Hi-res', /normal, charsize = cs, charthick = ct, color = axiscolor
xyouts, 0.8, 0.95, fname, /normal, charsize = cs, charthick = ct, color=axiscolor

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif


if keyword_set(write) then begin

	; Writing stitched spectra to disk
	
	writepath_coadd = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/stitched/coadd/'
	writepath_nod   = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/stitched/nods/'

	; Coadded

	forprint, det_sl1, wave_sl1, calflux_sl1, err_sl1, bit_sl1, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl1_cal.tbl', /silent

	forprint, det_sl2, wave_sl2, calflux_sl2, err_sl2, bit_sl2, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl2_cal.tbl', /silent

	forprint, det_sl3, wave_sl3, calflux_sl3, err_sl3, bit_sl3, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl3_cal.tbl', /silent

	forprint, det_ll1, wave_ll1, calflux_ll1, err_ll1, bit_ll1, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll1_cal.tbl', /silent

	forprint, det_ll2, wave_ll2, calflux_ll2, err_ll2, bit_ll2, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll2_cal.tbl', /silent

	forprint, det_ll3, wave_ll3, calflux_ll3, err_ll3, bit_ll3, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll3_cal.tbl', /silent

	forprint, det_sh, wave_sh, calflux_sh, err_sh, bit_sh, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sh_cal.tbl', /silent

	forprint, det_lh, wave_lh, calflux_lh, err_lh, bit_lh, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_lh_cal.tbl', /silent

	; Nods

	forprint, det_sl1_1p, wave_sl1_1p, calflux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl1_1p_cal.tbl', /silent

	forprint, det_sl1_2p, wave_sl1_2p, calflux_sl1_2p, err_sl1_2p, bit_sl1_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl1_2p_cal.tbl', /silent

	forprint, det_sl2_1p, wave_sl2_1p, calflux_sl2_1p, err_sl2_1p, bit_sl2_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl2_1p_cal.tbl', /silent

	forprint, det_sl2_2p, wave_sl2_2p, calflux_sl2_2p, err_sl2_2p, bit_sl2_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl2_2p_cal.tbl', /silent

	forprint, det_sl3_1p, wave_sl3_1p, calflux_sl3_1p, err_sl3_1p, bit_sl3_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl3_1p_cal.tbl', /silent

	forprint, det_sl3_2p, wave_sl3_2p, calflux_sl3_2p, err_sl3_2p, bit_sl3_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl3_2p_cal.tbl', /silent

	forprint, det_ll1_1p, wave_ll1_1p, calflux_ll1_1p, err_ll1_1p, bit_ll1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll1_1p_cal.tbl', /silent

	forprint, det_ll1_2p, wave_ll1_2p, calflux_ll1_2p, err_ll1_2p, bit_ll1_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll1_2p_cal.tbl', /silent

	forprint, det_ll2_1p, wave_ll2_1p, calflux_ll2_1p, err_ll2_1p, bit_ll2_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll2_1p_cal.tbl', /silent

	forprint, det_ll2_2p, wave_ll2_2p, calflux_ll2_2p, err_ll2_2p, bit_ll2_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll2_2p_cal.tbl', /silent

	forprint, det_ll3_1p, wave_ll3_1p, calflux_ll3_1p, err_ll3_1p, bit_ll3_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll3_1p_cal.tbl', /silent

	forprint, det_ll3_2p, wave_ll3_2p, calflux_ll3_2p, err_ll3_2p, bit_ll3_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll3_2p_cal.tbl', /silent

	forprint, det_sh_1p, wave_sh_1p, calflux_sh_1p, err_sh_1p, bit_sh_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sh_1p_cal.tbl', /silent

	forprint, det_sh_2p, wave_sh_2p, calflux_sh_2p, err_sh_2p, bit_sh_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sh_2p_cal.tbl', /silent

	forprint, det_lh_1p, wave_lh_1p, calflux_lh_1p, err_lh_1p, bit_lh_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_lh_1p_cal.tbl', /silent

	forprint, det_lh_2p, wave_lh_2p, calflux_lh_2p, err_lh_2p, bit_lh_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_lh_2p_cal.tbl', /silent

endif

if not keyword_set(quiet) then begin
	print, 'LL2 to LL1 scale factor: ', ll2_ll1_frac
	print, 'SL1 to LL2 scale factor: ', sl1_ll2_frac
	print, 'SL2 to SL1 scale factor: ', sl2_sl1_frac
	;print, 'SH to LH scale factor: ', hr_frac
	print, 'Completed STITCH_ORDERS on ',fname
endif

if not keyword_set(nospicehdr) then spicehdr, fname, readdir = 'stitched'

; Save the scaling fractions for later use

save, ll2_ll1_frac, sl1_ll2_frac, sl2_sl1_frac, $ ;,hr_frac
	filename='~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/scaling_frac/'+fname+'.sav'

; End program

if keyword_set(stop) then stop
end
