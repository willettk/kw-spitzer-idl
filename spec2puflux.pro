pro spec2puflux, fname, ps = ps, nostop = nostop, write = write, nocal = nocal
;+
; NAME:
;       SPEC2PUFLUX
;
; PURPOSE:
; 	Read in a trimmed, extracted spectrum along with the measured peakup flux
;	and calibrate spectrum to match PU
;
; INPUTS:
;
;	artificial_pu.csv - 	Excel spreadsheet containing measured artificial peakup
;				fluxes and true values for entire sample. 
;
; KEYWORDS:
;
;	PS - hard copies of overlaid raw and calibrated spectra
;
;	NOCAL - runs routine without calibration to add suffixes and read artificial fluxes in SMART
;
; EXAMPLE:
;	IDL> spec2puflux, 'mega001', /ps
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;	READCOL,pro
;
; NOTES:
;
;	To do: revise this so that the same stitching algorithm is applied as for the archival data
;	in SPEC2PUFLUX_ARCH. For the first run, this includes checking a given module for correlation
;	with the artificial PU flux and then stitching the rest; should do a secondary test at the end
;	to make sure that results are still kosher. 
;
;	Also want to store the peakup fluxes in an IDL structure, if possible, so that I don't have to 
;	keep reading spreadsheets. (?)
;
; REVISION HISTORY
;       Written by K. Willett                Jul 2007
; 	Major revisions - KW, Aug 2007
;	Added bonus order - KW, Aug 07
;-

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


if not keyword_set(nocal) then begin

; Read in the measured and artificial peakup fluxes

apupath = '~/Astronomy/Research/Spitzer/'

readcol, apupath+'artificial_pu.csv', atag, aobj, apu_blue_lo, apu_red_lo, apu_blue_hi, apu_red_hi, my_bluepu, my_redpu, $
	format = 'a,a,a,a,a,a,a,a', skipline = 3, /silent

; Restore the Spitzer PU fluxes from the IDL .sav file - CHANGE THIS - only applies to archived data!

;restore, '~/Astronomy/Research/Spitzer/archived/data/idl_sav/pu_fluxes.sav'
;if n_elements(archfiles) ne n_elements(atag) then begin
;	print,'Artificial PU list does not match IRS peakups'
;	stop
;endif

; Change units from mJy to Jy

tt = atag(0)
if strmid(tt,0,1) eq '"' then begin
	for i = 0, n_elements(atag) - 1 do begin
		junk1 = strmid(atag(i),1,strlen(atag(i)-2))
		atag(i) = junk1
		junk2 = strmid(aobj(i),1,strlen(aobj(i)-2))
		aobj(i) = junk2
	endfor
endif


ind1 = where(apu_blue_lo ne 'nodata')
ind2 = where(apu_red_lo ne 'nodata')
ind3 = where(apu_blue_hi ne 'nodata')
ind4 = where(apu_blue_hi ne 'nodata')
ind5 = where(my_bluepu ne 'nodata')
ind6 = where(my_redpu ne 'nodata')

tmp1 = fltarr(n_elements(apu_blue_lo))
tmp2 = fltarr(n_elements(apu_blue_lo))
tmp3 = fltarr(n_elements(apu_blue_lo))
tmp4 = fltarr(n_elements(apu_blue_lo))
tmp5 = fltarr(n_elements(apu_blue_lo))
tmp6 = fltarr(n_elements(apu_blue_lo))

tmp1(ind1) = float(apu_blue_lo(ind1)) * 1d-3
tmp2(ind2)  = float(apu_red_lo(ind2)) * 1d-3
tmp3(ind3) = float(apu_blue_hi(ind3)) * 1d-3
tmp4(ind4)  = float(apu_red_hi(ind4)) * 1d-3
tmp5(ind5)   = float(my_bluepu(ind5)) * 1d-3
tmp6(ind6)    = float(my_redpu(ind6)) * 1d-3

apu_blue_lo = tmp1 
apu_red_lo  = tmp2
apu_blue_hi = tmp3 
apu_red_hi  = tmp4
my_bluepu   = tmp5
my_redpu    = tmp6

; Locate spectrum of choice

tempind = where(fname eq atag)
if (tempind lt 0 or n_elements(tempind) ne 1) then print, 'Did not find fluxes'

apu_blue_lo = apu_blue_lo(tempind(0))
apu_red_lo = apu_red_lo(tempind(0))
apu_blue_hi = apu_blue_hi(tempind(0))
apu_red_hi = apu_red_hi(tempind(0))
my_bluepu = my_bluepu(tempind(0))
my_redpu = my_redpu(tempind(0))

if my_bluepu eq 0 then my_bluepu = apu_blue_lo
if my_redpu eq 0 then my_redpu = apu_red_lo

; Calibrated spectra

calflux_ll1 = flux_ll1 * my_redpu / apu_red_lo
calflux_ll2 = flux_ll2 * my_bluepu / apu_blue_lo
calflux_ll3 = flux_ll3 * my_redpu / apu_red_lo

calflux_sh = flux_sh
calflux_lh = flux_lh

; Calibrate SL1, SL2 to LL2

; Solution:
; - find the region of pixels that overlap in both orders
; - find each bin in one that corresponds most closely to another
; - match the bins
; - step through different scale heights until residuals are minimized

; - If no overlap exists between the orders, simply match the mean/median of some number of border pixels.

; Check for overlap

lastsl1_pixel = wave_sl1(n_elements(wave_sl1)-1)
firstll2_pixel= wave_ll2(0)

if firstll2_pixel gt lastsl1_pixel then begin

	npix = 2
	last_sl1 = flux_sl1(n_elements(flux_sl1)-npix:n_elements(flux_sl1)-1)
	first_ll2 = calflux_ll2(0:npix-1)
	ll2_sl1_frac = mean(first_ll2) / mean(last_sl1)
	
	calflux_sl1 = ll2_sl1_frac * flux_sl1
	calflux_sl2 = ll2_sl1_frac * flux_sl2
	calflux_sl3 = ll2_sl1_frac * flux_sl3		; Same fraction for bonus order calibration

endif else begin

	overlap_pixels = where(wave_sl1 ge firstll2_pixel)
	n_overlap = n_elements(overlap_pixels)
	overlap_bin = intarr(n_overlap)
	for i = 0, n_overlap - 1 do $
		overlap_bin(i) = where(abs(wave_sl1(overlap_pixels(0)+i) - wave_ll2) eq min(abs(wave_sl1(overlap_pixels(0)+i) - wave_ll2)))
	nsteps = 1000
	caltot = fltarr(nsteps)
	for i = 0, nsteps - 1 do begin
		tempcal = flux_sl1(overlap_pixels) * (0.5d + (1. * i / nsteps))
		caltot(i) = abs(total(tempcal - calflux_ll2(overlap_bin)))
	endfor

	junkindex = where(caltot eq min(caltot))
	ll2_sl1_frac = (0.5d + (1. * junkindex(0) / nsteps))

	calflux_sl1 = ll2_sl1_frac * flux_sl1
	calflux_sl2 = ll2_sl1_frac * flux_sl2
	calflux_sl3 = ll2_sl1_frac * flux_sl3

endelse

; Calibrate the hires to the lo-res spectra directly, instead of 
; using the artificial flux from SMART

; First, sample both the low and hires spectra at 0.5um intervals
; from 22 to 35 um

waveind_hi = intarr(2*(35-22)+1)
waveind_lo = intarr(2*(35-22)+1)

for i = 0, n_elements(waveind_hi) - 1 do begin
	waveind_hi(i) = where(abs(wave_lh - (i/2. +22)) eq min(abs(wave_lh - (i/2. + 22))))
	waveind_lo(i) = where(abs(wave_ll1 - (i/2. + 22)) eq min(abs(wave_ll1 - (i/2. + 22))))
endfor

; Step through a range of scale factors and find the SF
; that gives the minimum total residue for the flux. 

nsteps = 1000
caltot = fltarr(nsteps)
for i = 0, nsteps-1 do begin
	tempcal = flux_lh(waveind_hi) * (0.5d + (1. * i / nsteps))
	caltot(i) = abs(total(tempcal - calflux_ll1(waveind_lo)))
endfor

calmin_long = where(caltot eq min(caltot))
scalemin_long = (0.5d + (1. * calmin_long(0) / nsteps))

; Same process for lores, sampling from 15 to 19 um

waveind_hi = intarr(2*(19-15)+1)
waveind_lo = intarr(2*(19-15)+1)

for i = 0, n_elements(waveind_hi) - 1 do begin
	waveind_hi(i) = where(abs(wave_sh - (i/2. + 15)) eq min(abs(wave_sh - (i/2. + 15))))
	waveind_lo(i) = where(abs(wave_ll2 - (i/2. + 15)) eq min(abs(wave_ll2 - (i/2. + 15))))
endfor

; Step through a range of scale factors and find the SF
; that gives the minimum total residue for the flux. 

nsteps = 1000
caltot = fltarr(nsteps)
for i = 0, nsteps-1 do begin
	tempcal = flux_sh(waveind_hi) * (0.5d + (1. * i / nsteps))
	caltot(i) = abs(total(tempcal - calflux_ll2(waveind_lo)))
endfor

calmin_short = where(caltot eq min(caltot))
scalemin_short = (0.5d + (1. * calmin_short(0) / nsteps))

; Compare scale factors

print,''
print, 'Running ',fname
;print, 'Scaling to low-res: ', scalemin_long
;print, 'Scaling to red PU:  ', my_redpu / apu_red_hi

; Match the SH to the LH calibrated spectra

last_sh = flux_sh(n_elements(flux_sh)-11:n_elements(flux_sh)-1)
first_lh = calflux_lh(0:10)
lh_sh_frac_pu = median(first_lh) / median(last_sh)

last_sh = flux_sh(n_elements(flux_sh)-11:n_elements(flux_sh)-1)
first_lh = flux_lh(0:10) * scalemin_long
lh_sh_frac_lr = median(first_lh) / median(last_sh)

calflux_sh_alt_pu = lh_sh_frac_pu * flux_sh
calflux_sh_alt_lr = lh_sh_frac_lr * flux_sh

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
oplot, wave_sl3, flux_sl3, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll1, flux_ll1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll2, flux_ll2, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll3, flux_ll3, color = axiscolor, linestyle = 2, thick = ct

oplot, wave_ll2, calflux_ll2, color = blue, thick = ct
oplot, wave_ll3, calflux_ll3, color = orange, thick = ct
oplot, wave_ll1, calflux_ll1, color = red, thick = ct
oplot, wave_sl2, calflux_sl2, color = yellow, thick = ct
oplot, wave_sl3, calflux_sl3, color = orange, thick = ct
oplot, wave_sl1, calflux_sl1, color = green, thick = ct

oplot, rwave, rtrans / max(rtrans) * max(flux_ll1), linestyle = 1, color=red, thick = ct
oplot, bwave, btrans / max(btrans) * max(flux_ll1), linestyle = 1, color=blue, thick = ct


oplot, [16], [apu_blue_lo],  psym = 4, color = axiscolor,  symsize = 2
oplot, [22], [apu_red_lo],   psym = 4, color = axiscolor,  symsize = 2
oplot, [16], [my_bluepu], psym = 4, color = blue, symsize = 2
oplot, [22], [my_redpu],  psym = 4, color = red,  symsize = 2

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

oplot, wave_sh, flux_sh, color = axiscolor, linestyle = 0, thick = ct, psym = 10
oplot, wave_lh, flux_lh, color = axiscolor, linestyle = 0, thick = ct, psym = 10

;oplot, wave_sh, calflux_sh, color = blue, thick = ct
;oplot, wave_sh, flux_sh * scalemin_short, color = green, thick = ct
;oplot, wave_sh, calflux_sh_alt_pu, color = yellow, thick = ct
;oplot, wave_sh, calflux_sh_alt_lr, color = yellow, thick = ct
;oplot, wave_lh, calflux_lh, color = red, thick = ct
;oplot, wave_lh, flux_lh * scalemin_long, color = green, thick = ct
;oplot, wave_ll2, calflux_ll2, color = fsc_color("Hot Pink"), thick = ct
;oplot, wave_ll1, calflux_ll1, color = fsc_color("Hot Pink"), thick = ct

oplot, rwave, rtrans / max(rtrans) * max(flux_lh), linestyle = 1, color=red, thick = ct
oplot, bwave, btrans / max(btrans) * max(flux_lh), linestyle = 1, color=blue, thick = ct


;oplot, [16], [apu_blue_hi],  psym = 4, color = axiscolor,  symsize = 2
;oplot, [22], [apu_red_hi],   psym = 4, color = axiscolor,  symsize = 2
oplot, [16], [my_bluepu], psym = 4, color = blue, symsize = 2
oplot, [22], [my_redpu],  psym = 4, color = red,  symsize = 2

xyouts, 0.2, 0.35, 'Hi-res', /normal, charsize = cs, charthick = ct, color = axiscolor

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

endif		; NOCAL

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

; Calibrate to peakups

; Lo-res

calflux_sl1_1p = ll2_sl1_frac * flux_sl1_1p
calflux_sl2_1p = ll2_sl1_frac * flux_sl2_1p
calflux_sl3_1p = ll2_sl1_frac * flux_sl3_1p
calflux_ll1_1p = flux_ll1_1p * my_redpu / apu_red_lo
calflux_ll2_1p = flux_ll2_1p * my_bluepu / apu_blue_lo
calflux_ll3_1p = flux_ll3_1p * my_redpu / apu_red_lo

calflux_sl1_2p = ll2_sl1_frac * flux_sl1_2p
calflux_sl2_2p = ll2_sl1_frac * flux_sl2_2p
calflux_sl3_2p = ll2_sl1_frac * flux_sl3_2p
calflux_ll1_2p = flux_ll1_2p * my_redpu / apu_red_lo
calflux_ll2_2p = flux_ll2_2p * my_bluepu / apu_blue_lo
calflux_ll3_2p = flux_ll3_2p * my_redpu / apu_red_lo

; Hi-res

calflux_sh_1p = flux_sh_1p
calflux_lh_1p = flux_lh_1p

calflux_sh_alt_pu_1p = lh_sh_frac_pu * flux_sh_1p
calflux_sh_alt_lr_1p = lh_sh_frac_lr * flux_sh_1p
calflux_lh_alt_lt_1p = flux_lh_1p * scalemin_long

calflux_sh_2p = flux_sh_2p
calflux_lh_2p = flux_lh_2p

calflux_sh_alt_pu_2p = lh_sh_frac_pu * flux_sh_2p
calflux_sh_alt_lr_2p = lh_sh_frac_lr * flux_sh_2p
calflux_lh_alt_lt_2p = flux_lh_2p * scalemin_long

; Option for non-calibrated spectra (or when PUs are suspect)

endif

if keyword_set(nocal) then begin

	calflux_sl1 = flux_sl1
	calflux_sl2 = flux_sl2
	calflux_sl3 = flux_sl3
	calflux_ll1 = flux_ll1
	calflux_ll2 = flux_ll2
	calflux_ll3 = flux_ll3
	calflux_sh  = flux_sh
	calflux_lh  = flux_lh
	
	calflux_sl1_1p = flux_sl1_1p
	calflux_sl2_1p = flux_sl2_1p
	calflux_sl3_1p = flux_sl3_1p
	calflux_ll1_1p = flux_ll1_1p
	calflux_ll2_1p = flux_ll2_1p
	calflux_ll3_1p = flux_ll3_1p
	calflux_sh_1p  = flux_sh_1p
	calflux_lh_1p  = flux_lh_1p
	
	calflux_sl1_2p = flux_sl1_2p
	calflux_sl2_2p = flux_sl2_2p
	calflux_sl3_2p = flux_sl3_2p
	calflux_ll1_2p = flux_ll1_2p
	calflux_ll2_2p = flux_ll2_2p
	calflux_ll3_2p = flux_ll3_2p
	calflux_sh_2p  = flux_sh_2p
	calflux_lh_2p  = flux_lh_2p

endif

if keyword_set(write) then begin

	; Writing calibrated spectra to disk
	
	writepath_coadd = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
	writepath_nod   = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'

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

	forprint, det_sl1_1p, wave_sl1_1p, calflux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
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

print, 'Calibrated spectra to peakup fluxes for ',fname

if not keyword_set(nostop) then stop
end
