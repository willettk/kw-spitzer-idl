pro fullql_hr, fname, xr = xr, yr = yr, nods = nods, nolines = nolines, write = write, ps = ps, $
	bw = bw, notrim = notrim, no_offset = no_offset, noclean = noclean, stop = stop, $
	shsky = shsky, lhsky = lhsky, restframe = restframe, panel = panel, quiet = quiet
;+
; NAME:
;       FULLQL_HR
;
; PURPOSE:
; 	Display trimmed hi-res Spitzer IRS spectra
;
; OUTPUTS:
;	- Plots individual spectra, overplotting both optimal and regular extractions
;
; KEYWORDS:
;
; 	NODS - displays both nod positions with nod 2 artificially offset upwards. Default is
;		to display the coadded nods. 
;
;	NOLINES - removes indicators of common emission and absorption features in ULIRGs. Default
;			is to display them.
;
;	WRITE - write the spectra w/trimming to an ASCII file
;
;	PS - create a hard copy of spectra plot
;
;	NOTRIM - display raw spectra without trimming bad pixels at the edges
;
;	NO_OFFSET - overplot nods in true units rather than offsetting them
;
;	NOCLEAN - view data that has NOT gone through IRSCLEAN_MASK. Default is to view cleaned data. 
;
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;
; EXAMPLE:
;	IDL> fullql_hr, 'mega023'
;
; REVISION HISTORY
;       Written by K. Willett                Apr 2007
;	Added extensive commenting, rewrote file paths to work on morbo - KW, May 07
;	Added WRITE, PS keywords - KW, May 07
; 	Lookup object characteristics from TARGETS.pro - KW, May 07
;	Added NOTRIM keyword - KW, May 07
;	Added NOCLEAN keyword (IRSCLEAN_MASK now part of default pipeline) - KW, Jul 07
;	Added lookup tables for trim templates at order edges - KW, Jul 07
;	Fixed bug in SMART tables (sorting by wavelength vs. order) - KW, Aug 07
;	Negative pixels in only one nod are not used in the co-added spectrum - KW, Jul 09
;-

; Set device so that colors plotted for different orders appear on the screen

if not keyword_set(panel) then device, decomposed = 1, window_state = state

; Locate directory from which to read (either megamasers or CSOs)

tag, fname, dirtag

; Locate files and read them into IDL

	writedir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/'
	plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/'
	opath = '~/spice/output/hires/optimal/all/'

if keyword_set(noclean) then files_opt = opath+'*'+fname+'*p_spect*' $
	else files_opt = opath+'*'+fname+'*clean*spect.tbl'

ofiles = findfile(files_opt)
onum = n_elements(ofiles)

; Use the spectroscopic redshifts of the objects (from NED) to plot images in the target rest frame

targets, fname, redshift, obj

if keyword_set(restframe) then begin
	zzz = redshift
	redshift = 0.
endif

	; Read in files for optimal extraction (SH nod 1, SH nod 2, LH nod 1, LH nod 2)
	; Data read: spectral order, wavelength, flux, error, and bit type

	readcol, ofiles(0), odet0, owave0, oflux0, oerr0, obit0, format = 'i,f,f,f,i', /silent
	readcol, ofiles(1), odet1, owave1, oflux1, oerr1, obit1, format = 'i,f,f,f,i', /silent
	readcol, ofiles(2), odet2, owave2, oflux2, oerr2, obit2, format = 'i,f,f,f,i', /silent
	readcol, ofiles(3), odet3, owave3, oflux3, oerr3, obit3, format = 'i,f,f,f,i', /silent

	; Remove known persistent bad pixel at 20.07722 um in order 20

	badind2 = where(owave2 eq 20.07722)
	badind3 = where(owave3 eq 20.07722)
	
	if odet2(badind2(0)) ne 20 then begin
		goodind2 = where(indgen(n_elements(owave2)) ne badind2(0))
		odet2 = odet2(goodind2) & owave2 = owave2(goodind2) & oflux2 = oflux2(goodind2) & oerr2 = oerr2(goodind2) & obit2 = obit2(goodind2)
	endif

	if odet3(badind3(0)) ne 20 then begin
		goodind3 = where(indgen(n_elements(owave3)) ne badind3(0))
		odet3 = odet3(goodind3) & owave3 = owave3(goodind3) & oflux3 = oflux3(goodind3) & oerr3 = oerr3(goodind3) & obit3 = obit3(goodind3)
	endif

	; Templates for trimming pixels from the edges of orders

	; SH

	if keyword_set(shsky) then $
		shtrim_file = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/trim_templates/hrsky/'+fname+'_sh.trim') else $
		shtrim_file = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/trim_templates/'+fname+'_sh.trim')

	if shtrim_file(0) eq '' then begin

		sh_20 = [10,10]
		sh_19 = [12,20]
		sh_18 = [10,22]
		sh_17 = [12,20]
		sh_16 = [10,25]
		sh_15 = [10,30]
		sh_14 = [10,17]
		sh_13 = [22,15]
		sh_12 = [22,10]
		sh_11 = [10,35]

		if not keyword_set(quiet) then print,'Using standard template for SH modules on '+fname
		sh_trim = [[sh_20], [sh_19], [sh_18], [sh_17], [sh_16], [sh_15], [sh_14], [sh_13], [sh_12], [sh_11]] 

	endif else begin

		if not keyword_set(quiet) then print,'Using saved template for SH modules on '+fname
		readcol, shtrim_file(0), trimstart, trimend, format = 'x,i,i', skipline = 2, /silent

		sh_trim = [transpose(trimstart),transpose(trimend)]
		
	endelse

	oind0 = intarr(1)
	oind1 = intarr(1)

	; LH

	if keyword_set(lhsky) then $
		lhtrim_file = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/trim_templates/hrsky/'+fname+'_lh.trim') else $
		lhtrim_file = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/trim_templates/'+fname+'_lh.trim')

	if lhtrim_file(0) eq '' then begin

		lh_20 = [27,20]
		lh_19 = [10,26]
		lh_18 = [10,25]
		lh_17 = [10,25]
		lh_16 = [10,15]
		lh_15 = [22,10]
		lh_14 = [10,10]
		lh_13 = [10,22]
		lh_12 = [10,20]
		lh_11 = [0,42]

		if not keyword_set(quiet) then print,'Using standard template for LH modules on '+fname
		lh_trim = [[lh_20], [lh_19], [lh_18], [lh_17], [lh_16], [lh_15], [lh_14], [lh_13], [lh_12], [lh_11]] 

	endif else begin

		if not keyword_set(quiet) then print,'Using saved template for LH modules on '+fname
		readcol, lhtrim_file(0), trimstart, trimend, format = 'x,i,i', skipline = 2, /silent

		lh_trim = [transpose(trimstart),transpose(trimend)]
		
	endelse

	oind2 = intarr(1)
	oind3 = intarr(1)


	; SH trimming

	if keyword_set(notrim) then sh_trim = intarr(2,10)

	for i = 0,9 do begin
		junk0 = where(odet0 eq 20 - i)
		sortedind0 = junk0(sort(owave0(junk0)))
		junk0 = sortedind0(sh_trim(0,i):n_elements(sortedind0) - sh_trim(1,i) - 1)
		oind0 = [oind0,junk0]
		junk1 = where(odet1 eq 20 - i)
		sortedind1 = junk1(sort(owave1(junk1)))
		junk1 = sortedind1(sh_trim(0,i):n_elements(sortedind1) - sh_trim(1,i) - 1)
		oind1 = [oind1,junk1]
	endfor

	oind0 = oind0(1:n_elements(oind0)-1)
	owave0 = owave0(oind0) & oflux0 = oflux0(oind0) & odet0 = odet0(oind0) & oerr0 = oerr0(oind0) & obit0 = obit0(oind0)
	oind1 = oind1(1:n_elements(oind1)-1)
	owave1 = owave1(oind1) & oflux1 = oflux1(oind1) & odet1 = odet1(oind1) & oerr1 = oerr1(oind1) & obit1 = obit1(oind1)

	; LH trimming

	if keyword_set(notrim) then lh_trim = intarr(2,10)

	for i = 0,9 do begin
		junk2 = where(odet2 eq 20 - i)
		sortedind2 = junk2(sort(owave2(junk2)))
		junk2 = sortedind2(lh_trim(0,i):n_elements(sortedind2) - lh_trim(1,i) - 1)
;		junk2 = junk2(lh_trim(0,i):n_elements(junk2) - lh_trim(1,i) - 1)
		oind2 = [oind2,junk2]
		junk3 = where(odet3 eq 20 - i)
		sortedind3 = junk3(sort(owave3(junk3)))
		junk3 = sortedind3(lh_trim(0,i):n_elements(sortedind3) - lh_trim(1,i) - 1)
;		junk3 = junk3(lh_trim(0,i):n_elements(junk3) - lh_trim(1,i) - 1)
		oind3 = [oind3,junk3]
	endfor

	oind2 = oind2(1:n_elements(oind2)-1)
	owave2 = owave2(oind2) & oflux2 = oflux2(oind2) & odet2 = odet2(oind2) & oerr2 = oerr2(oind2) & obit2 = obit2(oind2)
	oind3 = oind3(1:n_elements(oind3)-1)
	owave3 = owave3(oind3) & oflux3 = oflux3(oind3) & odet3 = odet3(oind3) & oerr3 = oerr3(oind3) & obit3 = obit3(oind3)

; Coadd the nod positions using a weighted mean based on the spectral uncertainty

weights_sh_1p = 1d-6/oerr0^2
weights_sh_2p = 1d-6/oerr1^2

weights_lh_1p = 1d-6/oerr2^2
weights_lh_2p = 1d-6/oerr3^2

sh_coadd = (weights_sh_1p * oflux0 + weights_sh_2p * oflux1) / (weights_sh_1p + weights_sh_2p)
lh_coadd = (weights_lh_1p * oflux2 + weights_lh_2p * oflux3) / (weights_lh_1p + weights_lh_2p)

sh_coadd_err = sqrt(oerr0^2 + oerr1^2)
lh_coadd_err = sqrt(oerr2^2 + oerr3^2)

	; If pixel flux for one nod is negative, use only positive flux from the other nod

	posflux_sh_1p = where(oflux0 gt 0. and oflux1 lt 0.,sh_1p_poscount)
	posflux_sh_2p = where(oflux0 lt 0. and oflux1 gt 0.,sh_2p_poscount)
	if sh_1p_poscount gt 0 then begin
		sh_coadd[posflux_sh_1p] = oflux0[posflux_sh_1p]
		sh_coadd_err[posflux_sh_1p] = oerr0[posflux_sh_1p]
	endif
	if sh_2p_poscount gt 0 then begin
		sh_coadd[posflux_sh_2p] = oflux1[posflux_sh_2p]
		sh_coadd_err[posflux_sh_2p] = oerr1[posflux_sh_2p]
	endif
	
	posflux_lh_1p = where(oflux2 gt 0. and oflux3 lt 0.,lh_1p_poscount)
	posflux_lh_2p = where(oflux2 lt 0. and oflux3 gt 0.,lh_2p_poscount)
	if lh_1p_poscount gt 0 then begin
		lh_coadd[posflux_lh_1p] = oflux2[posflux_lh_1p]
		lh_coadd_err[posflux_lh_1p] = oerr2[posflux_lh_1p]
	endif
	if lh_2p_poscount gt 0 then begin
		lh_coadd[posflux_lh_2p] = oflux3[posflux_lh_2p]
		lh_coadd_err[posflux_lh_2p] = oerr3[posflux_lh_2p]
	endif


allflux_coadd = [sh_coadd, lh_coadd]
allwave_coadd = [owave0, owave2]
allerr_coadd = [sh_coadd_err, lh_coadd_err]
alldet_coadd = [odet0, odet2]
allbit_coadd = [obit0, obit2]

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

lines = [8.025, 8.6, 8.99138, 9.665, 10.511, 11.3, 12.279, 12.6, 12.814, 13.7, $
	14.0, 14.2, 14.322, 14.368, 15.0, 15.555, 16.4, 17.035, 17.4, 17.934, 18.713, 24.318, 25.890, $
	25.988, 28.218, 33.481, 34.815]

line_id = ['H!I2!N S(4)', 'PAH', 'ArIII', 'H!I2!N S(3)', 'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
	'C!I2!NH!I2!N', 'HCN', 'PAH', 'NeV', 'ClII', 'CO!I2!N', 'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', 'FeII', 'SIII', $
	'NeV', 'OIV', 'FeII', 'H!I2!N S(0)', 'SIII', 'SiII']

templines = ir_lines(/hr)
lines = templines(*,0) & line_id = templines(*,1)
if keyword_set(restframe) then lines = lines * (1. + zzz)

; Create indices for spectra arranged by increasing wavelength (removing any jumps due to overlapping spectra in different orders)

b0 = sort(owave0)
b1 = sort(owave1)
b2 = sort(owave2)
b3 = sort(owave3)

allwave_1p = [owave0(b0),owave2(b2)]
allflux_1p = [oflux0(b0),oflux2(b2)]
allerr_1p = [oerr0(b0),oerr2(b2)]
alldet_1p = [odet0(b0),odet2(b2)]
allbit_1p = [obit0(b0),obit2(b2)]

allwave_2p = [owave1(b1),owave3(b3)]
allflux_2p = [oflux1(b1),oflux3(b3)]
allerr_2p = [oerr1(b1),oerr3(b3)]
alldet_2p = [odet1(b1),odet3(b3)]
allbit_2p = [obit1(b1),obit3(b3)]

; Set plot ranges (if not done in command line)

if not keyword_set(xr) then xr = [8,18]
if not keyword_set(yr) then yr = [-0.1,1.]

;if state(0) eq 0 then window,0 else wset,0

; Hard copy option

if keyword_set(ps) then begin
	set_plot,'ps'
	device, /landscape, $
		filename = plotdir+strjoin(strsplit(obj,' ',/extract))+'_spect_hr.ps'
	defcolor = fsc_color('Black')
	lthick = 2
	cthick = 2
endif else begin
	defcolor = fsc_color('White')
	lthick = 1
	cthick = 1
endelse

; Sort the data by order for color coding (same as SMART)

colors = [fsc_color('Cyan'), $
	fsc_color('Magenta'), $
	defcolor, $
	fsc_color('Red'), $
	fsc_color('Green'), $
	fsc_color('Blue'), $
	fsc_color('Yellow'), $
	fsc_color('Cyan'), $
	fsc_color('Magenta'), $
	fsc_color('Goldenrod')]

if keyword_set(ps) then colors = replicate(defcolor,10)
if keyword_set(bw) then colors = replicate(defcolor,10)

; Plot the data in rest frame wavelengths 

if not keyword_set(panel) then !p.multi = [0,1,1]

plot, owave0(b0)/ (redshift + 1), oflux0(b0), $
	xrange = xr, $
	/xstyle, $
	yrange = yr, $
	/ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = obj, $
	charsize = 2, $
	color = defcolor, $
	thick = lthick, $
	xthick = lthick, $
	ythick = lthick, $
	charthick = cthick, $
	/nodata

; Plot coadded data

if not keyword_set(nods) then begin

	for i = 0, 9 do begin

		ajunk0 = where(odet0 eq (11 + i))
		cjunk0 = sort(owave0(ajunk0))

		ajunk2 = where(odet2 eq (11 + i))
		cjunk2 = sort(owave2(ajunk2))

		oplot, owave0(ajunk0(cjunk0))/(redshift + 1.), sh_coadd(ajunk0(cjunk0)), $
			color = colors(i), $
			psym = 10, thick = lthick

		oplot, owave2(ajunk2(cjunk2))/(redshift + 1.), lh_coadd(ajunk2(cjunk2)), $
			color = colors(i), $
			psym = 10, thick = lthick

	
	endfor
endif else begin
	
; Plot individual nod positions

	for i = 0, 9 do begin

		ajunk0 = where(odet0 eq (11 + i))
		cjunk0 = sort(owave0(ajunk0))
		ajunk1 = where(odet1 eq (11 + i))
		cjunk1 = sort(owave1(ajunk1))
		ajunk2 = where(odet2 eq (11 + i))
		cjunk2 = sort(owave2(ajunk2))
		ajunk3 = where(odet3 eq (11 + i))
		cjunk3 = sort(owave3(ajunk3))

		; Plot data for the 1st nod position

		oplot, owave0(ajunk0(cjunk0))/(redshift + 1.), oflux0(ajunk0(cjunk0)), $
			color = colors(i), $
			psym = 10, thick = lthick

		oplot, owave2(ajunk2(cjunk2))/(redshift + 1.), oflux2(ajunk2(cjunk2)), $
			color = colors(i), $
			psym = 10, thick = lthick

		; Plot data for the 2nd nod position (may be offset from the first)
	
		if keyword_set(no_offset) then offset = 0 else offset = 0.2

		oplot, owave1(ajunk1(cjunk1))/(redshift + 1.), oflux1(ajunk1(cjunk1)) + offset, $
			color = colors(i), $
			psym = 10, thick = lthick

		oplot, owave3(ajunk3(cjunk3))/(redshift + 1.), oflux3(ajunk3(cjunk3)) + offset, $
			color = colors(i), $
			psym = 10, thick = lthick
	endfor
endelse

; Overplot vertical lines with labels indicating locations of likely absorption and emission features

if not keyword_set(nolines) then begin
	for i = 0, n_elements(lines) - 1 do begin
		off = i mod 2
		ver, lines(i), linestyle = 1, color = defcolor
		xyouts, lines(i), 0.85*yr(1) + 0.05*off, line_id(i), orientation = 90, charsize = 1.5, /data, $
			color = defcolor, charthick = cthick
	endfor
endif

; Hard copy

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

; Write the new, trimmed spectra to file

if keyword_set(write) then begin

	; All modules

	forprint, alldet_1p, allwave_1p, allflux_1p, allerr_1p, allbit_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_hr_1p.tbl', /silent

	forprint, alldet_2p, allwave_2p, allflux_2p, allerr_2p, allbit_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_hr_2p.tbl', /silent

	forprint, alldet_coadd, allwave_coadd, allflux_coadd, allerr_coadd, allbit_coadd, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_hr_coadd.tbl', /silent

	; Individual modules

	forprint, odet0(b0), owave0(b0), oflux0(b0), oerr0(b0), obit0(b0), format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sh_1p.tbl', /silent

	forprint, odet2(b2), owave2(b2), oflux2(b2), oerr2(b2), obit2(b2), format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_lh_1p.tbl', /silent

	forprint, odet1(b1), owave1(b1), oflux1(b1), oerr1(b1), obit1(b1), format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sh_2p.tbl', /silent

	forprint, odet3(b3), owave3(b3), oflux3(b3), oerr3(b3), obit3(b3), format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_lh_2p.tbl', /silent

	forprint, odet0(b0), owave0(b0), sh_coadd(b0), oerr0(b0), obit0(b0), format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_sh_coadd.tbl', /silent

	forprint, odet2(b2), owave2(b2), lh_coadd(b2), oerr2(b2), obit2(b2), format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_lh_coadd.tbl', /silent

endif

if keyword_set(stop) then stop
end
