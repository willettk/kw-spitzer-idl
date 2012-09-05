pro lrcal, fname, bw = bw, nods=nods, ps = ps, nolines = nolines, $
	xr = xr, yr = yr
;+
; NAME:
;       LRCAL
;
; PURPOSE:
; 	Display trimmed, calibrated lo-res Spitzer IRS spectra
;
; INPUTS:
;	XR - plot range in x-direction (log scale)
;
;	YR - plot range in y-direction (log scale)
;
; OUTPUTS:
; 	- Plots low-res spectra sorted by order in rest frame
;
; KEYWORDS:
;
;	BW - displays spectra in black and white. Default is to display each order in a different color.
;
;	NOLINES - removes common mid-IR lines with labels. Default is to display them. 
;
;	NODS - displays individual nods vertically offset from each other (nod 2 on top). Default is
;		to display coadded spectra weighted by flux uncertainty in each point.
;
;	PS - plots the spectrum to a hardcopy postscript file
;
; EXAMPLE:
;	IDL> lrcal, 'mega023
;
; REQUIRES:
;
;	PAHFIT.pro (and associated routines)
;	READCOL.pro
;	TAG.pro
;	TARGETS.pro
;
; NOTES:
;
; REVISION HISTORY
;       Written by K. Willett                Sep 2007
;-

; Set device to read in colors

;if not keyword_set(ps) then device, decomposed = 1, window_state = state 
;device, window_state = state

; Locate directory from which to read (either megamasers or CSOs)

tag, fname, dirtag 
targets, fname, redshift, obj

; Read in the optimally extracted data from SPICE directory

plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibrated/'
copath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
nodpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'

co_flist = copath+'*'+fname+'*.tbl'
nod_flist = nodpath+'*'+fname+'*.tbl'

cofiles = findfile(co_flist)
conum = n_elements(cofiles)
nodfiles = findfile(nod_flist)
nodnum = n_elements(nodfiles)

if dirtag eq 'control' then begin
	if fname ne 'control004' or fname ne 'control008' or fname ne 'control013' then begin

		readcol, nodfiles(6), sl1_1p_order, sl1_1p_wave, sl1_1p_flux, sl1_1p_error, sl1_1p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(7), sl1_2p_order, sl1_2p_wave, sl1_2p_flux, sl1_2p_error, sl1_2p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(8), sl2_1p_order, sl2_1p_wave, sl2_1p_flux, sl2_1p_error, sl2_1p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(9), sl2_2p_order, sl2_2p_wave, sl2_2p_flux, sl2_2p_error, sl2_2p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(10), sl3_1p_order, sl3_1p_wave, sl3_1p_flux, sl3_1p_error, sl3_1p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(11), sl3_2p_order, sl3_2p_wave, sl3_2p_flux, sl3_2p_error, sl3_2p_bit, format = 'i,f,f,f,i', /silent
		
		readcol, nodfiles(0), ll1_1p_order, ll1_1p_wave, ll1_1p_flux, ll1_1p_error, ll1_1p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(1), ll1_2p_order, ll1_2p_wave, ll1_2p_flux, ll1_2p_error, ll1_2p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(2), ll2_1p_order, ll2_1p_wave, ll2_1p_flux, ll2_1p_error, ll2_1p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(3), ll2_2p_order, ll2_2p_wave, ll2_2p_flux, ll2_2p_error, ll2_2p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(4), ll3_1p_order, ll3_1p_wave, ll3_1p_flux, ll3_1p_error, ll3_1p_bit, format = 'i,f,f,f,i', /silent
		readcol, nodfiles(5), ll3_2p_order, ll3_2p_wave, ll3_2p_flux, ll3_2p_error, ll3_2p_bit, format = 'i,f,f,f,i', /silent
		
		readcol, cofiles(0), ll1_order, ll1_wave, ll1_flux, ll1_error, ll1_bit, format = 'i,f,f,f,i', /silent
		readcol, cofiles(1), ll2_order, ll2_wave, ll2_flux, ll2_error, ll2_bit, format = 'i,f,f,f,i', /silent
		readcol, cofiles(2), ll3_order, ll3_wave, ll3_flux, ll3_error, ll3_bit, format = 'i,f,f,f,i', /silent

		readcol, cofiles(3), sl1_order, sl1_wave, sl1_flux, sl1_error, sl1_bit, format = 'i,f,f,f,i', /silent
		readcol, cofiles(4), sl2_order, sl2_wave, sl2_flux, sl2_error, sl2_bit, format = 'i,f,f,f,i', /silent
		readcol, cofiles(5), sl3_order, sl3_wave, sl3_flux, sl3_error, sl3_bit, format = 'i,f,f,f,i', /silent
	endif
endif else begin

	; Read information in from tables for both nod positions
	
	readcol, nodfiles(10), sl1_1p_order, sl1_1p_wave, sl1_1p_flux, sl1_1p_error, sl1_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(11), sl1_2p_order, sl1_2p_wave, sl1_2p_flux, sl1_2p_error, sl1_2p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(12), sl2_1p_order, sl2_1p_wave, sl2_1p_flux, sl2_1p_error, sl2_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(13), sl2_2p_order, sl2_2p_wave, sl2_2p_flux, sl2_2p_error, sl2_2p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(14), sl3_1p_order, sl3_1p_wave, sl3_1p_flux, sl3_1p_error, sl3_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(15), sl3_2p_order, sl3_2p_wave, sl3_2p_flux, sl3_2p_error, sl3_2p_bit, format = 'i,f,f,f,i', /silent
	
	readcol, nodfiles(2), ll1_1p_order, ll1_1p_wave, ll1_1p_flux, ll1_1p_error, ll1_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(3), ll1_2p_order, ll1_2p_wave, ll1_2p_flux, ll1_2p_error, ll1_2p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(4), ll2_1p_order, ll2_1p_wave, ll2_1p_flux, ll2_1p_error, ll2_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(5), ll2_2p_order, ll2_2p_wave, ll2_2p_flux, ll2_2p_error, ll2_2p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(6), ll3_1p_order, ll3_1p_wave, ll3_1p_flux, ll3_1p_error, ll3_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(7), ll3_2p_order, ll3_2p_wave, ll3_2p_flux, ll3_2p_error, ll3_2p_bit, format = 'i,f,f,f,i', /silent
	
	readcol, cofiles(1), ll1_order, ll1_wave, ll1_flux, ll1_error, ll1_bit, format = 'i,f,f,f,i', /silent
	readcol, cofiles(2), ll2_order, ll2_wave, ll2_flux, ll2_error, ll2_bit, format = 'i,f,f,f,i', /silent
	readcol, cofiles(3), ll3_order, ll3_wave, ll3_flux, ll3_error, ll3_bit, format = 'i,f,f,f,i', /silent
	
	readcol, cofiles(5), sl1_order, sl1_wave, sl1_flux, sl1_error, sl1_bit, format = 'i,f,f,f,i', /silent
	readcol, cofiles(6), sl2_order, sl2_wave, sl2_flux, sl2_error, sl2_bit, format = 'i,f,f,f,i', /silent
	readcol, cofiles(7), sl3_order, sl3_wave, sl3_flux, sl3_error, sl3_bit, format = 'i,f,f,f,i', /silent

endelse

; Collapse orders into one large spectrum sorted by wavelength

allflux_1p   = [sl2_1p_flux, sl3_1p_flux, sl1_1p_flux, ll2_1p_flux, ll3_1p_flux, ll1_1p_flux]
allwave_1p   = [sl2_1p_wave, sl3_1p_wave, sl1_1p_wave, ll2_1p_wave, ll3_1p_wave, ll1_1p_wave]
allorder_1p  = [sl2_1p_order, sl3_1p_order, sl1_1p_order, ll2_1p_order, ll3_1p_order, ll1_1p_order]
allerror_1p  = [sl2_1p_error, sl3_1p_error, sl1_1p_error, ll2_1p_error, ll3_1p_error, ll1_1p_error]
allbit_1p    = [sl2_1p_bit, sl3_1p_bit, sl1_1p_bit, ll2_1p_bit, ll3_1p_bit, ll1_1p_bit]

allflux_2p   = [sl2_2p_flux, sl3_2p_flux, sl1_2p_flux, ll2_2p_flux, ll3_2p_flux, ll1_2p_flux]
allwave_2p   = [sl2_2p_wave, sl3_2p_wave, sl1_2p_wave, ll2_2p_wave, ll3_2p_wave, ll1_2p_wave]
allorder_2p  = [sl2_2p_order, sl3_2p_order, sl1_2p_order, ll2_2p_order, ll3_2p_order, ll1_2p_order]
allerror_2p  = [sl2_2p_error, sl3_2p_error, sl1_2p_error, ll2_2p_error, ll3_2p_error, ll1_2p_error]
allbit_2p    = [sl2_2p_bit, sl3_2p_bit, sl1_2p_bit, ll2_2p_bit, ll3_2p_bit, ll1_2p_bit]

allwave_1p  = allwave_1p(sort(allwave_1p))
allflux_1p  = allflux_1p(sort(allwave_1p))
allorder_1p = allorder_1p(sort(allwave_1p))
allerror_1p = allerror_1p(sort(allwave_1p))
allbit_1p   = allbit_1p(sort(allwave_1p))

allwave_2p  = allwave_2p(sort(allwave_2p))
allflux_2p  = allflux_2p(sort(allwave_2p))
allorder_2p = allorder_2p(sort(allwave_2p))
allerror_2p = allerror_2p(sort(allwave_2p))
allbit_2p   = allbit_2p(sort(allwave_2p))

; Coadd the nod positions using a weighted mean 

sl1_coadd = sl1_flux
sl2_coadd = sl2_flux
sl3_coadd = sl3_flux
ll1_coadd = ll1_flux
ll2_coadd = ll2_flux
ll3_coadd = ll3_flux

sl1_coadd_err = sl1_error
sl2_coadd_err = sl2_error
sl3_coadd_err = sl3_error
ll1_coadd_err = ll1_error
ll2_coadd_err = ll2_error
ll3_coadd_err = ll3_error

allflux_coadd = [sl2_coadd, sl3_coadd, sl1_coadd, ll2_coadd, ll3_coadd, ll1_coadd]
allerror_coadd = [sl2_coadd_err, sl3_coadd_err, sl1_coadd_err, ll2_coadd_err, ll3_coadd_err, ll1_coadd_err]

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

linelist = [5.511, 6.2, 6.909, 6.98, 7.7, 8.025, 8.6, 8.99138, 9.665, $
	 10.511, 11.3, 12.279, 12.6, 12.814, $
	  14.2,  14.322, 15.555, 16.4, 17.035, 17.4, $
	18.713,  24.318, 25.890, 28.218]

line_id = ['H!I2!N S(7)', 'PAH', 'H!I2!N S(5)', 'ArII', 'PAH', 'H!I2!N S(4)', 'PAH', 'ArIII', 'H!I2!N S(3)', $
	'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
	  'PAH',  'NeV', 'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', $
	'SIII', 'NeV', 'OIV', 'H!I2!N S(0)']

; Offset factor to plot 2nd nod position

offset = 4.

; Plot data

; Hard copy option

if keyword_set(ps) then begin
	set_plot,'ps'
	device, /landscape, filename = plotdir+strjoin(strsplit(obj,' ',/extract))+'_spect_lr.ps', /color
	defcolor = fsc_color('Black')
	lthick = 2
	cs = 1
endif else begin
	defcolor = fsc_color('White')
	lthick = 1
	cs = 2
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
green = fsc_color("Green")
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")

if not keyword_set(xr) then xr = [4,42]
if not keyword_set(yr) then yr = [1d-4,max(allflux_2p*offset)]

plot, sl1_1p_wave, sl1_1p_flux, $
	/xlog, $
	/ylog, $
;	xticks = 13, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
	xrange = xr, /xstyle, $
	yrange = yr, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Counts [Jy]', $
	title = obj, $
	charsize = cs, $
	color = defcolor, $
	thick = lthick, $
	charthick = lthick, $
	/nodata

if not keyword_set(bw) then begin
	if not keyword_set(nods) then begin
		oplot, sl1_1p_wave / (redshift + 1.),sl1_coadd, psym = 10, color = red, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.),sl2_coadd, psym = 10, color = blue, thick = lthick
		oplot, sl3_1p_wave / (redshift + 1.),sl3_coadd, psym = 10, color = orange, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.),ll1_coadd, psym = 10, color = yellow, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.),ll2_coadd, psym = 10, color = green, thick = lthick
		oplot, ll3_1p_wave / (redshift + 1.),ll3_coadd, psym = 10, color = orange, thick = lthick
;		xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
		legend, ['SL1', 'SL2', 'LL1', 'LL2', 'Bonus orders'], linestyle = [0,0,0,0,0], $
			color = [red, blue, yellow, green, orange], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
	endif else begin
		legend, ['SL1, pos 1', 'SL2, pos 1', 'LL1, pos 1', 'LL2, pos 1', 'Bonus orders, pos 1', $
			'SL1, pos 2', 'SL2, pos 2', 'LL1, pos 2', 'LL2, pos 2', 'Bonus orders, pos 2'], $
			linestyle = [0,0,0,0,0,1,1,1,1,1], $
			color = [red, blue, yellow, green, orange, red, blue, yellow, green, orange], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
		oplot, sl1_1p_wave / (redshift + 1.), sl1_1p_flux, psym = 10, color = red, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.), sl2_1p_flux, psym = 10, color = blue, thick = lthick
		oplot, sl3_1p_wave / (redshift + 1.), sl3_1p_flux, psym = 10, color = orange, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.), ll1_1p_flux, psym = 10, color = yellow, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_1p_flux, psym = 10, color = green, thick = lthick
		oplot, ll3_1p_wave / (redshift + 1.), ll3_1p_flux, psym = 10, color = orange, thick = lthick
		oplot, sl1_2p_wave / (redshift + 1.), sl1_2p_flux*offset, psym = 10, color = red, linestyle = 1, thick = lthick
		oplot, sl2_2p_wave / (redshift + 1.), sl2_2p_flux*offset, psym = 10, color = blue, linestyle = 1, thick = lthick
		oplot, sl3_2p_wave / (redshift + 1.), sl3_2p_flux*offset, psym = 10, color = orange, linestyle = 1, thick = lthick
		oplot, ll1_2p_wave / (redshift + 1.), ll1_2p_flux*offset, psym = 10, color = yellow, linestyle = 1, thick = lthick
		oplot, ll2_2p_wave / (redshift + 1.), ll2_2p_flux*offset, psym = 10, color = green, linestyle = 1, thick = lthick
		oplot, ll3_2p_wave / (redshift + 1.), ll3_2p_flux*offset, psym = 10, color = orange, linestyle = 1, thick = lthick
	endelse
endif else begin
	if not keyword_set(nods) then begin
		oplot, sl1_1p_wave / (redshift + 1.), sl1_coadd, psym = 10, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.), sl2_coadd, psym = 10, thick = lthick
		oplot, sl3_1p_wave / (redshift + 1.), sl3_coadd, psym = 10, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.), ll1_coadd, psym = 10, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_coadd, psym = 10, thick = lthick
		oplot, ll3_1p_wave / (redshift + 1.), ll3_coadd, psym = 10, thick = lthick
;		xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
	endif else begin
		oplot, sl1_1p_wave / (redshift + 1.), sl1_1p_flux, psym = 10, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.), sl2_1p_flux, psym = 10, thick = lthick
		oplot, sl3_1p_wave / (redshift + 1.), sl3_1p_flux, psym = 10, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.), ll1_1p_flux, psym = 10, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_1p_flux, psym = 10, thick = lthick
		oplot, ll3_1p_wave / (redshift + 1.), ll3_1p_flux, psym = 10, thick = lthick
		oplot, sl1_2p_wave / (redshift + 1.), sl1_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		oplot, sl2_2p_wave / (redshift + 1.), sl2_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		oplot, sl3_2p_wave / (redshift + 1.), sl3_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		oplot, ll1_2p_wave / (redshift + 1.), ll1_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		oplot, ll2_2p_wave / (redshift + 1.), ll2_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		oplot, ll3_2p_wave / (redshift + 1.), ll3_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		legend, ['Nod position 1', 'Nod position 2'], linestyle = [0,1], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
	endelse
endelse

if not keyword_set(nolines) then begin
	for i = 0, n_elements(linelist) - 1 do begin
		off = i mod 2
		plots, [linelist(i),linelist(i)],$
			[allflux_coadd(closeto(allwave_1p/(redshift + 1d),linelist(i))),$
			allflux_coadd(closeto(allwave_1p/(redshift + 1d),linelist(i)))*2],$
			color=defcolor, linestyle = 1
		lineend = (allflux_coadd(closeto(allwave_1p/(redshift + 1d),linelist(i)))*2)
		xyouts, linelist(i), lineend + 5*off*lineend, line_id(i), $
			orientation = 90, charsize = cs, /data

		; 6 um H2O absorption

		plots,[5.5,6.5],[allflux_coadd(closetomed(allwave_1p/(1 + redshift),allflux_coadd,6.0))*0.5,$
			allflux_coadd(closetomed(allwave_1p/(1 + redshift),allflux_coadd,6.0))*0.5]
		xyouts, 6.0, allflux_coadd(closetomed(allwave_1p/(1 + redshift),allflux_coadd,6.0))*0.3,'H!I2!NO', $
			charsize = cs

		; 9.7 um silicate absorption

		plots,[8,12],[allflux_coadd(closetomed(allwave_1p/(1 + redshift),allflux_coadd,9.7))*0.5,$
			allflux_coadd(closetomed(allwave_1p/(1 + redshift),allflux_coadd,9.7))*0.5]
		xyouts, 9.7, allflux_coadd(closetomed(allwave_1p/(1 + redshift),allflux_coadd,9.7))*0.3,'Silicate', $
			charsize = cs
		;xyouts, linelist(i), 0.35*yr(1) + 1.05*off, line_id(i), orientation = 90, charsize = cs, /data
	endfor
endif

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

stop
end
