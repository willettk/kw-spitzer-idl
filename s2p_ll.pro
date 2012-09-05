pro s2p_ll, fname, nostop = nostop, color = color, coadd = coadd, ps = ps, lines = lines, $
	xr = xr, yr = yr, pahfit = pahfit, write = write, notrim = notrim, ss = ss, noclean = noclean, $
	regular = regular
;+
; NAME:
;       S2P_LL
;
; PURPOSE:
; 	Display trimmed LL Spitzer IRS spectra with option of fitting ULIRG template via PAHFIT
;
; INPUTS:
;	XR - plot range in x-direction (log scale)
;
;	YR - plot range in y-direction (log scale)
;
; OUTPUTS:
; 	- Plots low-res spectra sorted by order in rest frame, with optional co-adding and line ID
;
; KEYWORDS:
;
;	NOSTOP - removes STOP statement at end of program
;
;	COLOR - plots individual orders in different colors. Default is to plot all orders stitched together. 
;
;	COADD - plots the coadded spectrum weighted by the uncertainty in each pixel. Default is
;		to plot both nod positions with one offset by an arbitrary multiplicative factor. 
;
;	PS - plots the spectrum to a hardcopy postscript file
;
;	LINES - overplots common IR emission lines with ID labels
;
;	NO_CLEAN - view data (if available) that has not gone through IRSCLEAN_MASK. Default
;			is to view cleaned data. 
;
;	REGULAR - view and write data extracted using nominal settings in SPICE (if present). Default
;			is to use the optimally extracted data. 
;
; EXAMPLE:
;	IDL> s2p_ll, 'mega023'
;
; REVISION HISTORY
;       Adapted from SPEC2PAHFIT.pro                KW, Aug 2007
;-

; Set device to read in colors

if not keyword_set(ps) then device, decomposed = 1, window_state = state 
device, window_state = state

; Locate directory from which to read (either megamasers or CSOs)

tag, fname, dirtag 

; Read in the optimally extracted data from SPICE directory

if keyword_set(regular) then begin
	writedir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/regular/'
	plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/regular/'
	opath = '~/spice/output/lores/regular/all/' 
	pahpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/regular/'
endif else begin
	writedir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/'
	plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/'
	opath = '~/spice/output/lores/optimal/all/'
	pahpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/'
endelse

if keyword_set(noclean) then files_opt = opath+'*'+fname+'*ch2*p_spect.tbl'$
	else files_opt = opath+'*'+fname+'*ch2*clean*spect.tbl'

;if keyword_set(clean) then files_opt = opath+'*'+fname+'*clean*spect.tbl' $
;	else files_opt = opath+'*'+fname+'*p_spect.tbl'

ofiles = findfile(files_opt)
onum = n_elements(ofiles)

; Stored object designations, redshifts for the Spitzer sample

targets, fname, redshift, obj

; Read information in from tables for both nod positions

readcol, ofiles(0), ll1_1p_order, ll1_1p_wave, ll1_1p_flux, ll1_1p_error, ll1_1p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(1), ll1_2p_order, ll1_2p_wave, ll1_2p_flux, ll1_2p_error, ll1_2p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(2), ll2_1p_order, ll2_1p_wave, ll2_1p_flux, ll2_1p_error, ll2_1p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(3), ll2_2p_order, ll2_2p_wave, ll2_2p_flux, ll2_2p_error, ll2_2p_bit, format = 'i,f,f,f,i', /silent

; Sort the pixels so that they are in order by wavelength

ll1_1p_index = sort(ll1_1p_wave)
ll2_1p_index = sort(ll2_1p_wave)

ll1_2p_index = sort(ll1_2p_wave)
ll2_2p_index = sort(ll2_2p_wave)

; Trim the ends of the modules - template based on SL spectrum for mega023 (IRAS 16255+2801)

trim_file = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/trim_templates/'+fname+'_lores.trim')
	if trim_file(0) eq '' then begin

		ll2_edges = [3,15]
		ll1_edges = [3,25]

		print,'Using standard template for lo-res modules on '+fname

	endif else begin

		print,'Using saved template for lo-res modules on '+fname
		readcol, trim_file(0), trimstart, trimend, format = 'x,i,i', skipline = 2, /silent

		ll2_edges = [trimstart(2),trimend(2)]	
		ll1_edges = [trimstart(3),trimend(3)]	

	endelse



if keyword_set(notrim) then begin
	ll1_edges = [0,0]
	ll2_edges = [0,0]
endif

ll1_1p_index = ll1_1p_index(ll1_edges(0):n_elements(ll1_1p_index)-ll1_edges(1)-1)
ll2_1p_index = ll2_1p_index(ll2_edges(0):n_elements(ll2_1p_index)-ll2_edges(1)-1)

ll1_2p_index = ll1_2p_index(ll1_edges(0):n_elements(ll1_2p_index)-ll1_edges(1)-1)
ll2_2p_index = ll2_2p_index(ll2_edges(0):n_elements(ll2_2p_index)-ll2_edges(1)-1)

; New orders with trimmed edges

ll1_1p_wave = ll1_1p_wave(ll1_1p_index) & ll1_1p_flux = ll1_1p_flux(ll1_1p_index) & ll1_1p_error = ll1_1p_error(ll1_1p_index) & $
	ll1_1p_order = ll1_1p_order(ll1_1p_index) & ll1_1p_bit = ll1_1p_bit(ll1_1p_index)
ll2_1p_wave = ll2_1p_wave(ll2_1p_index) & ll2_1p_flux = ll2_1p_flux(ll2_1p_index) & ll2_1p_error = ll2_1p_error(ll2_1p_index) & $
	ll2_1p_order = ll2_1p_order(ll2_1p_index) & ll2_1p_bit = ll2_1p_bit(ll2_1p_index)

ll1_2p_wave = ll1_2p_wave(ll1_2p_index) & ll1_2p_flux = ll1_2p_flux(ll1_2p_index) & ll1_2p_error = ll1_2p_error(ll1_2p_index) & $
	ll1_2p_order = ll1_2p_order(ll1_2p_index) & ll1_2p_bit = ll1_2p_bit(ll1_2p_index)
ll2_2p_wave = ll2_2p_wave(ll2_2p_index) & ll2_2p_flux = ll2_2p_flux(ll2_2p_index) & ll2_2p_error = ll2_2p_error(ll2_2p_index) & $
	ll2_2p_order = ll2_2p_order(ll2_2p_index) & ll2_2p_bit = ll2_2p_bit(ll2_2p_index)

; Collapse orders into one large spectrum sorted by wavelength

allflux_1p   = [ll2_1p_flux, ll1_1p_flux]
allwave_1p   = [ll2_1p_wave, ll1_1p_wave]
allorder_1p  = [ll2_1p_order, ll1_1p_order]
allerror_1p  = [ll2_1p_error, ll1_1p_error]
allbit_1p    = [ll2_1p_bit, ll1_1p_bit]

allflux_2p   = [ll2_2p_flux, ll1_2p_flux]
allwave_2p   = [ll2_2p_wave, ll1_2p_wave]
allorder_2p  = [ll2_2p_order, ll1_2p_order]
allerror_2p  = [ll2_2p_error, ll1_2p_error]
allbit_2p    = [ll2_2p_bit, ll1_2p_bit]

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

weights_ll1_1p = 1d-6/ll1_1p_error^2
weights_ll2_1p = 1d-6/ll2_1p_error^2

weights_ll1_2p = 1d-6/ll1_2p_error^2
weights_ll2_2p = 1d-6/ll2_2p_error^2

ll1_coadd = (weights_ll1_1p * ll1_1p_flux + weights_ll1_2p * ll1_2p_flux) / (weights_ll1_1p + weights_ll1_2p)
ll2_coadd = (weights_ll2_1p * ll2_1p_flux + weights_ll2_2p * ll2_2p_flux) / (weights_ll2_1p + weights_ll2_2p)

ll1_coadd_err = sqrt(ll1_1p_error^2 + ll1_2p_error^2)
ll2_coadd_err = sqrt(ll2_1p_error^2 + ll2_2p_error^2)

allflux_coadd = [ll2_coadd, ll1_coadd]
allerror_coadd = sqrt(allerror_1p^2 + allerror_2p^2)

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

linelist = [5.511, 6.1088, 6.2, 6.909, 6.98, 7.7, 8.025, 8.6, 9.665, $
	 10.511, 11.3, 12.279, 12.6, 12.814, $
	  14.2,  15.555, 16.4, 17.035, 17.4, $
	18.713,  28.218]

line_id = ['H!I2!N S(7)', 'H!I2!N S(6)', 'PAH', 'H!I2!N S(5)', 'ArII', 'PAH', 'H!I2!N S(4)', 'PAH', 'H!I2!N S(3)', $
	'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
	  'PAH',  'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', $
	'SIII', 'H!I2!N S(0)']

; Old lineset - includes many weaker medium-ionization lines usually not resolvable in lo-res spectra

;linelist = [5.34, 5.511, 6.2, 6.909, 6.98, 7.65, 7.7, 8.025, 8.6, 9.665, 10.511, 11.3, 12.279, 12.6, 12.814, 13.7, $
;	14.0, 14.2, 14.322, 14.368, 15.0, 15.555, 16.4, 17.035, 17.4, 17.934, 18.713, 24.318, 25.890, $
;	25.988, 28.218, 33.481, 34.815]

;line_id = ['FeII', 'H!I2!N S(7)', 'PAH', 'H!I2!N S(5)', 'ArII', 'NeVI', 'PAH', 'H!I2!N S(4)', 'PAH', 'H!I2!N S(3)', $
;	'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
;	'C!I2!NH!I2!N', 'HCN', 'PAH', 'NeV', 'ClII', 'CO!I2!N', 'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', 'FeII', 'SIII', $
;	'NeV', 'OIV', 'FeII', 'H!I2!N S(0)', 'SIII', 'SiII']

; Offset factor to plot 2nd nod position

offset = 4.

; Plot data

if state(0) eq 0 then window,0 else wset,0
!p.multi = [0,1,1]

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

if not keyword_set(xr) then xr = [4,42]
if not keyword_set(yr) then yr = [1d-4,max(allflux_2p*offset)]

plot, ll1_1p_wave, ll1_1p_flux, $
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

if keyword_set(color) then begin
	if keyword_set(coadd) then begin
		oplot, ll1_1p_wave / (redshift + 1.),ll1_coadd, psym = 10, color = fsc_color('Yellow'), thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.),ll2_coadd, psym = 10, color = fsc_color('Green'), thick = lthick
		xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
		legend, ['LL1', 'LL2'], linestyle = [0,0], $
			color = [fsc_color('Yellow'), fsc_color('Green')], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
	endif else begin
		legend, ['LL1, pos 1', 'LL2, pos 1', 'LL1, pos 2', 'LL2, pos 2'], $
			linestyle = [0,0,1,1], $
			color = [fsc_color('Yellow'), fsc_color('Green'),$
			fsc_color('Yellow'), fsc_color('Green')], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
		oplot, ll1_1p_wave / (redshift + 1.), ll1_1p_flux, psym = 10, color = fsc_color('Yellow'), thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_1p_flux, psym = 10, color = fsc_color('Green'), thick = lthick
		oplot, ll1_2p_wave / (redshift + 1.), ll1_2p_flux*offset, psym = 10, color = fsc_color('Yellow'), linestyle = 1, thick = lthick
		oplot, ll2_2p_wave / (redshift + 1.), ll2_2p_flux*offset, psym = 10, color = fsc_color('Green'), linestyle = 1, thick = lthick
	endelse
endif else begin
	if keyword_set(coadd) then begin
		oplot, ll1_1p_wave / (redshift + 1.), ll1_coadd, psym = 10, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_coadd, psym = 10, thick = lthick
		xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
	endif else begin
		oplot, ll1_1p_wave / (redshift + 1.), ll1_1p_flux, psym = 10, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_1p_flux, psym = 10, thick = lthick
		oplot, ll1_2p_wave / (redshift + 1.), ll1_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		oplot, ll2_2p_wave / (redshift + 1.), ll2_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
		legend, ['Nod position 1', 'Nod position 2'], linestyle = [0,1], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
	endelse
endelse

if keyword_set(lines) then begin
	for i = 0, n_elements(linelist) - 1 do begin
		off = i mod 2
		vline, linelist(i), linestyle = 1, color = defcolor, /data, /noerase, range = [yr(0),0.3*yr(1)]
		xyouts, linelist(i), 0.35*yr(1) + 1.05*off, line_id(i), orientation = 90, charsize = cs, /data
	endfor
endif

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif


; Run PAHFIT on coadded spectrum

if keyword_set(pahfit) then begin

	resolve_routine, 'pahfit'

	if state(2) eq 0 then window, 2, xsize = 800, ysize = 500 else wset,2
	
	; Add code to create report file if none exists (check)

	pahfile = pahpath+strjoin(strsplit(obj,' ',/extract))+'_pahfit.txt'
	findpahfile = findfile(pahfile)
	if findpahfile(0) eq '' then spawn,'touch '+pahfile

	fit = pahfit(allwave_1p, allflux_coadd, allerr_coadd, redshift = redshift, /plot_progress, $
		report = pahpath+strjoin(strsplit(obj,' ',/extract))+'_pahfit.txt', xsize = 800, ysize = 500)
	wset,0
endif

; Write the new, trimmed spectra to file

if keyword_set(write) then begin

	; All modules

	forprint, allorder_1p, allwave_1p, allflux_1p, allerror_1p, allbit_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]      Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_lr_1p.tbl', /silent

	forprint, allorder_2p, allwave_2p, allflux_2p, allerror_2p, allbit_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]      Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_lr_2p.tbl', /silent

	forprint, allorder_1p, allwave_1p, allflux_coadd, allerror_coadd, allbit_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]      Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_lr_coadd.tbl', /silent

	; Individual modules

	forprint, ll1_1p_order, ll1_1p_wave, ll1_1p_flux, ll1_1p_error, ll1_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll1_1p.tbl', /silent

	forprint, ll2_1p_order, ll2_1p_wave, ll2_1p_flux, ll2_1p_error, ll2_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll2_1p.tbl', /silent

	forprint, ll1_2p_order, ll1_2p_wave, ll1_2p_flux, ll1_2p_error, ll1_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll1_2p.tbl', /silent

	forprint, ll2_2p_order, ll2_2p_wave, ll2_2p_flux, ll2_2p_error, ll2_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll2_2p.tbl', /silent

	; Coadded

	forprint, ll1_1p_order, ll1_1p_wave, ll1_coadd, ll1_1p_error, ll1_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_ll1_coadd.tbl', /silent

	forprint, ll2_1p_order, ll2_1p_wave, ll2_coadd, ll2_1p_error, ll2_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Counts [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_ll2_coadd.tbl', /silent

endif

; Temporary section for assessing S/N in portions of coadded spectra

;redwave = ll1_1p_wave/(redshift+1.)
;sn_start = where(abs(redwave - 20.) eq min(abs(redwave - 20.)))
;sn_end   = where(abs(redwave - 30.) eq min(abs(redwave - 30.)))

;sn = alog(ll1_coadd(sn_start(0):sn_end(0)))
;snwave = redwave(sn_start(0):sn_end(0))
;snerr = ll1_coadd_err(sn_start(0):sn_end(0))
;expr = 'p[0] + p[1]*x'
;start = [0d,0.1d]
;result = mpfitexpr(expr, snwave, sn, snerr, start,/quiet)

;fitted_sn = sn - (result(0) + snwave*result(1))
;sn_dev = stddev(fitted_sn)

;print,''
;print,'S/N = ',sn_dev
;print,''

if not keyword_set(nostop) then stop
end

