pro shcal, fname, xr = xr, yr = yr, nods = nods, nolines = nolines, ps = ps, notrim = notrim, no_offset = no_offset

;+
; NAME:
;       SHCAL
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
;	PS - create a hard copy of spectra plot
;
;	NO_OFFSET - overplot nods in true units rather than offsetting them
;
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;
; EXAMPLE:
;	IDL> shcal, 'mega023'
;
; REVISION HISTORY
;       Written by K. Willett                Sep 2007
;-

; Set device so that colors plotted for different orders appear on the screen

device, decomposed = 1, window_state = state

; Locate directory from which to read (either megamasers or CSOs)

targets, fname, redshift, obj
tag, fname, dirtag


writedir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/'
plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/'
opath = '~/spice/output/hires/optimal/all/'

plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibrated/'
copath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
nodpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'

co_flist = copath+'*'+fname+'*.tbl'
nod_flist = nodpath+'*'+fname+'*.tbl'

cofiles = findfile(co_flist)
conum = n_elements(cofiles)
nodfiles = findfile(nod_flist)
nodnum = n_elements(nodfiles)

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

; Use the spectroscopic redshifts of the objects (from NED) to plot images in the target rest frame

	readcol, cofiles(4), sh_order, sh_wave, sh_flux, sh_err, sh_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(8), sh_1p_order, sh_1p_wave, sh_1p_flux, sh_1p_err, sh_1p_bit, format = 'i,f,f,f,i', /silent
	readcol, nodfiles(9), sh_2p_order, sh_2p_wave, sh_2p_flux, sh_2p_err, sh_2p_bit, format = 'i,f,f,f,i', /silent

; Coadd the nod positions using a weighted mean based on the spectral uncertainty

sh_coadd = sh_flux

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

linelist = [8.6, 8.99138, 9.665, 10.511, 11.3, 12.279, 12.6, 12.814, 13.7, $
	14.0, 14.2, 14.322, 14.368, 15.0, 15.555, 16.4, 17.035]

line_id = ['PAH', 'ArIII', 'H!I2!N S(3)', 'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
	'C!I2!NH!I2!N', 'HCN', 'PAH', 'NeV', 'ClII', 'CO!I2!N', 'NeIII', 'PAH', 'H!I2!N S(1)']

linestart = where(min(sh_wave/(redshift + 1)) gt linelist)
lineterm = where(max(sh_wave/(redshift + 1)) lt linelist)

if linestart(0) ne -1 then begin
       lsind = [linestart(0)+1,n_elements(linelist)-1]
       linelist = linelist(lsind(0):lsind(1)) 
       line_id  = line_id (lsind(0):lsind(1))
       print,'ls'
endif
if lineterm(0)  ne -1 then begin
	ltind = [0,lineterm(0)-1]
	linelist = linelist(ltind(0):ltind(1)) 
	line_id = line_id(ltind(0):ltind(1))
	print,'lt'
endif

; Create indices for spectra arranged by increasing wavelength (removing any jumps due to overlapping spectra in different orders)

b0 = sort(sh_1p_wave)
b1 = sort(sh_2p_wave)

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

; Plot the data in rest frame wavelengths 

if not keyword_set(xr) then xr = [min(sh_wave/(redshift+1))-0.5,max(sh_wave/(redshift+1))+0.5]
if not keyword_set(yr) then yr = [-0.1,max(sh_flux)+0.15]

plot, sh_1p_wave(b0)/ (redshift + 1), sh_1p_flux(b0), $
	xrange = xr, $
	/xstyle, $
	yrange = yr, $
	/ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Counts [Jy]', $
	title = obj, $
	charsize = 2, $
	color = defcolor, $
	thick = lthick, charthick = cthick, $
	/nodata

; Plot coadded data

if not keyword_set(nods) then begin

	wavesort = sort(sh_wave)
	oplot, sh_wave(wavesort)/(redshift + 1.), sh_flux(wavesort), $
		psym = 10, thick = lthick

endif else begin
	
; Plot individual nod positions

		; Plot data for the 1st nod position

		ws1 = sort(sh_1p_wave)
		oplot, sh_1p_wave(ws1)/(redshift + 1.), sh_1p_flux(ws1), $
			linestyle = 0, psym = 10, thick = lthick

		; Plot data for the 2nd nod position (may be offset from the first)
	
		if keyword_set(no_offset) then offset = 0 else offset = 0.2

		ws2 = sort(sh_1p_wave)
		oplot, sh_2p_wave(ws2)/(redshift + 1.), sh_2p_flux(ws2) + offset, $
			linestyle = 0, psym = 10, thick = lthick
endelse

; Overplot vertical lines with labels indicating locations of likely absorption and emission features

if not keyword_set(nolines) then begin
	for i = 0, n_elements(linelist) - 1 do begin
		off = i mod 2
		plots, [linelist(i),linelist(i)],$
			[sh_coadd(closeto(sh_wave/(redshift + 1d),linelist(i))),$
			sh_coadd(closeto(sh_wave/(redshift + 1d),linelist(i)))+0.1],$
			color=defcolor, linestyle = 1
		xyouts, linelist(i), (sh_coadd(closeto(sh_wave/(redshift + 1d),linelist(i)))+0.1) + 0.05*off, line_id(i), $
			orientation = 90, charsize = cs, /data
	endfor
endif

; Hard copy

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif


stop
end

