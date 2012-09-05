pro hrstitch, fname, xr = xr, yr = yr, nolines = nolines, ps = ps, stop = stop, hrsky = hrsky
;+
; NAME:
;       HRSTITCH
;
; PURPOSE:
; 	Display trimmed hi-res Spitzer IRS spectra
;
; OUTPUTS:
;	- Plots individual spectra, overplotting both optimal and regular extractions
;
; KEYWORDS:
;
;	NOLINES - removes indicators of common emission and absorption features in ULIRGs. Default
;			is to display them.
;
;	PS - create a hard copy of spectra plot
;
;	HRSKY - plots LH modules with sky background removed, if available
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;
; EXAMPLE:
;	IDL> hrstitch, 'mega023'
;
; REVISION HISTORY
;	Adapted from FULLQL_HR			KW, Feb 08
;	Added HRSKY keyword - KW, May 08
;-

;device, decomposed = 1

; Load the stitched, calibrated data from the IDL sav files

tag, fname, dirtag

if keyword_set(hrsky) then begin
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/hrsky/'+fname+'.sav'
endif else begin
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
endelse

plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibrated/'

; Identify IR lines for overplotting

templines = ir_lines(/hr)
lines = templines(*,0) & line_id = templines(*,1)

; Create indices for spectra arranged by increasing wavelength (removing any jumps due to overlapping spectra in different orders)

; Set default plot ranges (SH)

if not keyword_set(xr) then xr = [8,18]
if not keyword_set(yr) then yr = [-0.1,1.]

; Hard copy option

if keyword_set(ps) then begin
	set_plot,'ps'
	device, /landscape, /color, $
		filename = plotdir+strjoin(strsplit(sed.obj,' ',/extract))+'_spect_hr.ps'
	defcolor = fsc_color('Black')
	lthick = 5
	cthick = 4
	cs = 1.5
endif else begin
	defcolor = fsc_color('White')
	lthick = 1
	cthick = 1
	cs = 1.5
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

; Plot the data in rest frame wavelengths 

;!p.multi = [0,1,1]

plot, sed.wave_sh, sed.wave_lh, $
	xrange = xr, $
	/xstyle, $
	yrange = yr, $
	/ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = cs, $
	color = defcolor, $
	thick = lthick, charthick = cthick, $
	xthick = lthick, $
	ythick = lthick, $
	/nodata

; Plot coadded data

	for i = 0, 9 do begin

		; SH

		ajunk0 = where(sed.order_sh eq (11 + i))
		cjunk0 = sort(sed.wave_sh(ajunk0))
		
		; LH
		
		ajunk2 = where(sed.order_lh eq (11 + i))
		cjunk2 = sort(sed.wave_lh(ajunk2))

		oplot, sed.wave_sh(ajunk0(cjunk0)), sed.flux_sh(ajunk0(cjunk0)), $
			color = colors(i), psym = 10, thick = lthick

		oplot, sed.wave_lh(ajunk2(cjunk2)), sed.flux_lh(ajunk2(cjunk2)), $
			color = colors(i), psym = 10, thick = lthick
	endfor
	
; Overplot vertical lines with labels indicating locations of likely absorption and emission features

if not keyword_set(nolines) then begin
	for i = 0, n_elements(lines) - 1 do begin
		off = i mod 2
		ver, lines(i), linestyle = 1, color = defcolor
		xyouts, lines(i), 0.85*yr(1) + 0.05*off*yr(1), line_id(i), orientation = 90, charsize = 1.5, /data, $
			color = defcolor, charthick = cthick
	endfor
endif

; Hard copy

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif


if keyword_set(stop) then stop
end
