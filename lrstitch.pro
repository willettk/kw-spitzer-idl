pro lrstitch, fname, bw = bw, ps = ps, nolines = nolines, $
	xr = xr, yr = yr, bonus = bonus, stop = stop, hydlines = hydlines, $
	lin = lin, allpah = allpah
;+
; NAME:
;       LRSTITCH
;
; PURPOSE:
; 	Display stitched and calibrated lo-res Spitzer IRS spectra
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
;	HYDLINES - show hydrocarbon absorption locations.
;
;	PS - plots the spectrum to a hardcopy postscript file
;
;	BONUS - plots the bonus orders SL3 and LL3
;
;	LIN - displays the y-axis as a linear scale (default is log)
;
;	ALLPAH - vertical lines at all known major PAH locations
;
; EXAMPLE:
;	IDL> lrstitch, 'mega023'
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;
; NOTES:
;
; REVISION HISTORY
;       Adapted from SPEC2PAHFIT		- KW, Feb 08
; 	Added ALLPAH - Mar 08
;	Changed ALLPAH label height - Aug 09
;-

; Load the stitched, calibrated data from the IDL sav files

tag, fname, dirtag
restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

flux = sed.flux_lr
wave = sed.wave_lr
order = sed.order_lr

noneg = where(sed.flux_lr gt 0)
flux = flux[noneg]
wave = wave[noneg]
order = order[noneg]

plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibrated/'

; Identified IR lines for overplotting

templines = ir_lines(/lr)
lines = templines(*,0) & line_id = templines(*,1)

; Plot data

; Hard copy option

if keyword_set(ps) then begin
	set_plot,'ps'
	device, /landscape, filename = plotdir+strjoin(strsplit(sed.obj,' ',/extract))+'_spect_lr.ps', /color
	defcolor = fsc_color('Black')
	lthick = 5
	cthick = 5
	cs = 1.5
endif else begin
	defcolor = fsc_color('White')
	lthick = 1
	cthick = 1
	cs = 1.5
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
green = fsc_color("Green")
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")

if not keyword_set(xr) then xr = [4,42]
if not keyword_set(yr) then yr = [1d-4,max(flux)]

sl_bonus = where(order eq 3 and wave lt 10.)
ll_bonus = where(order eq 3 and wave gt 10.)

sl2_index = where((order eq 2) and (wave lt 10.))		; An ugly solution. There should be a better way of distinguishing modules.
sl1_index = where((order eq 1) and (wave lt 15.))
ll2_index = where((order eq 2) and (wave gt 10.))
ll1_index = where((order eq 1) and (wave gt 15.))

!p.multi=[0,1,1]

if not keyword_set(lin) then begin
	xlog = 1
	ylog = 1
endif else begin
	xlog = 0
	ylog = 0
endelse

plot, wave, flux, $
	xlog = xlog, $
	ylog = ylog, $
;	xticks = 13, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
	xrange = xr, /xstyle, $
	yrange = yr, /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = cs, $
	color = defcolor, $
	thick = lthick, $
	xthick = lthick, $
	ythick = lthick, $
	charthick = cthick, $
	/nodata

if not keyword_set(bw) then begin
		oplot, wave(sl1_index),flux(sl1_index), psym = 10, color = red, thick = lthick
		oplot, wave(sl2_index),flux(sl2_index), psym = 10, color = blue, thick = lthick
		oplot, wave(ll1_index),flux(ll1_index), psym = 10, color = yellow, thick = lthick
		oplot, wave(ll2_index),flux(ll2_index), psym = 10, color = green, thick = lthick
		if keyword_set(bonus) then begin
			oplot, wave(ll_bonus),flux(ll_bonus), psym = 10, color = orange, thick = lthick
			oplot, wave(sl_bonus),flux(sl_bonus), psym = 10, color = orange, thick = lthick
		endif
;		xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
;		legend, ['SL1', 'SL2', 'LL1', 'LL2', 'Bonus orders'], linestyle = [0,0,0,0,0], $
;			color = [red, blue, yellow, green, orange], $
;			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
endif else begin
		oplot, wave,flux, psym = 10, thick = lthick
;		xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
endelse

if not keyword_set(nolines) then begin
	for i = 0, n_elements(lines) - 1 do begin
		off = i mod 2
		ver, lines(i), linestyle = 1, color = defcolor
		xyouts, lines(i), 10d^((0.90 + 0.10*off)*alog10(yr[1]/yr[0]) + alog10(yr[0])), $
			line_id(i), orientation = 90, charsize = cs, /data
	endfor
	if keyword_set(hydlines) then begin
		ver,6.85,linestyle=2
		ver,7.25,linestyle=2
	endif
endif

if keyword_set(allpah) then begin
	pahlist = [5.3,5.7,6.2,7.7,8.6,11.3,12.7,14.2,16.4,17.1,17.4]
	for i=0, n_elements(pahlist) - 1 do begin
		off = i mod 2
		ver, pahlist[i], color=fsc_color("Red")
		xyouts, pahlist[i]*0.98, 10d^((0.80+ 0.10*off)*alog10(yr[1]/yr[0]) + alog10(yr[0])), $
			string(pahlist[i],format='(f4.1)'), orientation = 90, charsize = cs, /data, charthick = cthick
	endfor
endif

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif



if keyword_set(stop) then stop
end
