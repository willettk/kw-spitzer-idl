;+
; NAME:
;       CSO_PROP
;
; PURPOSE:
;
;	Plot overlaid lo-res spectra of CSOs
;
; INPUTS:
;
; OUTPUTS:
;
; KEYWORDS:
;
;	PS - hard copy of overlaid spectra
;
; EXAMPLE:
;
;	IDL> cso_prop, /ps
;
; REQUIRES:
;
; NOTES:
;
; REVISION HISTORY
;       Written by K. Willett                Oct 2007
;       Modified for CXC/HST proposal, Cycle 12		Mar 2010
;-

;pro cso_prop, ps = ps
ps = 1
; Set paths and restore the CSO data made with CSOMAKE

csodir = '~/Astronomy/Research/Spitzer/CSO/data/structures/'

restore, csodir+'cso001.sav'

; Hard copy option

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Chandra/cxc_spitzer_spectra.ps',/landscape,/color
	cs = 2
	ct = 5
	lightcolor=fsc_color("Magenta")
endif else begin
	lightcolor = fsc_color("Yellow")
endelse

; Plot data - each has a multiplicative scale factor (additive in log space) to allow the spectra to
; be easily distinguished

!p.multi=[0,1,1]

plot, sed.wave_lr, sed.flux_lr, $
	/nodata, $
	/xlog, $
	/ylog, $
	xrange = [4,51], /xstyle, $
	yrange = [1d-3, 1d3], /ystyle, $
	xtitle = '!7k!3!Irest!N [!7l!3m]', $
	ytitle = 'log f!I!7k!3!N (scaled)', $
	charsize = cs, $
	thick = ct, $
	xthick = ct, $
	ythick = ct, $
	charthick = ct, $
	xticks = 6, $
	xtickv = [5, 7, 10, 15, 20, 30, 40], $
	xtickname = strtrim([5, 7, 10, 15, 20, 30, 40],2), $
	yticks = 6, $
	ytickv = 10d^[-3,-2,-1,0,1,2,3], $
	ytickname = ['-3','-2','-1','0','1','2','3']


scale = 7.0
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Blue"), thick=ct
xyouts, 38, 3d-1,'NGC 3894', charthick=ct

restore, csodir+'cso002.sav'
scale = 6
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Magenta"), thick=ct
xyouts, 35, 4d-2,'4C +31.04', charthick=ct

restore, csodir+'cso004.sav'
scale = 1.5d2
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Orange Red"), thick=ct
xyouts, 37, 1d2,'OQ 208', charthick=ct

restore, csodir+'cso005.sav'
scale = 0.3
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Black"), thick=ct
xyouts, 31, 1d-2,'PKS 1413+135', charthick=ct

restore, csodir+'cso006.sav'
scale = 5d2
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Red"), thick=ct
xyouts, 35, 6d2,'4C +12.50', charthick=ct

restore, csodir+'cso007.sav'
scale = 1d3
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Green"), thick=ct
xyouts, 36, 9d0,'1946+708', charthick=ct

restore, csodir+'cso008.sav'
scale = 2d2
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Dark Green"), thick=ct
xyouts, 37, 4d0,'4C +37.11', charthick=ct

restore, csodir+'cso009.sav'
scale = 3d2
oplot, sed.wave_lr, sed.flux_lr*scale, color=fsc_color("Goldenrod"), thick=ct
xyouts, 36, 3d1,'NGC 6328', charthick=ct

; Show locations of the most common MIR PAH features

plots, [6.2,6.2],  [6d2,4d2],linestyle=0, thick=3
plots, [7.7,7.7],  [6d2,4d2],linestyle=0, thick=3
plots, [8.6,8.6],  [6d2,4d2],linestyle=0, thick=3
plots, [11.3,11.3],[6d2,4d2],linestyle=0, thick=3
plots, [12.7,12.7],[6d2,4d2],linestyle=0, thick=3

xyouts, 5.0, 4.5d2, 'PAHs', charthick=ct

;ver, 14.2, linestyle=2
;ver, 16.4, linestyle=2
;ver, 17.1, linestyle=2
;ver, 17.4, linestyle=2

; Mark common emission lines

plots, [15.555,15.555],[8d-3,3d-3]
xyouts, 15.2, 2d-3, 'NeIII', charthick=ct

plots, [18.713,18.713],[8d-3,3d-3]
xyouts, 18.2, 2d-3, 'SIII', charthick=ct

plots, [12.814,12.814],[8d-3,4d-3]
xyouts, 12.0, 2d-3, 'NeII', charthick=ct

plots, [17.04,17.04],[5d-3,7d-3]
xyouts, 16.85, 3d-3, 'H!I2!N', charthick=ct

plots, [6.91,6.91],[3.5d-3,5d-3]
xyouts, 6.9, 2d-3, 'H!I2!N', charthick=ct

; Mark locations of silicate absorption

xyouts, 9.0, 2d-3, 'Silicate', charthick=ct

; Close PS device

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

end
