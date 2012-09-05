pro lrcal3, fnames, ps = ps, nolines = nolines, xr = xr, yr = yr
;+
; NAME:
;       LRCAL3
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
;	NOLINES - removes common mid-IR lines with labels. Default is to display them. 
;
;	PS - plots the spectrum to a hardcopy postscript file
;
; EXAMPLE:
;	IDL> lrcal3, ['mega023','mega025','mega026']
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

; Hard copy option

nk = n_elements(fnames)

!p.multi=[0,1,4]

plotdir = '~/Astronomy/Comps2/figures/'
plotname = plotdir+strjoin(fnames)+'_lr.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, xs = 18, ys = 24, yoff = 2,/portrait
	lthick = 2
	cs = 1.5
	labsize = 0.5
endif else begin
	lthick = 1
	cs = 2
	labsize=2
endelse

for k = 0, nk-1 do begin

; Locate directory from which to read (either megamasers or CSOs)

tag, fnames(k), dirtag 
targets, fnames(k), redshift, obj

; Read in the optimally extracted data from saved structure

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fnames(k)+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

linelist = [5.511, 6.2, 6.85, 6.909, 6.98, 7.25, 7.7, 8.025, 8.6, 8.99138, 9.665, $
	 10.511, 11.3, 12.279, 12.6, 12.814, $
	  14.2,  14.322, 15.555, 16.4, 17.035, 17.4, $
	18.713,  24.318, 25.890, 28.218]

line_id = ['H!I2!N S(7)', 'PAH', 'C-H', 'H!I2!N S(5)', 'ArII', 'C-H','PAH', 'H!I2!N S(4)', 'PAH', 'ArIII', 'H!I2!N S(3)', $
	'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
	  'PAH',  'NeV', 'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', $
	'SIII', 'NeV', 'OIV', 'H!I2!N S(0)']

linestart = where(min(wave) gt linelist)
lineterm = where(max(wave) lt linelist)

if linestart(0) ne -1 then begin
       lsind = [linestart(n_elements(linestart)-1)+1,n_elements(linelist)-1]
       linelist = linelist(lsind(0):lsind(1)) 
       line_id  = line_id (lsind(0):lsind(1))
endif
if lineterm(0)  ne -1 then begin
	ltind = [0,lineterm(0)-1]
	linelist = linelist(ltind(0):ltind(1)) 
	line_id = line_id(ltind(0):ltind(1))
endif

; Plot data

if not keyword_set(xr) then xr = [4,35]
if not keyword_set(yr) then yr = [1d-4,max(flux)]
if sed.tag eq 'control008' then yr=[1d-2,10]

plot, wave, flux, $
	/xlog, $
	/ylog, $
	xticks = 13, $
	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
	xrange = xr, /xstyle, $
	yrange = yr, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = cs, $
	thick = lthick, $
	charthick = lthick, $
	/nodata

oplot, wave, flux, psym = 10, thick = lthick

if not keyword_set(nolines) then begin
	for i = 0, n_elements(linelist) - 1 do begin
		off = i mod 2
		plots, [linelist(i),linelist(i)],$
			[flux(closeto(wave,linelist(i)))*3,$
			flux(closeto(wave,linelist(i)))*6],$
			color=defcolor, linestyle = 0
		lineend = flux(closeto(wave,linelist(i)))*2
		xyouts, linelist(i), lineend + 5*off*lineend, line_id(i), $
			orientation = 90, charsize = labsize, /data, charthick = lthick

		; H2O absorption

		plots,[5.5,6.3],[flux(closetomed(wave,flux,6.0))*0.5, flux(closetomed(wave,flux,6.0))*0.5],linestyle=1
		xyouts, 5.7, flux(closetomed(wave,flux,6.0))*0.1,'H!I2!NO', charsize = labsize

		; Silicate absorption at 9.7 and 18.0 um

		plots,[8,12],[flux(closetomed(wave,flux,9.7))*0.5, flux(closetomed(wave,flux,9.7))*0.5],linestyle=1
		xyouts, 9, flux(closetomed(wave,flux,9.7))*0.1,'Silicate', charsize = labsize
		plots,[17,20],[flux(closetomed(wave,flux,18))*0.5, flux(closetomed(wave,flux,18))*0.5], linestyle=1
		xyouts, 17, flux(closetomed(wave,flux,18))*0.1,'Silicate', charsize = labsize

	endfor
endif

endfor

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

end
