pro lhcal3, fnames, nolines = nolines, ps = ps

;+
; NAME:
;       LHCAL3
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
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;
; EXAMPLE:
;	IDL> lhcal3, 'mega023'
;
; REVISION HISTORY
;       Written by K. Willett                Sep 2007
;-

nk = n_elements(fnames)

!p.multi=[0,1,4]

plotdir = '~/Astronomy/Comps2/figures/'
plotname = plotdir+strjoin(fnames)+'_lh.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, xs = 18, ys = 24, yoff = 2, /portrait
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

targets, fnames(k), redshift, obj
tag, fnames(k), dirtag

; Read in the optimally extracted data from saved structure

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fnames(k)+'.sav'

wave = sed.wave_lh
flux = sed.flux_lh

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

linelist = [5.511, 6.1088, 6.2, 6.909, 6.98, 7.7, 8.025, 8.6, 8.99138, 9.665, $
	 10.511, 11.3, 12.279, 12.6, 12.814, $
	  14.2,  14.322, 15.555, 16.4, 17.035, 17.4, $
	18.713,  24.318, 25.890, 28.218]

line_id = ['H!I2!N S(7)', 'H!I2!N S(6)', 'PAH', 'H!I2!N S(5)', 'ArII', 'PAH', 'H!I2!N S(4)', 'PAH', 'ArIII', 'H!I2!N S(3)', $
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

; Plot the data in rest frame wavelengths 

xr = [min(wave)-0.5,max(wave)+0.5]
yr = [-0.1,max(flux)+0.15]

plot, wave, flux, $
	xrange = xr, /xstyle, $
	yrange = yr, /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = cs, $
	thick = lthick, $
	charthick = cthick, $
	/nodata

oplot, wave, flux, psym = 10, thick = lthick

; Overplot vertical lines with labels indicating locations of likely absorption and emission features

if not keyword_set(nolines) then begin
	for i = 0, n_elements(linelist) - 1 do begin
		off = i mod 2
		plots, [linelist(i),linelist(i)],$
			[flux(closeto(wave,linelist(i)))+0.03, flux(closeto(wave,linelist(i)))+0.1], linestyle = 1
		xyouts, linelist(i), flux(closeto(wave,linelist(i)))+0.1 + 0.05*off, line_id(i), $
			orientation = 90, charsize = labsize, /data, charthick = lthick
	endfor
endif

endfor

; Hard copy

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

end

