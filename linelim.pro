function linelim, obj, ion, $
	linecen = linecen, lh = lh, width = width, $
	nolines = nolines, nsig = nsig, lores = lores, $
	noplot = noplot, string=string, quiet = quiet
;+
; NAME:
;       
;	LINELIM.pro
;
; PURPOSE:
;
;	Compute upper limit on line fluxes measured from HR IRS data
;
; INPUTS:
;
;	OBJ - 		string giving object tag (eg, 'mega001')
;
;	ION - 		string giving ion to match (eg, 'neII')
;
; OUTPUTS:
;
;	
;
; KEYWORDS:
;
;	LH - 		set to load data for LH spectra (default is for SH)
;
;	LORES - 	set to load data for LR spectra (default is for SH)
;
;	LINECEN - 	if ION is not specified, a float value giving location at which to measure limit
;
;	NOLINES - 	do not show vertical lines indicating transitions in the plot
;
;	NOPLOT - 	do not plot results to screen
;
;	NSIG - 		significance (in std. deviations) to which to compute the limit. Default is 3. 
;
;	WIDTH - 	float giving distance (in um) around line center to which to fit a baseline
;
; EXAMPLE:
;
;	IDL> result = linelim('mega001','neII')
;
; NOTES:
;
;	It would be nice if the program were adapted so that instead of an ion, one could
;		give a rest wavelength float and the program would compute the limit based on some
;		supplied FWHM. 
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;	Added NOPLOT, LORES keywords - Aug 08
;	For CSOs without any detections, try using the OHM ULIRGs to determine an average FWHM for the flux limit - Aug 09
;-

; Read in data

tagid = strmid(obj,0,3)
tag,tagid,dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'

if n_elements(ion) ne 0 then begin
	linecen = hrwavelength(ion)
	linecen = linecen(0)
endif

if not keyword_set(nsig) then nsig = 3.

if keyword_set(lores) then begin
	wave = sed.wave_lr
	flux = sed.flux_lr
	err  = sed.err_lr
	module = 'LR'
endif else if keyword_set(lh) then begin
	wave = sed.wave_lh
	flux = sed.flux_lh
	err  = sed.err_lh
	module = 'LH'
endif else begin
	wave = sed.wave_sh
	flux = sed.flux_sh
	err  = sed.err_sh
	module = 'SH'
endelse

if not keyword_set(width) then begin
	if not keyword_set(lores) then width = 0.2 else width = 1.0
endif

bind = closeto(wave,linecen-width)
eind = closeto(wave,linecen+width)

; Assuming ION is specified, find ensemble FWHM for measured lines

if dirtag eq 'CSO' then begin					; For CSOs without ANY detections, try the OHM ULIRG sample for the FWHM

	cso_line_exist = file_test('~/Astronomy/Research/Spitzer/cso/linedata/'+ion+'.sav')
	if cso_line_exist eq 1 then restore,'~/Astronomy/Research/Spitzer/cso/linedata/'+ion+'.sav' $
		else restore,'~/Astronomy/Research/Spitzer/linedata/'+ion+'.sav' 

endif else restore,'~/Astronomy/Research/Spitzer/linedata/'+ion+'.sav'

if n_elements(line.fwhm) gt 1 then avg_fwhm =   mean(line.fwhm) else avg_fwhm = line.fwhm
if n_elements(line.fwhm) gt 1 then sig_fwhm = stddev(line.fwhm) else sig_fwhm = line.fwhm

if bind eq eind then begin
	if not keyword_set(quiet) then message,strupcase(ion)+' line does not lie within spectral range of '+module+' module', /info
	area = -1
	return, area
endif

; Linear fit to continuum on either side of the expected line location

contindices = [bind,closeto(wave,linecen - avg_fwhm / 2),closeto(wave,linecen + avg_fwhm / 2), eind]	
contrange = fix([fillarr(1,contindices(0),contindices(1)),fillarr(1,contindices(2),contindices(3))])
expr = 'p[0] + x*p[1]'
start = [1d-2,1d-2]
result = mpfitexpr(expr,wave(contrange),flux(contrange),err(contrange),start,/quiet)

fit = result(0) + wave*result(1) ;+ wave^2 * result(2)
polyfit = flux - fit

meanflux = mean(polyfit(bind:eind))
sigflux  = stddev(polyfit(bind:eind))
	
; Plot spectrum and flux limits

if not keyword_set(noplot) then begin

!p.multi=[0,1,2]
	
	plot, wave, flux, $
		psym=10, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		xr=[linecen-1,linecen+1],/xstyle, $
		title='IRS spectrum with baseline fit'
	
	oplot, wave, fit, color=fsc_color("Cyan"), thick=2
	
	ver, linecen, color = fsc_color("Red")
	
	; Overplot vertical lines with labels indicating locations of likely absorption and emission features
	
	templines = ir_lines(/hr)
	lines = templines(*,0) & line_id = templines(*,1)
	lineindices = where(lines gt !x.crange(0) and lines lt !x.crange(1))
	
	if lineindices[0] ne -1 then begin
	
		lines = lines(lineindices) & line_id = line_id(lineindices)
		
		if not keyword_set(nolines) then begin
			for i = 0, n_elements(lines) - 1 do begin
				off = i mod 2
				ver, lines(i), linestyle = 1, color = defcolor
				xyouts, lines(i), 0.65*!y.crange(1) + 0.05*off*!y.crange(1), line_id(i), $
					orientation = 90, charsize = 1.5, /data, $
					color = defcolor, charthick = cthick
			endfor
		endif
	endif
	
	ver,wave(bind), linestyle = 2
	ver,wave(eind), linestyle = 2
	
	plot, wave, polyfit, $
		psym=10, $
		xtitle = 'Wavelength [!7l!3m]', $
		ytitle = 'Flux [Jy]', $
		title = 'Baseline subtracted', $
		xr=[linecen-1,linecen+1],/xstyle;, $
	;	yr=[0.5*lowersig,1.5*uppersig],/ystyle
	
	ver, linecen, color = fsc_color("Red")
	
	ver,wave(bind), linestyle = 2
	ver,wave(eind), linestyle = 2
	
	hor, nsig * sigflux, color=fsc_color("Yellow")
	hor, -1d * nsig * sigflux, color=fsc_color("Yellow")
		
endif

;print,'Sig. flux = ',sigflux
;print,'Avg. FWHM = ',avg_fwhm
;print,'Sig. FWHM = ',sig_fwhm

c = 299792.458d * 1d9
flux_cgs = nsig * sigflux / linecen^2 * c * 1d-30
sigwidth = avg_fwhm / (2d * sqrt(2d * alog(2d)))

;print,'Flux [cgs]',flux_cgs

area = flux_cgs * sigwidth * sqrt(2d * !dpi) / 1d-21		; Units of 10^-21 W/cm^2

wavearr = fillarr(1d-3,min(wave),max(wave))
fakeline = nsig * sigflux * exp(-1d * (wavearr - linecen)^2 / (2d * sig_fwhm)^2)

if not keyword_set(noplot) then begin
	oplot, wavearr, fakeline, thick=1.5, color=fsc_color("Green")
	!p.multi=[0,1,1]
endif

if keyword_set(string) then area = string(area,format='(f7.2)')

return, area

end
