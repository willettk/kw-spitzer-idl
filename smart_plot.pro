pro smart_plot, filename

;+
; NAME:
;       SMART_PLOT
;
; PURPOSE:
;       Plot LH data written out as an ASCII file after extraction by SMART
;
; CALLING SEQUENCE:
;       smart_plot
;
; OUTPUTS:
;
; EXAMPLE:
;	IDL> smart_plot
;
;
; REVISION HISTORY
;       Written by K. Willett                Feb 2007
;-

readcol, filename, wave, flux, stdev, det, line, sdir, scnt, status, sky, skye, flag, $
	format = 'f,f,f,i,i,i,i,i,f,f,l', /silent

colors = [fsc_color('Red'), $
	fsc_color('Orange'), $
	fsc_color('Yellow'), $
	fsc_color('Green'), $
	fsc_color('Blue'), $
	fsc_color('Purple'), $
	fsc_color('White'), $
	fsc_color('Hot Pink'), $
	fsc_color('Dodger Blue'), $
	fsc_color('Goldenrod')]

plot, wave, flux, /nodata, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
	title = 'SMART data plotted by order'

for i = 0, 9 do begin
	a = where(det eq 11 + i)
	c = sort(wave(a))
	oplot, wave(c+min(a)), flux(c+min(a)), color = colors(i)
endfor

stop
end
