pro pufilters, ps = ps

; Plot the transmission of IRS red and blue peakup filters as function of wavelength
;
; KW, Jul 07

device, decomposed = 1

fpath = '~/Astronomy/Research/Spitzer/spitzer/'

readcol, fpath+'bluePUtrans.txt', bwave, btrans, /silent

readcol, fpath+'redPUtrans.txt', rwave, rtrans, /silent

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'pufilters.ps', /landscape
	blue = fsc_color('White')
	red = blue
endif else begin
	blue = fsc_color('Blue')
	red = fsc_color('Red')
endelse


plot, bwave, btrans, /nodata, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'QE x gain x transmission', $
	title = 'Response of PU filters for IRS', $
	charsize = 2.0, $
	xrange = [10,35]

oplot, bwave, btrans, color=fsc_color('Blue'), thick = 2, linestyle = 0
oplot, rwave, rtrans, color=fsc_color('Red'),  thick = 2, linestyle = 2

legend, /top, /right, ['Blue PU', 'Red PU'], linestyle = [0,2], color=[fsc_color('Blue'),fsc_color('Red')]

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

end
