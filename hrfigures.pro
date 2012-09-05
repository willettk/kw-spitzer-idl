;+
; NAME:
;       
;	HRFIGURES
;
; PURPOSE:
;
;	Create figures displaying IRS HR spectra for paper
;
; INPUTS:
;
;
;
; OUTPUTS:
;
;
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Jun 08
;-


ohmtag = ohmdat('tag')
nohm = n_elements(ohmtag)

!p.multi=[0,1,4]

for m = 0,5 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/hrspectra'+strtrim(m+1,2)+'.ps', xs=18, ys=24

	for j = 4*m,4*m + 3 do begin
	
		tag, ohmtag(j), dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+ohmtag(j)+'.sav'
		
		flux_sh = sed.flux_sh
		wave_sh = sed.wave_sh
		order_sh = sed.order_sh
		
		noneg_sh = where(sed.flux_sh gt 0)
		flux_sh = flux_sh(noneg_sh)
		wave_sh = wave_sh(noneg_sh)
		order_sh = order_sh(noneg_sh)
		
		flux_lh = sed.flux_lh
		wave_lh = sed.wave_lh
		order_lh = sed.order_lh
		
		noneg_lh = where(sed.flux_lh gt 0)
		flux_lh = flux_lh(noneg_lh)
		wave_lh = wave_lh(noneg_lh)
		order_lh = order_lh(noneg_lh)
		
		; Identified IR lines for overplotting
		
		templines = ir_lines(/hr)
		lines = templines(*,0) & line_id = templines(*,1)
		
		; Plot data
		
		; Hard copy option
		
		lthick = 2
		cs = 1
		
		plot, wave, flux, $
			xrange = [8,33], /xstyle, $
			yrange = [0,max(flux_lh)], $
			xtitle = 'Wavelength (rest frame) [!7l!3m]', $
			ytitle = 'Flux density [Jy]', $
			title = sed.obj, $
			charsize = 1.1, $
			thick = lthick, $
			charthick = lthick, $
			/nodata
	
		oplot, wave_sh,flux_sh, psym = 10, thick = lthick
		oplot, wave_lh,flux_lh, psym = 10, thick = lthick
		
		qwer=0
		if qwer then begin
		for i = 0, n_elements(lines) - 1 do begin
			off = i mod 2
			ver, lines(i), linestyle = 1
			xyouts, $
				lines(i), $
				10d^((1.0 -  0.08*off)*alog10(!y.crange[1]/!y.crange[0]) + alog10(!y.crange[0])), $
				line_id(i), $
				orientation = 90, charsize = 0.5, /data
		endfor
		endif
		
	
	endfor	
	
	device,/close
	set_plot,'x'

endfor
	

end

