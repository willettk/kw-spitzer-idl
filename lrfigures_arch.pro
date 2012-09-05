
;+
; NAME:
;       
;	LRFIGURES_ARCH
;
; PURPOSE:
;
;	Create figures displaying IRS LR spectra for paper
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
;       Written by K. Willett                Aug 08
;-


archtag = archdat('tag')
narch = n_elements(archtag)

!p.multi=[0,2,4]

for m = 0,3 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lrspectra_arch'+strtrim(m+1,2)+'.ps', $
		xs=18, ys=24,/portrait,yoff=2

	for j = 8*m,8*m + 7 do begin
	
		tag, archtag(j), dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+archtag(j)+'.sav'
		
		flux = sed.flux_lr
		wave = sed.wave_lr
		order = sed.order_lr
		
			noneg = where(sed.flux_lr gt 0)
			flux = flux(noneg)
			wave = wave(noneg)
			order = order(noneg)

		; Identified IR lines for overplotting
		
		templines = ir_lines(/lr)
		lines = templines(*,0) & line_id = templines(*,1)
		
		; Plot data
		
		; Hard copy option
		
		lthick = 2
		cs = 1
		
		nobonus = where(order ne 3)
	
		plot, wave, flux, $
			/xlog, $
			/ylog, $
	;		xticks = 13, $
	;		xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
			xticks = 6, $
			xtickv = [5,10,15,20,25,30,35], $
			xrange = [4,36], /xstyle, $
	;		yrange = [1d-4,max(flux)], $
			xtitle = 'Wavelength (rest frame) [!7l!3m]', $
			ytitle = 'Flux density [Jy]', $
			title = sed.obj, $
			charsize = 1.1, $
			thick = lthick, $
			charthick = lthick, $
			/nodata
	
		oplot, wave(nobonus),flux(nobonus), psym = 10, thick = lthick
		
		qwer=0
		if qwer then begin
		for i = 0, n_elements(lines) - 1 do begin
			off = i mod 2
			ver, lines(i), linestyle = 1
			xyouts, lines(i), 10d^((1.0 -  0.08*off)*alog10(yr[1]/yr[0]) + alog10(yr[0])), line_id(i), $
				orientation = 90, charsize = 0.5, /data
		endfor
		endif
		
		; Instead of vertical lines, create bars with PAH positions
	
		pahlines = [6.2,7.7,8.6,11.3,12.7]
		bheight = 10d^(!y.crange[0]+0.15*(!y.crange[1] - !y.crange[0]))
		blevel  = 10d^(!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))
		blabel  = 10d^(!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
	
		plots, [pahlines(0),pahlines(4)], [blevel,blevel]
		for k=0,4 do plots, [pahlines(k),pahlines(k)], [blevel,bheight]
		xyouts,8.0,blabel,'PAH',charsize=0.5
	
		; Label the silicate features
	
		sillevel = 10d^(!y.crange[0]+0.7*(!y.crange(1) - !y.crange(0)))
		sillevel2 = 10d^(!y.crange[0]+0.75*(!y.crange(1) - !y.crange(0)))
		sillabel = 10d^(!y.crange[0]+0.85*(!y.crange(1) - !y.crange(0)))
		plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
		xyouts, 12, sillabel, 'Silicates', charsize = 0.7
	
		if m eq 3 and j eq 29 then j = 8*m + 7
	endfor	
	
	device,/close
	set_plot,'x'

endfor

end
