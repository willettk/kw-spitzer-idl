;+
; NAME:
;       CSO2_LRPLOT
;
; PURPOSE:
; 	Make quick 2x2 plot of new LR CSOs
;
;
; REVISION HISTORY
;       Adapted from LRSTITCH		- Jun 09
;-

names=['cso007','cso008','cso009','cso010']

!p.multi=[0,2,2]
		set_plot,'ps'
		device, /landscape, filename = plotdir+'cso2_lrplot.ps'
		defcolor = fsc_color('Black')
		lthick = 2
		cs = 1
	
for j=0,3 do begin

	; Load the stitched, calibrated data from the IDL sav files
	
	tag, names[j], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+names[j]+'.sav'
	
	flux = sed.flux_lr
	wave = sed.wave_lr
	order = sed.order_lr
	
	noneg = where(sed.flux_lr gt 0)
	flux = flux(noneg)
	wave = wave(noneg)
	order = order(noneg)
	
	plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibrated/'
	
	; Identified IR lines for overplotting
	
	templines = ir_lines(/lr)
	lines = templines(*,0) & line_id = templines(*,1)
	
	; Plot data
	
	; Hard copy option
	

	case j of
		0: begin
		xr=[4,40] 
		yr=[5d-4,2d-2]
		end
		1: begin
		xr=[4,40] 
		yr=[5d-4,1d-1]
		end
		2: begin
		xr=[4,40] 
		yr=[5d-3,5d-1]
		end
		3: begin
		xr=[4,40] 
		yr=[1d-4,1d-1]
		end
	endcase

	sl_bonus = where(order eq 3 and wave lt 10.)
	ll_bonus = where(order eq 3 and wave gt 10.)
	
	sl2_index = where((order eq 2) and (wave lt 10.))		; An ugly solution. There should be a better way of distinguishing modules.
	sl1_index = where((order eq 1) and (wave lt 15.))
	ll2_index = where((order eq 2) and (wave gt 10.))
	ll1_index = where((order eq 1) and (wave gt 15.))
	
		xlog = 1
		ylog = 1
	
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
		charthick = lthick, $
		/nodata
	
	oplot, wave,flux, psym = 10, thick = lthick
	
;		for i = 0, n_elements(lines) - 1 do begin
;			off = i mod 2
;			ver, lines(i), linestyle = 1, color = defcolor
;			xyouts, lines(i), 10d^((0.90 + 0.10*off)*alog10(yr[1]/yr[0]) + alog10(yr[0])), $
;				line_id(i), orientation = 90, charsize = cs, /data
;		endfor
;			ver,6.85,linestyle=2
;			ver,7.25,linestyle=2
	


endfor

		device, /close
		set_plot,'x'
	
end
