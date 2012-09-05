;pro lrfigures_cso
;+
; NAME:
;       
;	LRFIGURES_CSO
;
; PURPOSE:
;
;	Create figures displaying IRS LR spectra for CSO paper
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
;	IDL> .r lrfigures_cso
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;	Removed NGC 5793, 1245+676 - Nov 09
;-


csotag = csodat('tag')

!p.multi=[0,2,4]

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/cso/papers/lrspectra_cso.ps', xs=18, ys=24,/portrait

	; yrange

	yr = [$
		[4d-3,1d-1], $		; cso001
		[5d-4,5d-2], $		; cso002
		[1d-2,5d-0], $		; cso003
		[3d-2,1d-0], $		; cso004
		[7d-3,3d-1], $		; cso005
		[8d-3,3d-0], $		; cso006
		[3d-4,4d-2], $		; cso007
		[3d-4,1d-1], $		; cso008
		[3d-3,5d-1], $		; cso009
		[2d-5,5d-2]]		; cso010

	; Is there water ice?

	icearr = [0,1,1,0,0,0,0,0,0,0]

	; Sort by RA

	;bind = [1,7,0,9,5,3,4,2,8,6]
	bind = [1,7,0,5,3,4,8,6]
	
	csotag = csotag[bind]
	yr = yr[*,bind]
	icearr = icearr[bind]
	ncso = n_elements(csotag)

	; Begin loop

	for j = 0,ncso - 1 do begin
	
		tag, csotag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+csotag[j]+'.sav'
		
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
		
		lthick = 4
		cthick = 4
		labelthick = 3
		cs = 1
		
		nobonus = where(order ne 3)
	
		plot, wave, flux, $
			/xlog, $
			/ylog, $
	;		xticks = 13, $
	;		xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
			xticks = 7, $
			xtickv = [5,10,15,20,25,30,35,40], $
			xrange = [4,40], /xstyle, $
			yrange = yr[*,j], /ystyle, $
;			yrange = [1d-5,5d0], /ystyle, $
			xtitle = 'Wavelength (rest frame) [!7l!3m]', $
			ytitle = 'Flux density [Jy]', $
			title = sed.obj, $
			charsize = 1.5, $
			thick = lthick, $
			charthick = cthick, $
			xthick = lthick, $
			ythick = lthick, $
			/nodata
	
		oplot, wave(nobonus),flux(nobonus), psym = 10, thick = 3
		
		; Instead of vertical lines, create bars with PAH positions
	
		;case j of
		;	0: pahlines =             [7.7,8.6,11.3,12.7,     16.4,17.1]
		;	1: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
		;	2: pahlines = [5.3,5.7,6.2,7.7,8.6,11.3,12.7,14.2,16.4,     17.4]
		;	3: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
		;	4: pahlines = [0]
		;	5: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2]
		;	6: pahlines = [            7.7,8.6,11.3,12.7,14.2               ]
		;	7: pahlines = [            7.7,    11.3,12.7                    ]
		;	8: pahlines = [    5.7,6.2,7.7,8.6,11.3,12.7,     16.4,     17.4]
		;	9: pahlines = [0]
		;endcase

		;case j of
		;	0: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
		;	1: pahlines = [            7.7,    11.3,12.7                    ]
		;	2: pahlines =             [7.7,8.6,11.3,12.7,     16.4,17.1]
		;	3: pahlines = [0]
		;	4: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2]
		;	5: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
		;	6: pahlines = [0]
		;	7: pahlines = [5.3,5.7,6.2,7.7,8.6,11.3,12.7,14.2,16.4,     17.4]
		;	8: pahlines = [    5.7,6.2,7.7,8.6,11.3,12.7,     16.4,     17.4]
		;	9: pahlines = [            7.7,8.6,11.3,12.7,14.2               ]
		;endcase

		case j of
			0: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
			1: pahlines = [            7.7,    11.3,12.7                    ]
			2: pahlines =             [7.7,8.6,11.3,12.7,     16.4,17.1]
			3: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2]
			4: pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
			5: pahlines = [0]
			6: pahlines = [    5.7,6.2,7.7,8.6,11.3,12.7,     16.4,     17.4]
			7: pahlines = [            7.7,8.6,11.3,12.7,14.2               ]
		endcase

		if csotag[j] ne 'cso005' then begin

			bheight = 10d^(!y.crange[0]+0.15*(!y.crange[1] - !y.crange[0]))
			blevel  = 10d^(!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))
			blabel  = 10d^(!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
		
			plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
			for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
			xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=0.5, charthick = labelthick

		endif
	
		; Label the silicate features
	
		sillevel = 10d^(!y.crange[0]+0.7*(!y.crange(1) - !y.crange(0)))
		sillevel2 = 10d^(!y.crange[0]+0.75*(!y.crange(1) - !y.crange(0)))
		sillabel = 10d^(!y.crange[0]+0.85*(!y.crange(1) - !y.crange(0)))
		plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
		xyouts, 12, sillabel, 'Silicate', charsize = 0.7, charthick = labelthick
	
		; Label water ice feature

		if icearr[j] eq 1 then begin
			h2olevel = 10d^(!y.crange[0]+0.70*(!y.crange[1] - !y.crange[0]))
			h2olabel = 10d^(!y.crange[0]+0.75*(!y.crange[1] - !y.crange[0]))
			plots, [5.5,6.2], [h2olevel,h2olevel], linestyle=0, thick=lthick
			xyouts, 5.1, h2olabel, 'H!I2!NO ice', charsize = 0.7, charthick = labelthick
		endif
	
	endfor	
	
	device,/close
	set_plot,'x'

end
