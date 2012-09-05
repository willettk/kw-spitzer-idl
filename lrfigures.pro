;pro lrfigures, stop=stop

;+
; NAME:
;       
;	LRFIGURES
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
;       Written by K. Willett                Jun 08
;	Added all OHMs and control sample; correct spacing of PAH and silicate lines - Oct 08
;	Added H2O water ice feature - Oct 08
; 	Added alternate (common) titles for certain objects, re HWWS suggestion - Dec 08
;	Added individual PAH bars - Mar 09
; 	Switched to 4x1 plots, MULTIPLOT - Mar 2010
;	Modified for 51 OHMs in total - Jun 2010
;	Removed final 1x3 plot from FOR loop; added ygap=0.005; switched to log(flux) - Aug 2010
;-

multiplot, /reset

; Retreive all OHM lo-res spectra, sorted by IRAS number

ohmtag = ohmdat('tag')
archtag = archdat('tag')

ohmobj = ohmdat('obj')
archobj = archdat('obj')

allobj = [transpose(ohmobj),transpose(archobj)]
alltag = [transpose(ohmtag),transpose(archtag)]

objsort = sort(allobj)

allobj = allobj[objsort]
alltag = alltag[objsort]

;;;;;;;;;;
;; OHMs ;;
;;;;;;;;;;

iceobjects = ['mega'+string([4,8,17,23,27,28,33],format='(i03)'),$
	'arch'+string([4,5,7,8,10,13,17,18,24,26,29,48,30,31,32,33,35,36,39,40],format='(i03)')]

!p.font=-1

lthick = 5
titlesize = 1.5
labelsize = 1.1
cthick = 5

for m = 0,11 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lrspectra'+strtrim(m+1,2)+'.ps', xs=22, ys=24, /portrait, /encap

	multiplot, [1,4], $
		ytickformat='(f4.1)', $
		ygap = 0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'log S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick

	for j = 4*m,4*m + 3 do begin
	
		tag, alltag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+alltag[j]+'.sav'
		
		flux = sed.flux_lr
		wave = sed.wave_lr
		order = sed.order_lr
		
		noneg = where(sed.flux_lr gt 0)
		flux =   flux[noneg]
		wave =   wave[noneg]
		order = order[noneg]

		flux = alog10(flux)

		; Identified IR lines for overplotting
		
		templines = ir_lines(/lr)
		lines = templines[*,0] & line_id = templines[*,1]
		
		altlist = 'IRAS '+['01418+1651','09320+6134','12540+5708','13428+5608','15327+2340']
		alt_title =['III Zw 35','UGC 5101','Mrk 231','Mrk 273','Arp 220'] 
		title = sed.obj
		tt = where(title eq altlist,count)
		if count eq 1 then title = sed.obj+' ('+alt_title[tt]+')'

		; Plot data
		
		nobonus = where(order ne 3)

		if alltag[j] ne 'mega034' then begin

			plot, wave, flux, $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
				charsize = titlesize, $
				thick = lthick, $
;				yr = [10.^(floor(alog10(min(flux)))),0.9*10^ceil(alog10(max(flux)))], /ystyle, $
				yr = [floor(min(flux)),0.9*ceil(max(flux))], /ystyle, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
		
			oplot, wave[nobonus],flux[nobonus], psym = 10, thick = lthick 

		endif else begin

 			wavenb = wave[nobonus]
			fluxnb = flux[nobonus]
			llmods = where(wavenb gt 12.49)

			plot, wavenb[llmods], fluxnb[llmods], $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
				charsize = titlesize, $
;				yr = [10.^(floor(alog10(min(fluxnb[llmods])))),0.9*10^ceil(alog10(max(fluxnb[llmods])))], /ystyle, $
				yr = [floor(min(flux)),0.9*ceil(max(flux))], /ystyle, $
				thick = lthick, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
	
			oplot, wavenb[llmods], fluxnb[llmods], psym = 10, thick = lthick

		endelse

;		xyouts, 4.4, 10d^(!y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])), title, charsize = 1.0, charthick = lthick
		xyouts, 4.4, !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0]), title, charsize = 1.0, charthick = lthick

		; Instead of vertical lines, create bars with PAH positions. 17.1 um features are likely to
		; 	be the H2 S(1) line, rather than a PAH feature.
	
		case j of
			 0:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch045  17.1 
			 1:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch003  17.1 
			 2:  pahlines = [5.3,        6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;mega001  17.1,
			 3:  pahlines = [5.3,        6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega002  17.1 
			 4:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;arch004  17.1 
			 5:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega004  17.1 
			 6:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch005  17.1 
			 7:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega005  17.1 
			 8:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;mega006  17.1,
			 9:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega007       
			10:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4,    17.4]   ;mega008  17.1,
			11:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;mega009  17.1 
			12:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega010       
			13:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;arch007  17.1 
			14:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch009  17.1,
			15:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega013  17.1 
			16:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch010  17.1 
			17:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch012  17.1 
			18:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega014       
			19:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch013       
			20:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch014  17.1 
			21:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega016  17.1 
			22:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega017  17.1 
			23:  pahlines = [5.3,        6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;mega018  17.1,
			24:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch016       
			25:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch017  17.1 
			26:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;arch018  17.1 
			27:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch020       
			28:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch023       
			29:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch024  17.1,
			30:  pahlines = [5.3,        6.2,  7.7,        11.3,  12.7                       ]   ;arch025       
			31:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega020  17.1 
			32:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch026  17.1,
			33:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;mega022       
			34:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch029  17.1 
			35:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch048       
			36:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega023       
			37:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch030       
			38:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch032  17.1,
			39:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega026  17.1 
			40:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega027       
			41:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;arch033       
			42:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega028  17.1 
			43:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega029  17.1 
			44:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;arch035  17.1 
			45:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega031  17.1 
			46:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega032  17.1 
			47:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch036  17.1,
			48:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega034       
			49:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch039       
			50:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch040  17.1,
		endcase

;		bheight = 10d^(!y.crange[0]+0.15*(!y.crange[1] - !y.crange[0]))
;		blevel  = 10d^(!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))
;		blabel  = 10d^(!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
		bheight = !y.crange[0]+0.15*(!y.crange[1] - !y.crange[0])
		blevel  = !y.crange[0]+0.11*(!y.crange[1] - !y.crange[0])
		blabel  = !y.crange[0]+0.04*(!y.crange[1] - !y.crange[0])
	
		plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
		for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
		xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=labelsize, charthick = lthick

		; Label the silicate features
	
;		sillevel  = 10d^(!y.crange[0]+0.7*(!y.crange(1) - !y.crange(0)))
;		sillevel2 = 10d^(!y.crange[0]+0.75*(!y.crange(1) - !y.crange(0)))
;		sillabel  = 10d^(!y.crange[0]+0.85*(!y.crange(1) - !y.crange(0)))
		sillevel  = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
		sillevel2 = !y.crange[0]+0.75*(!y.crange[1] - !y.crange[0])
		sillabel  = !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])
		plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
		xyouts, 12, sillabel, 'Silicate', charsize = labelsize, charthick = lthick
	
		; Label water ice feature

		junk = where(alltag[j] eq iceobjects,jcount)
		if jcount eq 1 then begin
;			h2olevel = 10d^(!y.crange[0]+0.65*(!y.crange[1] - !y.crange[0]))
;			h2olabel = 10d^(!y.crange[0]+0.70*(!y.crange[1] - !y.crange[0]))
			h2olevel = !y.crange[0]+0.65*(!y.crange[1] - !y.crange[0])
			h2olabel = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
			plots, [5.5,6.2], [h2olevel,h2olevel], linestyle=0, thick=lthick
			xyouts, 5.1, h2olabel, 'H!I2!NO ice', charsize = labelsize, charthick = lthick
		endif
	
		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

endfor

; Last plot - lrspectra13.ps (make it 1x3)

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lrspectra'+strtrim(13,2)+'.ps', xs=22, ys=24, /portrait, /encap

	multiplot, [1,3], $
		ytickformat='(f4.1)', $
		ygap = 0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'log S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick

	for j = 48, 50 do begin
	
		tag, alltag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+alltag[j]+'.sav'
		
		flux = sed.flux_lr
		wave = sed.wave_lr
		order = sed.order_lr
		
		noneg = where(sed.flux_lr gt 0)
		flux =   flux[noneg]
		wave =   wave[noneg]
		order = order[noneg]

		flux = alog10(flux)

		; Identified IR lines for overplotting
		
		templines = ir_lines(/lr)
		lines = templines[*,0] & line_id = templines[*,1]
		
		altlist = 'IRAS '+['01418+1651','09320+6134','12540+5708','13428+5608','15327+2340']
		alt_title =['III Zw 35','UGC 5101','Mrk 231','Mrk 273','Arp 220'] 
		title = sed.obj
		tt = where(title eq altlist,count)
		if count eq 1 then title = sed.obj+' ('+alt_title[tt]+')'

		; Plot data
		
		nobonus = where(order ne 3)

		if alltag[j] ne 'mega034' then begin

			plot, wave, flux, $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
				charsize = titlesize, $
				thick = lthick, $
;				yr = [10.^(floor(alog10(min(flux)))),0.9*10^ceil(alog10(max(flux)))], /ystyle, $
				yr = [floor(min(flux)),0.9*ceil(max(flux))], /ystyle, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
		
			oplot, wave[nobonus],flux[nobonus], psym = 10, thick = lthick 

		endif else begin

 			wavenb = wave[nobonus]
			fluxnb = flux[nobonus]
			llmods = where(wavenb gt 12.49)

			plot, wavenb[llmods], fluxnb[llmods], $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
				charsize = titlesize, $
;				yr = [10.^(floor(alog10(min(fluxnb[llmods])))),0.9*10^ceil(alog10(max(fluxnb[llmods])))], /ystyle, $
				yr = [floor(min(fluxnb[llmods])),0.9*ceil(max(fluxnb[llmods]))], /ystyle, $
				thick = lthick, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
	
			oplot, wavenb[llmods], fluxnb[llmods], psym = 10, thick = lthick

		endelse

		xyouts, 4.4, !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0]), title, charsize = 1.0, charthick = lthick

		; Instead of vertical lines, create bars with PAH positions. 17.1 um features are likely to
		; 	be the H2 S(1) line, rather than a PAH feature.
	
		case j of
			48:  pahlines = [                                     12.7                       ]   ;mega034       
			49:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch039       
			50:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch040  17.1,
		endcase

		bheight = !y.crange[0]+0.15*(!y.crange[1] - !y.crange[0])
		blevel  = !y.crange[0]+0.11*(!y.crange[1] - !y.crange[0])
		blabel  = !y.crange[0]+0.04*(!y.crange[1] - !y.crange[0])
	
		plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
		for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
		xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=labelsize, charthick = lthick

		; Label the silicate features
	
		sillevel  = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
		sillevel2 = !y.crange[0]+0.75*(!y.crange[1] - !y.crange[0])
		sillabel  = !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])
		plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
		xyouts, 12, sillabel, 'Silicate', charsize = labelsize, charthick = lthick
	
		; Label water ice feature

		junk = where(alltag[j] eq iceobjects,jcount)
		if jcount eq 1 then begin
			h2olevel = !y.crange[0]+0.65*(!y.crange[1] - !y.crange[0])
			h2olabel = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
			plots, [5.5,6.2], [h2olevel,h2olevel], linestyle=0, thick=lthick
			xyouts, 5.1, h2olabel, 'H!I2!NO ice', charsize = labelsize, charthick = lthick
		endif
	
		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

;;;;;;;;;;;;;;;;;;;;
;; Control sample ;;
;;;;;;;;;;;;;;;;;;;;

contag = condat('tag')
conobj = condat('obj')

consort = sort(conobj)

contag = contag[consort]
conobj = conobj[consort]

iceobjects = ['control'+string([4,35,34],format='(i03)')]

for m = 0,2 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lrspectra_con'+strtrim(m+1,2)+'.ps', xs=22, ys=24, /portrait, /encap

	multiplot, [1,4], $
		ytickformat='(f4.1)', $
		ygap = 0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'log S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick

	for j = 4*m,4*m + 3 do begin
	
		tag, contag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+contag[j]+'.sav'
		
		flux = sed.flux_lr
		wave = sed.wave_lr
		order = sed.order_lr

			noneg = where(sed.flux_lr gt 0)
			flux =   flux[noneg]
			wave =   wave[noneg]
			order = order[noneg]

		flux = alog10(flux)

		; Identified IR lines for overplotting
		
		templines = ir_lines(/lr)
		lines = templines[*,0] & line_id = templines[*,1]
		
		; Plot data
		
		if contag[j] ne 'control033' then begin

			nobonus = where(order ne 3)
		
			title = sed.obj
			if title eq 'IRAS 01572+0009' then title = sed.obj+' (Mrk 1014)'

			plot, wave, flux, $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
;				yr = [10.^(floor(alog10(min(flux)))),0.9*10^ceil(alog10(max(flux)))], /ystyle, $
				yr = [floor(min(flux)),0.9*ceil(max(flux))], /ystyle, $
				charsize = titlesize, $
				thick = lthick, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
		
				oplot, wave[nobonus],flux[nobonus], psym = 10, thick = lthick 

		endif else begin

 			wavenb = wave[nobonus]
			fluxnb = flux[nobonus]
			llmods = where(wavenb gt 11.90)

			title = sed.obj

			plot, wavenb[llmods], fluxnb[llmods], $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
				charsize = titlesize, $
;				yr = [10.^(floor(alog10(min(fluxnb[llmods])))),0.9*10^ceil(alog10(max(fluxnb[llmods])))], /ystyle, $
				yr = [floor(min(fluxnb[llmods])),0.9*ceil(max(fluxnb[llmods]))], /ystyle, $
				thick = lthick, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
	
			oplot, wavenb[llmods], fluxnb[llmods], psym = 10, thick = lthick

		endelse
	
;		xyouts, 4.4, 10d^(!y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])), title, charsize = 1.0, charthick = lthick
		xyouts, 4.4, !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0]), title, charsize = 1.0, charthick = lthick

		; Instead of vertical lines, create bars with PAH positions
	
		case j of
			0:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control023  17.1, 
			1:   pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2,          17.4]  ;control037        
			2:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control026  17.1, 
			3:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control024  17.1, 
			4:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control004  17.1, 
			5:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4,   17.4]  ;control025  17.1, 
			6:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control028  17.1, 
			7:   pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,         14.2               ]  ;control035        
			8:   pahlines = [0   							        ]  ;control008        
			9:   pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control034  17.1, 
			10:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4,   17.4]  ;control013  17.1, 
			11:  pahlines = [                                     12.7                      ]  ;control033        
		endcase 

;		bheight = 10d^(!y.crange[0]+0.15*(!y.crange[1] - !y.crange[0]))
;		blevel  = 10d^(!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))
;		blabel  = 10d^(!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
		bheight = !y.crange[0]+0.15*(!y.crange[1] - !y.crange[0])
		blevel  = !y.crange[0]+0.11*(!y.crange[1] - !y.crange[0])
		blabel  = !y.crange[0]+0.04*(!y.crange[1] - !y.crange[0])
	
		if j ne 8 then begin
			plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
			for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
			xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=labelsize, charthick = lthick
		endif

		; Label the silicate features
	
;		sillevel  = 10d^(!y.crange[0]+0.7*(!y.crange(1) - !y.crange(0)))
;		sillevel2 = 10d^(!y.crange[0]+0.75*(!y.crange(1) - !y.crange(0)))
;		sillabel  = 10d^(!y.crange[0]+0.85*(!y.crange(1) - !y.crange(0)))
		sillevel  = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
		sillevel2 = !y.crange[0]+0.75*(!y.crange[1] - !y.crange[0])
		sillabel  = !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])
		if contag[j] ne 'control033' then plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
		xyouts, 12, sillabel, 'Silicate', charsize = labelsize, charthick = lthick
	
		; Label water ice feature

		junk = where(contag[j] eq iceobjects,jcount)
		if jcount eq 1 then begin
;			h2olevel = 10d^(!y.crange[0]+0.65*(!y.crange[1] - !y.crange[0]))
;			h2olabel = 10d^(!y.crange[0]+0.70*(!y.crange[1] - !y.crange[0]))
			h2olevel = !y.crange[0]+0.65*(!y.crange[1] - !y.crange[0])
			h2olabel = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
			plots, [5.5,6.2], [h2olevel,h2olevel], linestyle=0, thick=lthick
			xyouts, 5.1, h2olabel, 'H!I2!NO ice', charsize = labelsize, charthick = lthick
		endif

		multiplot
		if j eq 14 then j = 15
	

	endfor	
	
	device,/close
	set_plot,'x'

endfor

; Last plot - make it 1x3

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lrspectra_con'+strtrim(4,2)+'.ps', xs=22, ys=24, /portrait, /encap

	multiplot, [1,3], $
		ytickformat='(f4.1)', $
		ygap = 0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'log S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick

	for j = 12,14 do begin
	
		tag, contag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+contag[j]+'.sav'
		
		flux = sed.flux_lr
		wave = sed.wave_lr
		order = sed.order_lr

			noneg = where(sed.flux_lr gt 0)
			flux =   flux[noneg]
			wave =   wave[noneg]
			order = order[noneg]

		flux = alog10(flux)

		; Identified IR lines for overplotting
		
		templines = ir_lines(/lr)
		lines = templines[*,0] & line_id = templines[*,1]
		
		; Plot data
		
		if contag[j] ne 'control033' then begin

			nobonus = where(order ne 3)
		
			title = sed.obj
			if title eq 'IRAS 01572+0009' then title = sed.obj+' (Mrk 1014)'

			plot, wave, flux, $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
;				yr = [10.^(floor(alog10(min(flux)))),0.9*10^ceil(alog10(max(flux)))], /ystyle, $
				yr = [floor(min(flux)),0.9*ceil(max(flux))], /ystyle, $
				charsize = titlesize, $
				thick = lthick, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
		
				oplot, wave[nobonus],flux[nobonus], psym = 10, thick = lthick 

		endif else begin

 			wavenb = wave[nobonus]
			fluxnb = flux[nobonus]
			llmods = where(wavenb gt 11.90)

			title = sed.obj

			plot, wavenb[llmods], fluxnb[llmods], $
				/xlog, $
;				/ylog, $
				xticks = 6, $
				xtickv = [5,10,15,20,25,30,35], $
				xrange = [4,36], /xstyle, $
				charsize = titlesize, $
;				yr = [10.^(floor(alog10(min(fluxnb[llmods])))),0.9*10^ceil(alog10(max(fluxnb[llmods])))], /ystyle, $
				yr = [floor(min(fluxnb[llmods])),0.9*ceil(max(fluxnb[llmods]))], /ystyle, $
				thick = lthick, $
				charthick = lthick, $
				xthick = lthick, $
				ythick = lthick, $
				/nodata
	
			oplot, wavenb[llmods], fluxnb[llmods], psym = 10, thick = lthick

		endelse
	
;		xyouts, 4.4, 10d^(!y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])), title, charsize = 1.0, charthick = lthick
		xyouts, 4.4, !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0]), title, charsize = 1.0, charthick = lthick

		; Instead of vertical lines, create bars with PAH positions
	
		case j of
			12:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,          17.4]  ;control039  17.1, 
			13:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,   17.4]  ;control040  17.1, 
			14:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2               ]  ;control036  17.1  
		endcase 

		bheight = !y.crange[0]+0.15*(!y.crange[1] - !y.crange[0])
		blevel  = !y.crange[0]+0.11*(!y.crange[1] - !y.crange[0])
		blabel  = !y.crange[0]+0.04*(!y.crange[1] - !y.crange[0])
	
		plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
		for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
		xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=labelsize, charthick = lthick

		; Label the silicate features
	
		sillevel  = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
		sillevel2 = !y.crange[0]+0.75*(!y.crange[1] - !y.crange[0])
		sillabel  = !y.crange[0]+0.85*(!y.crange[1] - !y.crange[0])
		if contag[j] ne 'control033' then plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
		xyouts, 12, sillabel, 'Silicate', charsize = labelsize, charthick = lthick
	
		; Label water ice feature

		junk = where(contag[j] eq iceobjects,jcount)
		if jcount eq 1 then begin
			h2olevel = !y.crange[0]+0.65*(!y.crange[1] - !y.crange[0])
			h2olabel = !y.crange[0]+0.70*(!y.crange[1] - !y.crange[0])
			plots, [5.5,6.2], [h2olevel,h2olevel], linestyle=0, thick=lthick
			xyouts, 5.1, h2olabel, 'H!I2!NO ice', charsize = labelsize, charthick = lthick
		endif

		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

if keyword_set(stop) then stop

end
