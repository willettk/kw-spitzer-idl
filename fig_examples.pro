;+
; NAME:
;       
;	FIG_EXAMPLES
;
; PURPOSE:
;
;	Create figures displaying example IRS LR, SH, LH spectra
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
;       Written by K. Willett                Mar 2010
;-

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

;!p.multi=[0,1,4]

;;;;;;;;;;;;;;;;
;; OHM sample ;;
;;;;;;;;;;;;;;;;


iceobjects = ['mega'+string([4,8,17,23,27,28,33],format='(i03)'),$
	'arch'+string([4,5,7,8,10,13,17,18,24,26,29,48,30,31,32,33,35,36,39,40],format='(i03)')]

set_plot,'ps'
device, filename='~/Astronomy/Research/Spitzer/papers/apj_submission/paperI/fig2.ps', xs=22, ys=24

!p.font=-1

lthick = 5
titlesize = 1.5
labelsize = 1.1
cthick = 5

multiplot, [1,4], $
	mxtitle = 'Rest wavelength [!7l!3m]', $
	mytitle = 'Flux density [Jy]', $
	mxtitsize = titlesize, $
	mytitsize = titlesize, $
	mxtitoffset = 1.2, $
	mytitoffset = 1.2, $
	mcharthick = cthick, $
	mxcharthick = cthick, $
	mycharthick = cthick

for j = 0, 3 do begin

	tag, alltag[j], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+alltag[j]+'.sav'
	
	flux = sed.flux_lr  * 1d3
	wave = sed.wave_lr
	order = sed.order_lr
	
		noneg = where(sed.flux_lr gt 0)
		flux =   flux[noneg]
		wave =   wave[noneg]
		order = order[noneg]

	; Identified IR lines for overplotting
	
	templines = ir_lines(/lr)
	lines = templines[*,0] & line_id = templines[*,1]
	
	; Plot data
	
		nobonus = where(order ne 3)
	
		title = sed.obj
		altlist = 'IRAS '+['01418+1651','09320+6134','12540+5708','13428+5608','15327+2340']
		alt_title =['III Zw 35','UGC 5101','Mrk 231','Mrk 273','Arp 220'] 
		tt = where(title eq altlist,count)
		if count eq 1 then title = sed.obj+' ('+alt_title[tt]+')'

		plot, wave, flux, $
			/xlog, $
			/ylog, $
			xticks = 6, $
			xtickv = [5,10,15,20,25,30,35], $
			xrange = [4,36], /xstyle, $
			yr = [10.^(floor(alog10(min(flux)))),0.9*10^ceil(alog10(max(flux)))], /ystyle, $
			charsize = cs, $
			thick = lthick, $
			charthick = lthick, $
			xthick = lthick, $
			ythick = lthick, $
			/nodata
	
		oplot, wave[nobonus],flux[nobonus], psym = 10, thick = lthick 

		xyouts, 4.4, 10d^(!y.crange[0]+0.8*(!y.crange[1] - !y.crange[0])), title, charsize = 1.0, charthick = lthick

	; Instead of vertical lines, create bars with PAH positions. 17.1 um features are likely to
	; 	be the H2 S(1) line, rather than a PAH feature.

	case j of
		0:   pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch045  17.1 
		1:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch003  17.1 
		2:   pahlines = [5.3,        6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;mega001  17.1,
		3:   pahlines = [5.3,        6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega002  17.1 
		4:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;arch004  17.1 
		5:   pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega004  17.1 
		6:   pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch005  17.1 
		7:   pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega005  17.1 
		8:   pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;mega006  17.1,
		9:   pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega007       
		10:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4,    17.4]   ;mega008  17.1,
		11:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;mega009  17.1 
		12:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega010       
		13:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;arch007  17.1 
		14:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;arch008  17.1 
		15:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch009  17.1,
		16:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega013  17.1 
		17:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch010  17.1 
		18:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch012  17.1 
		19:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega014       
		20:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch013       
		21:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch014  17.1 
		22:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega016  17.1 
		23:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega017  17.1 
		24:  pahlines = [5.3,        6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;mega018  17.1,
		25:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch016       
		26:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch017  17.1 
		27:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;arch018  17.1 
		28:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch020       
		29:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch023       
		30:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch024  17.1,
		31:  pahlines = [5.3,        6.2,  7.7,        11.3,  12.7                       ]   ;arch025       
		32:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega020  17.1 
		33:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch026  17.1,
		34:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,         16.4         ]   ;mega022       
		35:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch029  17.1 
		36:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch048       
		37:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega023       
		38:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch030       
		39:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch031  17.1,
		40:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;arch032  17.1,
		41:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega025  17.1 
		42:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega026  17.1 
		43:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega027       
		44:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;arch033       
		45:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;mega028  17.1 
		46:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch034       
		47:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega029  17.1 
		48:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2                ]   ;arch035  17.1 
		49:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4         ]   ;mega031  17.1 
		50:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega032  17.1 
		51:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch036  17.1,
		52:  pahlines = [5.3,  5.7,  6.2,  7.7,  8.6,  11.3,  12.7,  14.2,  16.4,    17.4]   ;mega033  17.1,
		53:  pahlines = [            6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;mega034       
		54:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7                       ]   ;arch039       
		55:  pahlines = [      5.7,  6.2,  7.7,  8.6,  11.3,  12.7,                  17.4]   ;arch040  17.1,
	endcase

;	pahlines = [6.2,7.7,8.6,11.3,12.7]

	bheight = 10d^(!y.crange[0]+0.15*(!y.crange[1] - !y.crange[0]))
	blevel  = 10d^(!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))
	blabel  = 10d^(!y.crange[0]+0.03*(!y.crange[1] - !y.crange[0]))

	plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
	for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
	xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=labelsize, charthick = lthick

	; Label the silicate features

	sillevel = 10d^(!y.crange[0]+0.7*(!y.crange(1) - !y.crange(0)))
	sillevel2 = 10d^(!y.crange[0]+0.75*(!y.crange(1) - !y.crange(0)))
	sillabel = 10d^(!y.crange[0]+0.85*(!y.crange(1) - !y.crange(0)))
	plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
	plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
	xyouts, 12, sillabel, 'Silicate', charsize = labelsize, charthick = lthick

	; Label water ice feature

	junk = where(alltag[j] eq iceobjects,jcount)
	if jcount eq 1 then begin
		h2olevel = 10d^(!y.crange[0]+0.70*(!y.crange[1] - !y.crange[0]))
		h2olabel = 10d^(!y.crange[0]+0.75*(!y.crange[1] - !y.crange[0]))
		plots, [5.5,6.2], [h2olevel,h2olevel], linestyle=0, thick=lthick
		xyouts, 5.1, h2olabel, 'H!I2!NO ice', charsize = labelsize, charthick = lthick
	endif

	multiplot

endfor	

device,/close
set_plot,'x'

multiplot, /reset

; SH

set_plot,'ps'
device, filename='~/Astronomy/Research/Spitzer/papers/apj_submission/paperI/fig3.ps', xs=22, ys=24

!p.font=-1

lthick = 5
titlesize = 1.5
labelsize = 0.9
cthick = 5

multiplot, [1,4], $
	mxtitle = 'Rest wavelength [!7l!3m]', $
	mytitle = 'Flux density [Jy]', $
	mxtitsize = titlesize, $
	mytitsize = titlesize, $
	mxtitoffset = 1.2, $
	mytitoffset = 1.2, $
	mcharthick = cthick, $
	mxcharthick = cthick, $
	mycharthick = cthick

for j = 0, 3 do begin

	tag, alltag(j), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+alltag(j)+'.sav'
	
	flux_sh = sed.flux_sh
	wave_sh = sed.wave_sh
	order_sh = sed.order_sh
	
	noneg_sh = where(sed.flux_sh gt 0)
	flux_sh = flux_sh(noneg_sh)
	wave_sh = wave_sh(noneg_sh)
	order_sh = order_sh(noneg_sh)
	
	; Identified IR lines for overplotting
	
	templines = ir_lines(/hr)
	lines = templines(*,0) & line_id = templines(*,1)
	
	; Plot data
	
	ymin = min(flux_sh)
	ymax = max(flux_sh) * 1.1

	plot, wave_sh, flux_sh, $
		xrange = [floor(min(wave_sh)),ceil(max(wave_sh))+0.5], /xstyle, $
		yrange = [(ymin - 0.15 * ymax)/0.85,ymax], /ystyle, $
;		xtitle = 'Wavelength (rest frame) [!7l!3m]', $
;		ytitle = 'Flux density [Jy]', $
;		title = title, $
		charsize = 1.1, $
		thick = lthick, $
		charthick = cthick, $
		xthick = lthick, $
		ythick = lthick, $
		/nodata

	oplot, wave_sh,flux_sh, psym = 10, thick = lthick
	
	altlist = 'IRAS '+['01418+1651','09320+6134','12540+5708','13428+5608','15327+2340']
	alt_title =['III Zw 35','UGC 5101','Mrk 231','Mrk 273','Arp 220'] 
	title = sed.obj
	tt = where(title eq altlist,count)
	if count eq 1 then title = sed.obj+' ('+alt_title[tt]+')'

	xyouts, (!x.crange[0]+0.05*(!x.crange[1] - !x.crange[0])), (!y.crange[0]+0.6*(!y.crange[1] - !y.crange[0])), $
		title, charsize = 1.1, charthick = lthick

	; Label detected lines

	thelines = whichlines(alltag[j])
	fcount = n_elements(thelines)

	if fcount gt 0 then begin

		linewaves = dblarr(fcount)
		for i = 0, fcount - 1 do linewaves[i] = hrwavelength(thelines[i])

		; Remove H2

		noh2 = where(strmid(thelines,0,2) ne 'h2',noh2count)
		h2 = where(strmid(thelines,0,2) eq 'h2',h2count)

		if noh2count gt 0 then begin

			fsid = thelines[noh2]
			fslines = linewaves[noh2]
			for i = 0, n_elements(fsid) -1 do fsid[i] = strupcase(strmid(fsid[i],0,1)) + strmid(fsid[i],1,strlen(fsid[i])-1)

			hualpha = where(fsid eq 'HI76',hucount)
			if hucount gt 0 then fsid[hualpha] = 'Hu-!7a!3'

			chlorine = where(fsid eq 'ClII',clcount)
			if clcount gt 0 then fsid[chlorine] = 'Cl II'

			fsind = where(fslines lt max(wave_sh) and fslines gt min(wave_sh))
			fslines = fslines[fsind] & fsid = fsid[fsind]
			bheight = (!y.crange[0] + 0.15*(!y.crange[1] - !y.crange[0]))
			blevel = (!y.crange[0]+0.08*(!y.crange[1] - !y.crange[0]))
			blabellev = (!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))

			plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
			for k=0,n_elements(fslines)-1 do begin
				plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
				if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.4,blabellev,fsid[k],charsize=labelsize, charthick=lthick $
					else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
					then xyouts,fslines[k]-0.3,blabellev,fsid[k],charsize=labelsize, charthick=lthick $
					else xyouts,fslines[k]+0.1,blabellev,fsid[k],charsize=labelsize, charthick=lthick
			endfor

		endif

		if h2count gt 0 then begin

			h2id = 'S('+strmid(thelines[h2],3,4)+')'
			h2lines = linewaves[h2]

			h2ind = where(h2lines lt max(wave_sh),h2count2)
			if h2count2 gt 0 then begin

				h2lines = h2lines[h2ind] & h2id = h2id[h2ind]
				bheight = (!y.crange[1]-0.15*(!y.crange[1] - !y.crange[0]))
				blevel = (!y.crange[0]+0.90*(!y.crange[1] - !y.crange[0]))
				blabel = (!y.crange[0]+0.92*(!y.crange[1] - !y.crange[0]))
				blabellev = (!y.crange[0]+0.83*(!y.crange[1] - !y.crange[0]))
		
				plots, [h2lines(0),h2lines(n_elements(h2lines)-1)], [blevel,blevel], thick = lthick
				for k=0,n_elements(h2lines)-1 do begin
					plots, [h2lines(k),h2lines(k)], [blevel, bheight], thick = lthick
					xyouts,h2lines[k]+0.1,blabellev,h2id[k],charsize=labelsize, charthick=lthick
				endfor
				xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=labelsize, charthick=lthick

			endif

		endif

	endif

	multiplot

endfor	

device,/close
set_plot,'x'

multiplot, /reset

; LH

set_plot,'ps'
device, filename='~/Astronomy/Research/Spitzer/papers/apj_submission/paperI/fig4.ps', xs=22, ys=24

!p.font=-1

lthick = 5
titlesize = 1.5
labelsize = 0.9
cthick = 5

multiplot, [1,4], $
	mxtitle = 'Rest wavelength [!7l!3m]', $
	mytitle = 'Flux density [Jy]', $
	mxtitsize = titlesize, $
	mytitsize = titlesize, $
	mxtitoffset = 1.2, $
	mytitoffset = 1.2, $
	mcharthick = cthick, $
	mxcharthick = cthick, $
	mycharthick = cthick

for j = 0, 3 do begin

	tag, alltag(j), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+alltag(j)+'.sav'
	
	flux_lh = sed.flux_lh
	wave_lh = sed.wave_lh
	order_lh = sed.order_lh
	
	if alltag[j] ne 'mega023' and alltag[j] ne 'mega027' then begin
		noneg_lh = where(sed.flux_lh gt 0)
		flux_lh = flux_lh(noneg_lh)
		wave_lh = wave_lh(noneg_lh)
		order_lh = order_lh(noneg_lh)
	endif
	
	; Identified IR lines for overplotting
	
	templines = ir_lines(/hr)
	lines = templines(*,0) & line_id = templines(*,1)
	
	; Plot data
	
	ymin = min(flux_lh)
	ymax = max(flux_lh) * 1.1

		altlist = 'IRAS '+['01418+1651','09320+6134','12540+5708','13428+5608','15327+2340']
		alt_title =['III Zw 35','UGC 5101','Mrk 231','Mrk 273','Arp 220'] 
		title = sed.obj
		tt = where(title eq altlist,count)
		if count eq 1 then title = sed.obj+' ('+alt_title[tt]+')'

	plot, wave_lh, flux_lh, $
		xrange = [floor(min(wave_lh)),ceil(max(wave_lh))], /xstyle, $
		yrange = [(ymin - 0.15 * ymax)/0.85,ymax], /ystyle, $
;		xtitle = 'Wavelength (rest frame) [!7l!3m]', $
;		ytitle = 'Flux density [Jy]', $
;		title = title, $
		charsize = 1.1, $
		thick = lthick, $
		charthick = lthick, $
		xthick = lthick, $
		ythick = lthick, $
		/nodata

	oplot, wave_lh,flux_lh, psym = 10, thick = lthick
	
	xyouts, (!x.crange[0]+0.05*(!x.crange[1] - !x.crange[0])), (!y.crange[0]+0.6*(!y.crange[1] - !y.crange[0])), $
		title, charsize = 1.1, charthick = lthick

	; Label detected lines

	thelines = whichlines(alltag[j])
	thelines_label = thelines
	fcount = n_elements(thelines)

	; Adjust labels

	nev24 = where(thelines eq 'neV24', cnev24)
	if cnev24 gt 0 then thelines_label[nev24] = 'neV'
	feII26 = where(thelines eq 'feII26', cfeII26)
	if cfeII26 gt 0 then thelines_label[feII26] = 'feII'
	siii33 = where(thelines eq 'sIII33', csiii33)
	if csiii33 gt 0 then thelines_label[siii33] = 'sIII'
	siii34 = where(thelines eq 'siII34', csiii34)
	if csiii34 gt 0 then thelines_label[siii34] = 'siII'

	if fcount gt 0 then begin

		linewaves = dblarr(fcount)
		for i = 0, fcount - 1 do linewaves[i] = hrwavelength(thelines[i])

		; Remove H2

		noh2 = where(strmid(thelines,0,2) ne 'h2',noh2count)
		h2 = where(strmid(thelines,0,2) eq 'h2',h2count)

		if noh2count gt 0 then begin

			fsid = thelines[noh2]
			fsid_label = thelines_label[noh2]
			fslines = linewaves[noh2]
			for i = 0, n_elements(fsid) -1 do begin
				fsid[i] = strupcase(strmid(fsid[i],0,1)) + strmid(fsid[i],1,strlen(fsid[i])-1)
				fsid_label[i] = strupcase(strmid(fsid_label[i],0,1)) + $
					strmid(fsid_label[i],1,strlen(fsid_label[i])-1)
			endfor

			hualpha = where(fsid eq 'HI76',hucount)
			if hucount gt 0 then fsid[hualpha] = 'Hu-!7a!3'

			chlorine = where(fsid eq 'ClII',clcount)
			if clcount gt 0 then fsid[chlorine] = 'Cl II'

			fsind = where(fslines lt max(wave_lh) and fslines gt min(wave_lh),fscount)
			if fscount gt 0 then begin

				fslines = fslines[fsind] & fsid = fsid[fsind] & fsid_label = fsid_label[fsind]
				bheight = (!y.crange[0] + 0.15*(!y.crange[1] - !y.crange[0]))
				blevel = (!y.crange[0]+0.1*(!y.crange[1] - !y.crange[0]))
				blabel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
				blabellev = (!y.crange[0]+0.12*(!y.crange[1] - !y.crange[0]))
	
				plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
				for k=0,n_elements(fslines)-1 do begin
					plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
					if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.4,blabellev,fsid[k],$
						charsize=labelsize, charthick=lthick $
					else if fsid[k] eq 'OIV' then xyouts,fslines[k]-0.4,blabellev,fsid[k],$
						charsize=labelsize, charthick=lthick $
					else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
						then xyouts,fslines[k]-0.3,blabellev,fsid_label[k],$
							charsize=labelsize, charthick=lthick $
					else xyouts,fslines[k]+0.1,blabellev,fsid_label[k],$
						charsize=labelsize, charthick=lthick
				endfor
				
			endif

		endif

		if h2count gt 0 then begin

			h2id = 'S('+strmid(thelines[h2],3,4)+')'
			h2lines = linewaves[h2]

			h2ind = where(h2lines lt max(wave_lh) and h2lines gt min(wave_lh),h2count2)
			if h2count2 gt 0 then begin

				h2lines = h2lines[h2ind] & h2id = h2id[h2ind]
				bheight = (!y.crange[1]-0.15*(!y.crange[1] - !y.crange[0]))
				blevel = (!y.crange[0]+0.90*(!y.crange[1] - !y.crange[0]))
				blabel = (!y.crange[0]+0.92*(!y.crange[1] - !y.crange[0]))
				blabellev = (!y.crange[0]+0.83*(!y.crange[1] - !y.crange[0]))
			
				plots, [h2lines(0),h2lines(n_elements(h2lines)-1)], [blevel,blevel], thick = lthick
				for k=0,n_elements(h2lines)-1 do begin
					plots, [h2lines(k),h2lines(k)], [blevel, bheight], thick = lthick
					xyouts,h2lines[k]+0.1,blabellev,h2id[k],charsize=labelsize, charthick=lthick
				endfor
				xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',$
					charsize=labelsize, charthick=lthick

			endif

		endif

	endif

	multiplot

endfor	

device,/close
set_plot,'x'

multiplot, /reset

end
