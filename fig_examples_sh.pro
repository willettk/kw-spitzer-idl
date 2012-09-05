; SH

ohmtag = ohmdat('tag')
archtag = archdat('tag')

ohmobj = ohmdat('obj')
archobj = archdat('obj')

allobj = [transpose(ohmobj),transpose(archobj)]
alltag = [transpose(ohmtag),transpose(archtag)]

objsort = sort(allobj)

allobj = allobj[objsort]
alltag = alltag[objsort]

set_plot,'ps'
device, filename='~/Astronomy/Research/Spitzer/papers/apj_submission/paperI/fig3.ps', xs=22, ys=24

!p.font = -1

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

		altlist = 'IRAS '+['01418+1651','09320+6134','12540+5708','13428+5608','15327+2340']
		alt_title =['III Zw 35','UGC 5101','Mrk 231','Mrk 273','Arp 220'] 
		title = sed.obj
		tt = where(title eq altlist,count)
		if count eq 1 then title = sed.obj+' ('+alt_title[tt]+')'

	plot, wave_sh, flux_sh, $
		xrange = [floor(min(wave_sh)),ceil(max(wave_sh))+1], /xstyle, $
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

	oplot, wave_sh,flux_sh, psym = 10, thick = lthick
	
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
			blabel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
			blabellev = (!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))

			plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
			for k=0,n_elements(fslines)-1 do begin
				plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
				if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.4,$
					blabellev,fsid[k],charsize=labelsize, charthick=lthick $
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
				xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',$
					charsize=labelsize, charthick=lthick

			endif

		endif

	endif

	multiplot

endfor	

device,/close
set_plot,'x'

multiplot,/reset

end
