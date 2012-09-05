;pro lhfigures
;+
; NAME:
;       
;	LHFIGURES
;
; PURPOSE:
;
;	Create figures displaying IRS LH spectra for paper
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
;	Added archived and control sample objects - Oct 08
; 	Added alternate (common) titles for certain objects, re HWWS suggestion - Dec 08
;	Modified for 51 total OHMs - Jun 10
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

for m = 0,11 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lhspectra'+string(m+1,format='(i02)')+'.ps', xs=22, ys=24,/portrait, /encap

	!p.font=-1
	
	lthick = 5
	titlesize = 1.5
	labelsize = 0.9
	cthick = 5
	
	multiplot, [1,4], $
		ytickformat='(f4.1)', $
		ygap=0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick


	for j = 4*m,4*m + 3 do begin
	
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
			charsize = titlesize, $
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
					bheight = (!y.crange[0] + 0.10*(!y.crange[1] - !y.crange[0]))
					blevel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
					blabellev = (!y.crange[0]+0.07*(!y.crange[1] - !y.crange[0]))
		
					plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
					for k=0,n_elements(fslines)-1 do begin
						plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
						if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.6,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'OIV' then xyouts,fslines[k]-0.5,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
							then xyouts,fslines[k]-0.5,blabellev,fsid_label[k],$
								charsize=labelsize, charthick=lthick $
						else xyouts,fslines[k]+0.06,blabellev,fsid_label[k],$
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
						xyouts,h2lines[k]+0.06,blabellev,h2id[k],charsize=labelsize, charthick=lthick
					endfor
					xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=labelsize, charthick=lthick
	
				endif
	
			endif

		endif

		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

endfor

; Last plot 1x3

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lhspectra'+string(13,format='(i02)')+'.ps', xs=22, ys=24,/portrait, /encap

	!p.font=-1
	
	lthick = 5
	titlesize = 1.5
	labelsize = 0.9
	cthick = 5
	
	multiplot, [1,3], $
		ytickformat='(f4.1)', $
		ygap=0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick


	for j = 48,50 do begin
	
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
			charsize = titlesize, $
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
					bheight = (!y.crange[0] + 0.10*(!y.crange[1] - !y.crange[0]))
					blevel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
					blabellev = (!y.crange[0]+0.07*(!y.crange[1] - !y.crange[0]))
		
					plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
					for k=0,n_elements(fslines)-1 do begin
						plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
						if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.6,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'OIV' then xyouts,fslines[k]-0.5,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
							then xyouts,fslines[k]-0.5,blabellev,fsid_label[k],$
								charsize=labelsize, charthick=lthick $
						else xyouts,fslines[k]+0.06,blabellev,fsid_label[k],$
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
						xyouts,h2lines[k]+0.06,blabellev,h2id[k],charsize=labelsize, charthick=lthick
					endfor
					xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=labelsize, charthick=lthick
	
				endif
	
			endif

		endif

		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

; Retreive all control LH spectra, sorted by IRAS number

contag = condat('tag')
conobj = condat('obj')

consort = sort(conobj)

contag = contag[consort]
conobj = conobj[consort]

for m = 0,2 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lhspectra_con'+string(m+1,format='(i02)')+'.ps', xs=22, ys=24,/portrait, /encap

	multiplot, [1,4], $
		ytickformat='(f4.1)', $
		ygap=0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick


	for j = 4*m,4*m + 3 do begin
	
		tag, contag(j), dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+contag(j)+'.sav'
		
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
		
		ymin = min(flux_lh)
		ymax = max(flux_lh) * 1.1

			title = sed.obj
			if title eq 'IRAS 01572+0009' then title = sed.obj+' (Mrk 1014)'

		plot, wave_lh, flux_lh, $
			xrange = [floor(min(wave_lh)),ceil(max(wave_lh))], /xstyle, $
			yrange = [(ymin - 0.15 * ymax)/0.85,ymax], /ystyle, $
			charsize = titlesize, $
			thick = lthick, $
			charthick = lthick, $
			xthick = lthick, $
			ythick = lthick, $
			/nodata
	
		oplot, wave_lh,flux_lh, psym = 10, thick = lthick
		
	xyouts, (!x.crange[0]+0.05*(!x.crange[1] - !x.crange[0])), (!y.crange[0]+0.6*(!y.crange[1] - !y.crange[0])), $
		title, charsize = 1.1, charthick = lthick

		; Label detected lines

		thelines = whichlines(contag[j])
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

			; Plot lines that are not H2

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
					bheight = (!y.crange[0] + 0.10*(!y.crange[1] - !y.crange[0]))
					blevel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
					blabellev = (!y.crange[0]+0.07*(!y.crange[1] - !y.crange[0]))
		
					plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
					for k=0,n_elements(fslines)-1 do begin
						plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
						if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.6,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'OIV' then xyouts,fslines[k]-0.5,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
							then xyouts,fslines[k]-0.5,blabellev,fsid_label[k],$
								charsize=labelsize, charthick=lthick $
						else xyouts,fslines[k]+0.06,blabellev,fsid_label[k],$
							charsize=labelsize, charthick=lthick
					endfor
					
				endif

			endif
	
			; Plot H2 lines

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
						xyouts,h2lines[k]+0.06,blabellev,h2id[k],charsize=labelsize, charthick=lthick
					endfor
					xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=labelsize, charthick=lthick
	
				endif
	
			endif

		endif

		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

endfor
	
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lhspectra_con'+string(4,format='(i02)')+'.ps', xs=22, ys=24,/portrait, /encap

	multiplot, [1,3], $
		ytickformat='(f4.1)', $
		ygap=0.005, $
		mxtitle = 'Rest wavelength [!7l!3m]', $
		mytitle = 'S!I!7k!3!N [Jy]', $
		mxtitsize = titlesize, $
		mytitsize = titlesize, $
		mxtitoffset = 1.2, $
		mytitoffset = 1.2, $
		mcharthick = cthick, $
		mxcharthick = cthick, $
		mycharthick = cthick


	for j = 12,14 do begin
	
		tag, contag(j), dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+contag(j)+'.sav'
		
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
		
		ymin = min(flux_lh)
		ymax = max(flux_lh) * 1.1

			title = sed.obj
			if title eq 'IRAS 01572+0009' then title = sed.obj+' (Mrk 1014)'

		plot, wave_lh, flux_lh, $
			xrange = [floor(min(wave_lh)),ceil(max(wave_lh))], /xstyle, $
			yrange = [(ymin - 0.15 * ymax)/0.85,ymax], /ystyle, $
			charsize = titlesize, $
			thick = lthick, $
			charthick = lthick, $
			xthick = lthick, $
			ythick = lthick, $
			/nodata
	
		oplot, wave_lh,flux_lh, psym = 10, thick = lthick
		
	xyouts, (!x.crange[0]+0.05*(!x.crange[1] - !x.crange[0])), (!y.crange[0]+0.6*(!y.crange[1] - !y.crange[0])), $
		title, charsize = 1.1, charthick = lthick

		; Label detected lines

		thelines = whichlines(contag[j])
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

			; Plot lines that are not H2

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
					bheight = (!y.crange[0] + 0.10*(!y.crange[1] - !y.crange[0]))
					blevel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
					blabellev = (!y.crange[0]+0.07*(!y.crange[1] - !y.crange[0]))
		
					plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
					for k=0,n_elements(fslines)-1 do begin
						plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
						if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.6,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'OIV' then xyouts,fslines[k]-0.5,blabellev,fsid[k],$
							charsize=labelsize, charthick=lthick $
						else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
							then xyouts,fslines[k]-0.5,blabellev,fsid_label[k],$
								charsize=labelsize, charthick=lthick $
						else xyouts,fslines[k]+0.06,blabellev,fsid_label[k],$
							charsize=labelsize, charthick=lthick
					endfor
					
				endif

			endif
	
			; Plot H2 lines

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
						xyouts,h2lines[k]+0.06,blabellev,h2id[k],charsize=labelsize, charthick=lthick
					endfor
					xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=labelsize, charthick=lthick
	
				endif
	
			endif

		endif

		multiplot

	endfor	
	
	device,/close
	set_plot,'x'

end
