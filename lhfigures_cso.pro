;pro lhfigures_cso, stop = stop
;+
; NAME:
;       
;	LHFIGURES_CSO
;
; PURPOSE:
;
;	Create figures displaying IRS LH spectra for CSO paper
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
;	Removed NGC 5793, 1245+676 - Nov 09
; 	Switched to single page, MULTIPLOT - Nov 09
;-

csotag = csodat('tag')

;bind = [1,7,0,9,5,3,4,2,8,6]
bind = [1,7,0,5,3,4,8,6]
csotag = csotag[bind]
ncso = n_elements(csotag)

!p.font = -1
;!p.multi=[0,1,4]

;for m = 0,1 do begin

	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/cso/papers/lhspectra_cso.ps', xs=18, ys=24,/portrait
;	!p.multi=[0,1,ncso]

erase 

multiplot, [1,ncso], $
	mxtitle = 'Wavelength (rest frame) [!7l!3m]', $
	mytitle = 'Flux density [Jy]', $
	mxtitsize = 1.2, $
	mytitsize = 1.5, $
	mcharthick = 4, $
	mxcharthick = 4, $
	mycharthick = 4

for j = 0,ncso-1 do begin

;	set_plot,'ps'
;	device, filename='~/Astronomy/Research/Spitzer/cso/papers/lhspectra_cso'+strtrim(m+1,2)+'.ps', xs=18, ys=24,/portrait
;	for j = 4*m,4*m + 3 do begin
;	for j = 0,5 do begin

		tag, csotag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+csotag[j]+'.sav'
		
		flux_lh = sed.flux_lh
		wave_lh = sed.wave_lh ;* (1d + sed.redshift)
		order_lh = sed.order_lh
		
		noneg_lh = where(sed.flux_lh gt 0)

		; Kluge for the low S/N object VII Zw 485

		if csotag[j] ne 'cso010' then begin
			flux_lh = flux_lh(noneg_lh)
			wave_lh = wave_lh(noneg_lh)
			order_lh = order_lh(noneg_lh)
		endif
		
		; Identified IR lines for overplotting
		
		templines = ir_lines(/hr)
		lines = templines(*,0) & line_id = templines(*,1)
		
		; Plot data
		
		; Hard copy option
		
		lthick = 3
		specthick = 2
		cs = 1
		
		ymin = min(flux_lh)
		ymax = max(flux_lh) * 1.1

		plot, wave_lh, flux_lh, $
;			xrange = [floor(min(wave_lh)),ceil(max(wave_lh))], /xstyle, $
			yrange = [(ymin - 0.15 * ymax)/0.85,ymax], /ystyle, $
			xrange = [15,36], /xstyle, $
			;xtitle = 'Wavelength (rest frame) [!7l!3m]', $
			;ytitle = 'Flux density [Jy]', $
			;title = sed.obj, $
;			ytickv = indgen(nyticks) , $
			yminor = 3, $
			ytickformat='(f4.2)', $
			charsize = 1.0, $
			xthick = 4, $
			ythick = 4, $
			charthick = 4, $
			/nodata
	
		oplot, wave_lh,flux_lh, psym = 10, thick = specthick

		pdotposition = !p.position
		xyouts, 0.15, pdotposition[3] - 0.025, sed.obj, /normal, charthick = 4
		
		; Label detected lines

		thelines = whichlines(csotag[j])
		fcount = n_elements(thelines)

		if fcount gt 0 and strtrim(thelines[0],2) ne '0' then begin

			linewaves = dblarr(fcount)
			for i = 0, fcount - 1 do linewaves[i] = hrwavelength(thelines[i]) ;* (1d + sed.redshift)

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

				numion = where(fsid eq 'NeV24',nicount)
				if nicount gt 0 then fsid[numion] = 'NeV'

				numion = where(fsid eq 'SIII33',nicount)
				if nicount gt 0 then fsid[numion] = 'SIII'

				numion = where(fsid eq 'SiII34',nicount)
				if nicount gt 0 then fsid[numion] = 'SiII'

				numion = where(fsid eq 'FeII26',nicount)
				if nicount gt 0 then fsid[numion] = 'FeII'

				fsind = where(fslines lt max(wave_lh) and fslines gt min(wave_lh),fsindcount)
				
				if fsindcount gt 0 then begin
					fslines = fslines[fsind] & fsid = fsid[fsind]
					bheight = (!y.crange[0] + 0.15*(!y.crange[1] - !y.crange[0]))
					blevel = (!y.crange[0]+0.1*(!y.crange[1] - !y.crange[0]))
					blabel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
					blabellev = (!y.crange[0]+0.12*(!y.crange[1] - !y.crange[0]))
		
					plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
					for k=0,n_elements(fslines)-1 do begin
						plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
						if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.4,blabellev,fsid[k],charsize=0.5, charthick=lthick $
							else if fsid[k] eq 'OIV' $
							then xyouts,fslines[k]-0.4,blabellev,fsid[k],charsize=0.5, charthick=lthick $
							else xyouts,fslines[k]+0.1,blabellev,fsid[k],charsize=0.5, charthick=lthick
					endfor
				endif

			endif
	
			if h2count gt 0 then begin
	
				h2id = 'S('+strmid(thelines[h2],3,4)+')'
				h2lines = linewaves[h2]
	
				h2ind = where(h2lines lt max(wave_lh) and h2lines gt min(wave_lh),h2count)
				if h2count gt 0 then begin
					h2lines = h2lines[h2ind] & h2id = h2id[h2ind]
					bheight = (!y.crange[1]-0.15*(!y.crange[1] - !y.crange[0]))
					blevel = (!y.crange[0]+0.90*(!y.crange[1] - !y.crange[0]))
					blabel = (!y.crange[0]+0.92*(!y.crange[1] - !y.crange[0]))
					blabellev = (!y.crange[0]+0.83*(!y.crange[1] - !y.crange[0]))
				
					plots, [h2lines(0),h2lines(n_elements(h2lines)-1)], [blevel,blevel], thick = lthick
					for k=0,n_elements(h2lines)-1 do begin
						plots, [h2lines(k),h2lines(k)], [blevel, bheight], thick = lthick
						xyouts,h2lines[k]+0.1,blabellev,h2id[k],charsize=0.5, charthick=lthick
					endfor
					xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=0.6, charthick=lthick
				endif
	
			endif

		endif
	
	multiplot

endfor	

multiplot, /reset

device,/close
set_plot,'x'

;endfor
	
if keyword_set(stop) then stop

end
