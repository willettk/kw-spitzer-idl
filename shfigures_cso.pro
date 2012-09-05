;pro shfigures_cso, stop = stop
;+
; NAME:
;       
;	SHFIGURES_CSO
;
; PURPOSE:
;
;	Create figures displaying IRS SH spectra for CSO paper
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
;	IDL> shfigures_cso
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;	Removed NGC 5793, 1245+676		Nov 09
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
	device, filename='~/Astronomy/Research/Spitzer/cso/papers/shspectra_cso.ps', xs=18, ys=24,/portrait
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


for j = 0, ncso-1 do begin

;	set_plot,'ps'
;	device, filename='~/Astronomy/Research/Spitzer/cso/papers/shspectra_cso'+strtrim(m+1,2)+'.ps', xs=18, ys=24,/portrait
;;	device, filename='~/Astronomy/Research/Spitzer/cso/papers/shspectra_cso.ps', xs=18, ys=24,/portrait
;	for j = 4*m,4*m + 3 do begin
;;	for j = 0,ncso-1 do begin
	
		tag, csotag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+csotag[j]+'.sav'
		
		flux_sh = sed.flux_sh
		wave_sh = sed.wave_sh ;* (1d + sed.redshift)
		order_sh = sed.order_sh
		
		noneg_sh = where(sed.flux_sh gt 0)
;		flux_sh = flux_sh(noneg_sh)
;		wave_sh = wave_sh(noneg_sh)
;		order_sh = order_sh(noneg_sh)
		
		; Plot data
		
		; Hard copy option
		
		lthick = 3
		cs = 1
		specthick = 2
		
		ymin = min(flux_sh)
		ymax = max(flux_sh) * 1.1

		plot, wave_sh, flux_sh, $
;			xrange = [floor(min(wave_sh)),ceil(max(wave_sh))], /xstyle, $
			xrange = [7.5,20], /xstyle, $
;			xrange = [9.5,20], /xstyle, $
			yrange = [(ymin - 0.15 * ymax)/0.85,ymax], /ystyle, $
;			xtitle = 'Wavelength (rest frame) [!7l!3m]', $
;			ytitle = 'Flux density [Jy]', $
;			title = sed.obj, $
			yminor = 3, $
			ytickformat='(f4.2)', $
			charsize = 1.0, $
			xthick = 4, $
			ythick = 4, $
			charthick = 4, $
			/nodata
	
		oplot, wave_sh,flux_sh, psym = 10, thick = specthick
		
		pdotposition = !p.position
		xyouts, 0.15, pdotposition[3] - 0.025, sed.obj, /normal, charthick = 4

		; Label detected lines

		thelines = whichlines(csotag[j])
		fcount = n_elements(thelines)

		if fcount gt 0 and strtrim(thelines[0],2) ne '0' then begin

			linewaves = dblarr(fcount)
			for i = 0, fcount - 1 do linewaves[i] = hrwavelength(thelines[i]) ;* (1d + sed.redshift)

			; Remove H2 from line array

			noh2 = where(strmid(thelines,0,2) ne 'h2',noh2count)
			h2 = where(strmid(thelines,0,2) eq 'h2',h2count)

			if noh2count gt 0 then begin

				; Locate fine-structure lines

				fsid = thelines[noh2]
				fslines = linewaves[noh2]
				for i = 0, n_elements(fsid) -1 do $
					fsid[i] = strupcase(strmid(fsid[i],0,1)) + strmid(fsid[i],1,strlen(fsid[i])-1)

				; Manual labels for HI 7-6, CLII

				hualpha = where(fsid eq 'HI76',hucount)
				if hucount gt 0 then fsid[hualpha] = 'Hu-!7a!3'

				chlorine = where(fsid eq 'ClII',clcount)
				if clcount gt 0 then fsid[chlorine] = 'Cl II'

				; Find lines within SH module

				fsind = where(fslines lt max(wave_sh) and fslines gt min(wave_sh))
				fslines = fslines[fsind] & fsid = fsid[fsind]

				; Define locations on plot in normalized CS

				bheight = (!y.crange[0] + 0.15*(!y.crange[1] - !y.crange[0]))
				blevel = (!y.crange[0]+0.1*(!y.crange[1] - !y.crange[0]))
				blabel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
				blabellev = (!y.crange[0]+0.12*(!y.crange[1] - !y.crange[0]))
				pahlabellev = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))

				; Add PAH lines for visible features in HR (request by Perlman)
				;	Lines from LRFIGURES_CSO, all lines longward of 9 um

				;case j of
				;	0: pahlines =             [        11.3,12.7,     16.4,17.1]
				;	1: pahlines =         [            11.3,12.7,14.2,16.4]
				;	2: pahlines = [                    11.3,12.7,14.2,16.4,     17.4]
				;	3: pahlines =         [            11.3,12.7,14.2,16.4]
				;	4: pahlines = [0]
				;	5: pahlines =         [            11.3,12.7,14.2]
				;	6: pahlines = [                    11.3,12.7,14.2               ]
				;	7: pahlines = [                    11.3,12.7                    ]
				;	8: pahlines = [                    11.3,12.7,     16.4,     17.4]
				;	9: pahlines = [0]                                                
				;endcase

				;case j of
				;	0: pahlines =         [            11.3,12.7,14.2,16.4]
				;	1: pahlines = [                    11.3,12.7                    ]
				;	2: pahlines =             [        11.3,12.7,     16.4,17.1]
				;	3: pahlines = [0]                                                
				;	4: pahlines =         [            11.3,12.7,14.2]
				;	5: pahlines =         [            11.3,12.7,14.2,16.4]
				;	6: pahlines = [0]
				;	7: pahlines = [                    11.3,12.7,14.2,16.4,     17.4]
				;	8: pahlines = [                    11.3,12.7,     16.4,     17.4]
				;	9: pahlines = [                    11.3,12.7,14.2               ]
				;endcase

				case j of
					0: pahlines =         [            11.3,12.7,14.2,16.4]
					1: pahlines = [                    11.3,12.7                    ]
					2: pahlines =             [        11.3,12.7,     16.4,17.1]
					3: pahlines =         [            11.3,12.7,14.2]
					4: pahlines =         [            11.3,12.7,14.2,16.4]
					5: pahlines = [0]
					6: pahlines = [                    11.3,12.7,     16.4,     17.4]
					7: pahlines = [                    11.3,12.7,14.2               ]
				endcase

				if pahlines[0] ne 0 then begin
					fslines = [fslines,pahlines]
					fsid = [fsid,replicate(' ',n_elements(pahlines))]
				endif

				; Plot labels for fine structure lines

				plots, [min(fslines),max(fslines)], [blevel,blevel], thick = lthick
				for k=0,n_elements(fslines)-1 do begin
					plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
					if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.4,blabellev,fsid[k],charsize=0.5, charthick=lthick $
						else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
						then xyouts,fslines[k]-0.3,blabellev,fsid[k],charsize=0.5, charthick=lthick $
						else xyouts,fslines[k]+0.1,blabellev,fsid[k],charsize=0.5, charthick=lthick
				endfor

				pahind = where(fsid eq ' ', npah)
				if npah gt 0 then begin
					for k = 0, npah - 1 do begin
						xyouts, fslines[pahind[k]] - 0.1, pahlabellev, 'PAH', charsize = 0.4, charthick = lthick
					endfor
				endif

			endif

			; Plot H2 line labels
	
			if h2count gt 0 then begin
	
				h2id = 'S('+strmid(thelines[h2],3,4)+')'
				h2lines = linewaves[h2]
	
				h2ind = where(h2lines lt max(wave_sh))
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

	multiplot

endfor	

multiplot, /reset

device,/close
set_plot,'x'

if keyword_set(stop) then stop

end
