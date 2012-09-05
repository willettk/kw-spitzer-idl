pro h2temp_cso, verbose = verbose, ps = ps
;+
; NAME:
;       
;	H2TEMP
;
; PURPOSE:
;
;	Compute H2 excitation temperature of gas using H2 lines for CSOs
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
;	VERBOSE - 	prints data to screen
;
;	PS - 		make hard copies of figures
;
; EXAMPLE:
;
;	IDL> h2temp_cso, /ps
;
; NOTES:
;
;	Should be changed so that it shows limits for all objects
;
; REVISION HISTORY
;	Adapted from H2TEMP.pro	- Aug 08
;	Added label on y-axis 		Sep 09
;	Removed NGC 5793, 1245+676 - Nov 09
;-

if not keyword_set(verbose) then quiet = 1

; Data lists

taglist = csodat('tag',/ra) & taglist = taglist[1,*]
objlist = csodat('obj',/ra) & objlist = objlist[1,*]
dllist = csodat('dl',/ra)   & dllist  = dllist[1,*]

; Remove NGC 5793, 1245+676

csoind = [0,1,2,4,5,6,8,9]

taglist = taglist[csoind]
objlist = objlist[csoind]
dllist  = dllist[csoind]

ncso = n_elements(taglist)

h2 = {obj:strarr(ncso), tag:strarr(ncso), $
	temp_warm:fltarr(ncso), $
	temp_hot:fltarr(ncso), $
	mass_warm:fltarr(ncso), $
	mass_hot:fltarr(ncso), $
	temp_warmerr:fltarr(ncso), $
	temp_hoterr:fltarr(ncso), $
	;mass_warmerr:fltarr(ncso), $
	;mass_hoterr:fltarr(ncso), $
	warm_avg:0d, hot_avg:0d}
 
; Physical constants for H2 transitions taken from 						
; http://www.not.iac.es/instruments/notcam/ReferenceInfo/h2_lines.html

h2_indices = [0,1,2,3,4,5,6,7]						; All H2 transitions seen in sample(s)
ej = [510, 1015, 1682, 2504, 3474, 4586, 5829, 7197]			; Energy of upper state
aj = [3d-4, 4.8d-3, 2.76d-2, 9.84d-2, 0.264, 0.588, 1.14, 2.00]		; Einstein-A coefficient
lambda = [28.221,17.035,12.279,9.6649,8.0258,6.9091,6.1088,5.5115]	; Rest wavelength of transition
gj = [5,21,9,33,13,45,17,57]						; Statistical weight of upper level

plancks = 6.6261d-27					; Planck's constant [cgs]
c = 3d10						; Speed of light [cm/s]
del_ej = plancks * c / lambda				; Energy of transition

; Linear fit for systems with a single excitation temperature

expr = 'p[0] + x*p[1]'
start = [10, 3d-3]

filenames = '~/Astronomy/Research/Spitzer/cso/papers/h2temp_cso.ps'

!p.multi=[0,4,2]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename=filenames, /portrait, xs=18, ys=12;, xoff=1, yoff=1
	tempcs = 1.0
	axescs = 1.0
	labelcs = 0.5
	defthick = 4
	tempthick = 4
	cthick = 4
	psymsize = 0.5
	arrowsize = 100
endif else begin
	set_plot,'x'
	tempcs = 1
	axescs = 1
	labelcs = 1
	defthick = 1
	tempthick = 1
	cthick = 1
	psymsize = 1
	arrowsize = 13
endelse

!p.font=-1

multiplot, [4,2], $
	mxtitle = 'E!IJ!N/k!Ib!N [K]', $
	mytitle = 'ln(N!IJ!N / g!IJ!N)', $
	mxtitsize = 1.2, $
	mytitsize = 1.2, $
	mcharthick = 4, $
	mxcharthick = 4, $
	mycharthick = 4

for j = 0, ncso-1 do begin

	fname = taglist[j]
	obj = objlist(where(taglist eq fname))

	; Empty plot window
	
	xleft = 0.05
	xright = 0.95
	ybottom = 0.05
	ytop = 0.95

	y0 = [replicate(ybottom + (ytop-ybottom)/2.,4),replicate(ybottom,4)]
	y1 = [replicate(ytop,4),replicate(ybottom + (ytop-ybottom)/2.,4)]

	nxticks = 3
	nyticks = 5
	if j/4 eq 0 then begin
		xtitle = ' '
		xticknames = replicate(' ',nxticks)
		titlelevel = -2
	endif else begin
		xtitle = 'E!IJ!N/k!Ib!N [K]'
		xticknames = [strtrim(string(lindgen(nxticks - 1) * 10000/(nxticks-1),format='(i4)'),2), ' ']
		titlelevel = -2
	endelse

	if j mod 4 eq 0 then ytitle = 'ln(N!IJ!N / g!IJ!N)' else ytitle = ''
	if j mod 4 ne 0 then yticknames = replicate(' ',nyticks) else $
		yticknames = reverse([strtrim(string((-1) * indgen(nyticks - 1) * 20/(nyticks-1),format='(i4)'),2), ' '])

	plot, indgen(10), $
		/nodata, $
;		xtitle = xtitle, $
;		ytitle = ytitle, $
;		ytitle = 'ln(N!IJ!N / g!IJ!N)', $
;		title = obj, $
;		position = [xleft + (j mod 4) * (xright - xleft) / 4., y0[j], $
;			(3 * xleft + xright)/4. + (j mod 4) * (xright - xleft) / 4., y1[j]], $
		charsize = axescs, $
		charthick = cthick, $
		thick = defthick, $
		xthick = defthick, $
		ythick = defthick, $
		xtickv = indgen(nxticks) * (10000 / (nxticks - 1)), $
		xticks = nxticks - 1, $
		xtickname = xticknames, $
		ytickv = (-1) * reverse(indgen(nyticks) * 20/(nyticks - 1)), $
		yticks = nyticks - 1, $
		ytickname = yticknames, $
		xr = [0,9000], $
		yr = [-20,0]

	xyouts, /data, 2500, titlelevel, obj, charsize = 0.8, charthick=cthick

	h2.tag[j] = fname
	h2.obj[j] = obj

	temp_hot = 0.

	; Load in the line fluxes and errors
	
	orthoh2 = ['h2s0','h2s1','h2s2','h2s3','h2s4_lr','arIIh2s5_lr','h2s6_lr','h2s7_lr']

	h2fluxarr = dblarr(n_elements(orthoh2))
	h2errarr  = dblarr(n_elements(orthoh2))
	h2limarr = dblarr(n_elements(orthoh2))
	
	for i = 0, n_elements(h2fluxarr)-1 do begin
		
		h2_query = getlineflux(fname, orthoh2[i],/quiet)
		if h2_query[0] ne -1 then begin
;			if keyword_set(verbose) then print,strupcase(orthoh2[i])+' line found for '+fname
			h2fluxarr[i] = h2_query[0]
			h2errarr[i] = h2_query[1]
		endif else begin
			if i gt 3 then limit = linelim(fname,orthoh2[i],/lores,/noplot,/quiet) * 1d-21 $
				else begin
					limit = linelim(fname,orthoh2[i],/noplot,/quiet) * 1d-21
					if limit eq -1d-21 then limit = linelim(fname,orthoh2[i],/lh,/noplot,/quiet) * 1d-21
				endelse
			h2limarr[i] = limit
		endelse
	
	endfor
	
	; Plot the excitation diagrams					
	
	hlabels = ['S(0)', 'S(1)', 'S(2)', 'S(3)', 'S(4)', 'S(5)', 'S(6)', 'S(7)']
	
	; Plot all hot points that exist w/error bars
	; For S(5) and other limits, plot them with an arrow
	; Check if hot temperature exists - if so, compute
	; Subtract hot component from warm lines and plot the warm points w/error bars
	; Plot warm limits with arrows
	; Compute warm temperature

		; Plot points that exist from S(4) through S(7)

		hhot = where(h2fluxarr[4:7] ne 0) + 4

		if hhot[0] gt 3 then begin

			nj_hot		 = h2fluxarr[hhot] / (aj[hhot] * del_ej[hhot])
			njerr_hot	 = h2errarr[hhot]  / (aj[hhot] * del_ej[hhot])
	
			oplot, ej[hhot], alog(nj_hot / gj[hhot]), psym = symcat(16), symsize = psymsize
			xyouts, ej[hhot]+250, alog(nj_hot / gj[hhot]), hlabels[hhot], charsize = labelcs, charthick = cthick

			oploterror, ej[hhot], alog(nj_hot / gj[hhot]), alog(nj_hot / gj[hhot]) $
				- alog((nj_hot - njerr_hot) / gj[hhot]), psym = 3, /lobar
			oploterror, ej[hhot], alog(nj_hot / gj[hhot]), alog(nj_hot / gj[hhot]) $
				- alog((nj_hot + njerr_hot) / gj[hhot]), psym = 3, /hibar
		
			if h2fluxarr[5] ne 0 then begin
				nj5    = h2fluxarr[5] / (aj[5] * del_ej[5])
				njerr5 = h2errarr[5]  / (aj[5] * del_ej[5])

				if nj5 lt njerr5 then begin
					print,'S(5) error is larger than flux for '+fname
					njerr5 = 0.5 * nj5
				endif
				oploterror, ej[5], alog(nj5 / gj[5]), $
					alog(nj5 / gj[5]) - alog((nj5 - njerr5) / gj[5]), $
					psym = 3, /lobar
				oploterror, ej[5], alog(nj5 / gj[5]), $
					alog(nj5 / gj[5]) - alog((nj5 + njerr5) / gj[5]), $
					psym = 3, /hibar
			endif
	
		endif

		; Upper limits on H2 S(5)

		; S(5)

		if h2fluxarr[5] ne 0 then $
			arrow, ej[5], alog((h2fluxarr[5]/ (aj[5] * del_ej[5])) / gj[5]), $
				ej[5], alog((h2fluxarr[5]/ (aj[5] * del_ej[5])) / gj[5]) - 2, $
				/data, thick=defthick, hsize = arrowsize else begin

				nj5    = h2limarr[5] / (aj[5] * del_ej[5])
				oplot, [ej[5]], [alog(nj5 / gj[5])], psym = symcat(16), symsize = psymsize
				xyouts, ej[5]+250, alog(nj5 / gj[5]), hlabels[5], charsize = labelcs, charthick = cthick
				arrow, ej[5], alog((h2limarr[5]/ (aj[5] * del_ej[5])) / gj[5]), $
					ej[5], alog((h2limarr[5]/ (aj[5] * del_ej[5])) / gj[5]) - 2, $
					/data, thick=defthick, hsize = arrowsize
			endelse

		; If H2 S(7) exists, fit a temperature between S(3) and S(7) hot data points

		if h2fluxarr[3] ne 0 and h2fluxarr[7] ne 0 then begin

			hfi_temp = where(h2fluxarr[3:7] ne 0) + 3
			hotfitind = hfi_temp[where(hfi_temp ne 5)]

			nj_hotfit	 = h2fluxarr[hotfitind] / (aj[hotfitind] * del_ej[hotfitind])
			njerr_hotfit	 = h2errarr[hotfitind]  / (aj[hotfitind] * del_ej[hotfitind])

			hot_result = mpfitexpr(expr, ej[hotfitind], alog(nj_hotfit / gj[hotfitind]), $
				alog(njerr_hotfit / gj[hotfitind]), start, /quiet, perr = hot_err)
			hotarr = fillarr(1,ej[hotfitind[0]],ej[hotfitind[n_elements(hotfitind)-1]])
			oplot, hotarr, (hot_result[0] + hotarr*hot_result[1]), linestyle = 2, thick = tempthick
			temp_hot = -1d / hot_result[1]
			temp_hoterr = 1d / hot_err[1]
;			xyouts, 1500, -18, 'T!Ihot!N='+string(temp_hot, format = '(i3)')+' K', charsize = tempcs, charthick = cthick
			h2.temp_hot[j] = temp_hot
			h2.temp_hoterr[j] = temp_hoterr

		endif

		; Warm lines from S(0) - S(3)

		hwarm = where(h2fluxarr[0:3] ne 0)

		if hwarm[0] ne -1 then begin
			
			; If a hot gas component was detected, subtract expected flux contribution from the warm lines
	
			warmh2 = where(h2fluxarr[0:3] ne 0)
			if temp_hot ne 0 and warmh2[0] ne -1 then begin
				;print,''
				;print,h2errarr
				;print,h2fluxarr
				;print,fname
				;print,''
				orig_flux = h2fluxarr[warmh2]
				ytemp = (hot_result[0] + hot_result[1] * ej[warmh2])
				ysub = aj[warmh2] * del_ej[warmh2] * gj[warmh2] * exp(ytemp)
				h2fluxarr[warmh2] = (h2fluxarr[warmh2] - ysub) + orig_flux[n_elements(orig_flux)-1]


				; Scale errors to new warm flux

				h2errarr[warmh2] = h2errarr[warmh2] * h2fluxarr[warmh2] / orig_flux	

				;print,''
				;print,h2errarr
				;print,h2fluxarr
				;print,fname
				;print,''
			endif
	
			nj_warm     = h2fluxarr[hwarm] / (aj[hwarm] * del_ej[hwarm])
			njerr_warm = h2errarr[hwarm]  / (aj[hwarm] * del_ej[hwarm])
	
			oplot, ej[hwarm], alog(nj_warm / gj[hwarm]), psym = symcat(16), symsize = psymsize
			xyouts, ej[hwarm]+250, alog(nj_warm / gj[hwarm]), hlabels[hwarm], charsize = labelcs, charthick = cthick
	
			oploterror, ej[hwarm], alog(nj_warm / gj[hwarm]), alog(nj_warm / gj[hwarm]) $
				- alog((nj_warm - njerr_warm) / gj[hwarm]), psym = 3, /lobar
			oploterror, ej[hwarm], alog(nj_warm / gj[hwarm]), alog(nj_warm / gj[hwarm]) $
				- alog((nj_warm + njerr_warm) / gj[hwarm]), psym = 3, /hibar
			
			if n_elements(hwarm) gt 1 then begin
				warm_result = mpfitexpr(expr, ej[hwarm], alog(nj_warm / gj[hwarm]), $
					alog(njerr_warm / gj[hwarm]), start, /quiet, perr=warmerr)
				warmarr = fillarr(1,ej[hwarm[0]],ej[hwarm[n_elements(hwarm)-1]])
				oplot, warmarr, (warm_result[0] + warmarr*warm_result[1]), linestyle = 1, thick = tempthick
				temp_warm = -1d / warm_result[1]
				temp_warmerr = 1d / warmerr[1]
				;xyouts, 3000, -5, 'T!Iwarm!N='+string(temp_warm, format = '(i3)')+' K', charsize = tempcs, charthick = cthick
				xyouts, 4000, -7, string(temp_warm, format = '(i3)')+' K', charsize = tempcs, charthick = cthick
	
				h2.temp_warm[j] = temp_warm
				h2.temp_warmerr[j] = temp_warmerr
			endif
	
		endif

		; Plot limits on warm lines

		nohwarm = where(h2limarr[0:3] ne 0)

		if nohwarm[0] ne -1 then begin
			nj_nowarm    = h2limarr[nohwarm] / (aj[nohwarm] * del_ej[nohwarm])
			oplot, ej[nohwarm], alog(nj_nowarm / gj[nohwarm]), psym = symcat(16), symsize = psymsize
			xyouts, ej[nohwarm]+250, alog(nj_nowarm / gj[nohwarm]), hlabels[nohwarm], $
				charsize = labelcs, charthick = cthick
			arrow, ej[nohwarm], alog((h2limarr[nohwarm]/ (aj[nohwarm] * del_ej[nohwarm])) / gj[nohwarm]), $
				ej[nohwarm], alog((h2limarr[nohwarm]/ (aj[nohwarm] * del_ej[nohwarm])) / gj[nohwarm]) - 2, $
				/data, thick=defthick, hsize = arrowsize
		endif

		; Print the mass of H2 in M_sun

		dl = dllist(where(taglist eq fname))
		if n_elements(temp_warm) gt 0 and h2fluxarr[1] ne 0 then begin
			mass_warm = h2mass(h2fluxarr[1], temp_warm, dl, 1, /msun, /quiet) 
			if keyword_set(verbose) then begin
				print,''
				print,'Mass of warm H2 for '+fname+' [10^7 M_sun] = ',mass_warm
			endif
			h2.mass_warm[j] = mass_warm
		endif

		dl = dllist(where(taglist eq fname))
		if n_elements(temp_hot) gt 0 and h2fluxarr[3] ne 0 then begin
			if temp_hot gt 0 then begin
				mass_hot = h2mass(h2fluxarr[3], temp_hot, dl, 3, /msun, /quiet) 
				if keyword_set(verbose) then begin
					print,''
					print,'Mass of hot H2 '+fname+' [10^7 M_sun] = ',mass_hot
				endif
				h2.mass_hot[j] = mass_hot
			endif
		endif

;		print,'Completed H2TEMP_CSO for '+fname

		multiplot

endfor

multiplot,/reset

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

h2.warm_avg = mean(h2.temp_warm[where(h2.temp_warm ne 0.)])

save, h2, filename='~/Astronomy/Research/Spitzer/cso/h2_cso.sav'

end
