pro h2temp_2x2, verbose = verbose, ps = ps, silent = silent
;+
; NAME:
;       
;	H2TEMP_2x2
;
; PURPOSE:
;
;	Create 2x2 diagram of example H2 diagrams for Paper II
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
;	IDL> h2temp_2x2, /ps
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett			Nov 08
;-

if not keyword_set(verbose) then quiet = 1

; Set plot type

!p.multi=[0,2,2]

; Data lists (both MEGA and ARCH OHMs)

ohmtag = ohmdat('tag')
ohmobj = ohmdat('obj')
ohmdl = ohmdat('dl')

archtag = archdat('tag')
archobj = archdat('obj')
archdl = archdat('dl')

contag = condat('tag')
conobj = condat('obj')
condl = condat('dl')

taglist = ['control037','mega010','arch007','mega016']
objlist = [conobj[where(contag eq taglist[0])], ohmobj[where(ohmtag eq taglist[1])], $
	archobj[where(archtag eq taglist[2])], ohmobj[where(ohmtag eq taglist[3])]] 
dllist = [condl[where(contag eq taglist[0])], ohmdl[where(ohmtag eq taglist[1])], $
	archdl[where(archtag eq taglist[2])], ohmdl[where(ohmtag eq taglist[3])]] 

nobj = n_elements(objlist)

h2 = {obj:strarr(nobj), tag:strarr(nobj), $
	temp_warm:fltarr(nobj), temp_hot:fltarr(nobj), $
	mass_warm:fltarr(nobj), mass_hot:fltarr(nobj), $
	warm_avg:0d, hot_avg:0d}
 
; Physical constants for H2 transitions taken from 						
; http://www.not.iac.es/instruments/notcam/ReferenceInfo/h2_lines.html

h2_indices = [0,1,2,3,5,7]					; All H2 transitions seen in sample(s)
ej = [510, 1015, 1682, 2504, 4586, 7197]			; Energy of upper state
aj = [3d-4, 4.8d-3, 2.76d-2, 9.84d-2, 0.588, 2.00]		; Einstein-A coefficient
lambda = [28.221,17.035,12.279,9.6649,6.9091,5.5115]		; Rest wavelength of transition
gj = [5,21,9,33,45,57]						; Statistical weight of upper level

plancks = 6.6261d-27					; Planck's constant
c = 3d10						; Speed of light [cm/s]
del_ej = plancks * c / lambda				; Energy of transition

; Linear fit for systems with a single excitation temperature

expr = 'p[0] + x*p[1]'
start = [10, 3d-3]

	h2filename = '~/Astronomy/Research/Spitzer/papers/h2temp_2x2.ps'

	if keyword_set(ps) then begin
		set_plot,'ps'
		device, filename=h2filename, /portrait
		tempcs = 1.0
		axescs = 0.8
		labelcs = 0.9
		defthick = 3
		tempthick = 3
		cthick = 3
		psymsize = 0.6
		hsize = 300
	endif else begin
		set_plot,'x'
		tempcs = 1
		axescs = 1
		labelcs = 1
		defthick = 1
		tempthick = 1
		cthick = 1
		psymsize = 1
		hsize = 10
	endelse
	
	for j=0,3 do begin
	
		kind = j
	
		fname = taglist[kind]
		obj = objlist(where(taglist eq fname))
	
		h2.tag[kind] = fname
		h2.obj[kind] = obj
	
		temp_hot = 0.
	
		; Load in the line fluxes and errors
		
		orthoh2 = ['h2s0','h2s1','h2s2','h2s3',$
			'arIIh2s5_lr','h2s7_lr']
		h2fluxarr = dblarr(n_elements(orthoh2))
		h2errarr  = dblarr(n_elements(orthoh2))
		
		for i = 0, n_elements(h2fluxarr)-1 do begin
			
			h2_query = getlineflux(fname, orthoh2[i],quiet=quiet)
			if h2_query[0] ne -1 then begin
				if keyword_set(verbose) then print,strupcase(orthoh2[i])+' line found for '+fname
				h2fluxarr[i] = h2_query[0]
				h2errarr[i] = h2_query[1]
			endif
		
		endfor
		
		; Plot the excitation diagrams					
		
		hlabels = ['S(0)', 'S(1)', 'S(2)', 'S(3)', 'S(5)', 'S(7)']
		
		if n_elements(where(h2fluxarr ne 0)) gt 1 then begin
	
			; Plots where multiple lines are detected and at least one temperature can be fit 
	
			plot, ej, alog(gj), $
				/nodata, $
				xtitle = 'E!IJ!N/k!Ib!N [K]', $
				ytitle = 'ln(N!IJ!N / g!IJ!N)', $
				charsize = axescs, $
				charthick = cthick, $
				thick = defthick, $
				xthick = defthick, $
				ythick = defthick, $
				title = obj, $
				xr = [0,9000], $
				yr = [-20,1], /ystyle
	
			; Plot points (if any) from S(4) through S(7)
	
			hhot = where(h2fluxarr[4:5] ne 0) + 4
	
			if hhot[0] gt 3 then begin
	
				nj_hot		 = h2fluxarr[hhot] / (aj[hhot] * del_ej[hhot])
				njerr_hot	 = h2errarr[hhot]  / (aj[hhot] * del_ej[hhot])
		
				oplot, ej[hhot], alog(nj_hot / gj[hhot]), psym = symcat(16), symsize = psymsize
				xyouts, ej[hhot]+250, alog(nj_hot / gj[hhot]), hlabels[hhot], charsize = labelcs, charthick = cthick
		
				if h2fluxarr[5] ne 0 then begin
					nj5    = h2fluxarr[5] / (aj[5] * del_ej[5])
					njerr5 = h2errarr[5]  / (aj[5] * del_ej[5])
	
					if nj5 lt njerr5 then begin
						print,'S(5) error is larger than flux for '+fname
						njerr5 = 0.5 * nj5
					endif
					oploterror, ej[5], alog(nj5 / gj[5]), $
						alog(nj5 / gj[5]) - alog((nj5 - njerr5) / gj[5]), $
						psym = 3, /lobar, thick = 2
					oploterror, ej[5], alog(nj5 / gj[5]), $
						alog(nj5 / gj[5]) - alog((nj5 + njerr5) / gj[5]), $
						psym = 3, /hibar, thick = 2
				endif
		
			endif
	
			; For S(5), plot an upper limit
	
			if h2fluxarr[4] ne 0 then $
				arrow, ej[4], alog((h2fluxarr[4]/ (aj[4] * del_ej[4])) / gj[4]), $
					ej[4], alog((h2fluxarr[4]/ (aj[4] * del_ej[4])) / gj[4]) - 4, $
					/data, thick=defthick, hsize = hsize
	
			; If H2 S(7) exists, fit a temperature between S(3) and S(7)
	
			if h2fluxarr[3] ne 0 and h2fluxarr[5] ne 0 then begin
				hotind = [3,5]
				nj_hotfit	 = h2fluxarr[hotind] / (aj[hotind] * del_ej[hotind])
				njerr_hotfit	 = h2errarr[hotind]  / (aj[hotind] * del_ej[hotind])
	
				hot_result = mpfitexpr(expr, ej[hotind], alog(nj_hotfit / gj[hotind]), $
					alog(njerr_hotfit / gj[hotind]), start, /quiet)
				hotarr = fillarr(1,ej[3],ej[5])
				oplot, hotarr, (hot_result[0] + hotarr*hot_result[1]), $
					linestyle = 2, thick = tempthick
				temp_hot = -1d / hot_result[1]
				xyouts, 5200, -9, 'T!Ihot!N = '+string(temp_hot, format = '(i4)')+' K', $
					charsize = tempcs, charthick = cthick
				h2.temp_hot[kind] = temp_hot
			endif
	
			; Fit warm lines from S(0) - S(3)
	
			hwarm = where(h2fluxarr[0:3] ne 0)
	
			; If a hot gas component was detected, subtract expected flux contribution from the warm lines
	
			warmh2 = where(h2fluxarr[0:3] ne 0)
			if temp_hot ne 0 and warmh2[0] ne -1 then begin
				orig_flux = h2fluxarr[warmh2]
				ytemp = (hot_result[0] + hot_result[1] * ej[warmh2])
				ysub = aj[warmh2] * del_ej[warmh2] * gj[warmh2] * exp(ytemp)
				h2fluxarr[warmh2] = (h2fluxarr[warmh2] - ysub) + orig_flux[n_elements(orig_flux)-1]
				h2errarr[warmh2] = h2errarr[warmh2] * h2fluxarr[warmh2] / orig_flux	; Scale errors
			endif
	
			nj_warm     = h2fluxarr[hwarm] / (aj[hwarm] * del_ej[hwarm])
			njerr_warm = h2errarr[hwarm]  / (aj[hwarm] * del_ej[hwarm])
	
	
			oplot, ej[hwarm], alog(nj_warm / gj[hwarm]), psym = symcat(16), symsize = psymsize
			xyouts, ej[hwarm]+250, alog(nj_warm / gj[hwarm]), hlabels[hwarm], charsize = labelcs, charthick = cthick
	
			if fname eq 'arch030' then begin

				ej4 = 3474d
				aj4 = 0.264d
				l4  = 8.0258d
				gj4 = 13d

				h2s4arr = getlineflux('arch030','h2s4')
				h2s4 = h2s4arr[0]
				h2s4err = h2s4arr[1]

				delej4 = plancks * c / l4
				nj4 = h2s4 / (aj4 * delej4)
				nj4err = h2s4err / (aj4 * delej4)

				ej4 = [ej4] & nj4 = [nj4] & gj4 = [gj4] & nj4err = [nj4err]

				oplot, ej4, alog(nj4 / gj4), psym = symcat(16), symsize = psymsize
				xyouts, ej4+250, alog(nj4 / gj4), 'S(4)', $
					charsize = labelcs, charthick = cthick

				oploterror, ej4, alog(nj4 / gj4), alog(nj4 / gj4) $
					- alog((nj4 - nj4err) / gj[hwarm]), psym = 3, /lobar
				oploterror, ej4, alog(nj4 / gj4), alog(nj4 / gj4) $
					- alog((nj4 + nj4err) / gj4), psym = 3, /hibar

			endif

			oploterror, ej[hwarm], alog(nj_warm / gj[hwarm]), alog(nj_warm / gj[hwarm]) $
				- alog((nj_warm - njerr_warm) / gj[hwarm]), psym = 3, /lobar
			oploterror, ej[hwarm], alog(nj_warm / gj[hwarm]), alog(nj_warm / gj[hwarm]) $
				- alog((nj_warm + njerr_warm) / gj[hwarm]), psym = 3, /hibar
			
			if n_elements(hwarm) gt 1 then begin
				warm_result = mpfitexpr(expr, ej[hwarm], alog(nj_warm / gj[hwarm]), $
					alog(njerr_warm / gj[hwarm]), start, /quiet)
				warmarr = fillarr(1,ej[hwarm[0]],ej[hwarm[n_elements(hwarm)-1]])
				oplot, warmarr, (warm_result[0] + warmarr*warm_result[1]), linestyle = 1, thick = tempthick
				temp_warm = -1d / warm_result[1]
				xyouts, 3000, -5, 'T!Iwarm!N = '+string(temp_warm, format = '(i3)')+' K', $
					charsize = tempcs, charthick = cthick
	
				h2.temp_warm[kind] = temp_warm
			endif
	
			; Print the mass of H2 in M_sun
	
			dl = dllist(where(taglist eq fname))
			if n_elements(temp_warm) gt 0 and h2fluxarr[1] ne 0 then begin
				mass_warm = h2mass(h2fluxarr[1], temp_warm, dl, 1, /msun, quiet=quiet) 
				if keyword_set(verbose) then begin
					print,''
					print,'Mass of warm H2 [10^7 M_sun] = ',mass_warm
				endif
				h2.mass_warm[kind] = mass_warm
			endif
	
			dl = dllist(where(taglist eq fname))
			if n_elements(temp_hot) gt 0 and h2fluxarr[3] ne 0 then begin
				if temp_hot gt 0 then begin
					mass_hot = h2mass(h2fluxarr[3], temp_hot, dl, 3, /msun, quiet=quiet) 
					if keyword_set(verbose) then begin
						print,''
						print,'Mass of hot H2 [10^7 M_sun] = ',mass_hot
					endif
					h2.mass_hot[kind] = mass_hot
				endif
			endif
	
			print,'Completed H2TEMP for '+fname
	
		endif else if total(where(h2fluxarr ne 0)) eq 1 then begin
	
			; Plot diagrams where only a single measurement of H2 is made
	
			onlypoint = where(h2fluxarr ne 0)
			nj     = h2fluxarr[onlypoint] / (aj[onlypoint] * del_ej[onlypoint])
			nj_err = h2errarr[onlypoint]  / (aj[onlypoint] * del_ej[onlypoint])
	
			plot, ej[onlypoint], alog(nj / gj[onlypoint]), $
				xtitle = 'E!IJ!N/k!Ib!N [K]', $
				ytitle = 'ln(N!IJ!N / g!IJ!N)', $
				psym = symcat(16), $
				symsize = psymsize, $
				charsize = axescs, $
				charthick = cthick, $
				thick = defthick, $
				xthick = defthick, $
				ythick = defthick, $
				title = obj, $
				xr = [0,9000], $
				yr = [-20,0]
	
			oploterror, ej[onlypoint], alog(nj / gj[onlypoint]), alog(nj / gj[onlypoint]) $
				- alog((nj - nj_err) / gj[onlypoint]), psym = 3, /lobar, thick = 2
			oploterror, ej[onlypoint], alog(nj / gj[onlypoint]), -alog(nj / gj[onlypoint]) $
				+ alog((nj + nj_err) / gj[onlypoint]), psym = 3, /hibar, thick = 2
	
			xyouts, ej[onlypoint]+250, alog(nj / gj[onlypoint]), hlabels[onlypoint], $
				charsize = labelcs, charthick = cthick
	
			if not keyword_set(silent) then print,'Completed H2TEMP for '+fname
	
		endif else begin
	
			; No H2 lines found - do not plot 
	
			if not keyword_set(silent) then print,'No H2 lines found for '+fname
		endelse
	
	endfor
	
	if keyword_set(ps) then begin
		device,/close
		set_plot,'x'
	endif
	
h2.warm_avg = mean(h2.temp_warm[where(h2.temp_warm ne 0.)])
h2.hot_avg = mean(h2.temp_hot[where(h2.temp_hot ne 0.)])

end
