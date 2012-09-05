pro h2temp_con, verbose = verbose, ps = ps
;+
; NAME:
;       
;	H2TEMP_CON
;
; PURPOSE:
;
;	Compute excitation temperature of gas using H2 lines and create figures for publication
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
;	IDL> h2temp_con, /ps
;
; NOTES:
;
;
;
; REVISION HISTORY
; 	Adapted from H2TEMP.pro - Aug 08
;-

if not keyword_set(verbose) then quiet = 1

; Set plot type

!p.multi=[0,3,5]

; Data lists

taglist = condat('tag')
objlist = condat('obj')
dllist = condat('dl')

sortind = sort(objlist)

objlist = objlist[sortind]
taglist = taglist[sortind]
dllist  = dllist[sortind]

; Find which targets have at least one confirmed H2 line and trim the list

h2 = ['h2s0','h2s1','h2s2','h2s3','h2s7_lr']
nh2 = n_elements(h2)

ind = [0]	; Not perfect, since this assumes the zeroth element is present.

for j=0, nh2 - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/linedata/'+h2[j]+'.sav'
	match, taglist, line.tag, a, b
	ind = setunion(ind,a)

endfor

objlist = objlist[ind]
taglist = taglist[ind]
dllist = dllist[ind]

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

	filenames = ['~/Astronomy/Research/Spitzer/papers/h2temp_con.ps']
	
	if keyword_set(ps) then begin
		set_plot,'ps'
		device, filename=filenames[0], /portrait
		tempcs = 0.6
		axescs = 1
		labelcs = 0.5
		defthick = 3
		tempthick = 2
		cthick = 2
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
	
	for j=0,nobj-1 do begin
	
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
				charsize = axes, $
				charthick = cthick, $
				thick = defthick, $
				xthick = 2, $
				ythick = 2, $
				title = obj, $
				xr = [0,9000], $
				yr = [-20,0]
	
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
						psym = 3, /lobar
					oploterror, ej[5], alog(nj5 / gj[5]), $
						alog(nj5 / gj[5]) - alog((nj5 + njerr5) / gj[5]), $
						psym = 3, /hibar
				endif
		
			endif
	
			; For S(5), plot an upper limit
	
			if h2fluxarr[4] ne 0 then $
				arrow, ej[4], alog((h2fluxarr[4]/ (aj[4] * del_ej[4])) / gj[4]), $
					ej[4], alog((h2fluxarr[4]/ (aj[4] * del_ej[4])) / gj[4]) - 4, $
					/data, thick=defthick, hsize = arrowsize
	
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
				xyouts, 6500, -11, 'T!Ihot!N = '+string(temp_hot, format = '(i4)')+' K', $
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
	
			print,'Completed H2TEMP_CON for '+fname
	
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
				title = obj, $
				xr = [0,9000], $
				yr = [-20,0]
	
			oploterror, ej[onlypoint], alog(nj / gj[onlypoint]), alog(nj / gj[onlypoint]) $
				- alog((nj - nj_err) / gj[onlypoint]), psym = 3, /lobar
			oploterror, ej[onlypoint], alog(nj / gj[onlypoint]), -alog(nj / gj[onlypoint]) $
				+ alog((nj + nj_err) / gj[onlypoint]), psym = 3, /hibar
	
			xyouts, ej[onlypoint]+250, alog(nj / gj[onlypoint]), hlabels[onlypoint], $
				charsize = labelcs, charthick = cthick
	
			print,'Completed H2TEMP_CON for '+fname
	
		endif else begin
	
			; No H2 lines found - do not plot 
	
			print,'No H2 lines found for '+fname
		endelse
	
;		if k eq 1 and j eq 13 then j = 15		; Kluge to make two plots of 3x5 apiece
	endfor
	
	if keyword_set(ps) then begin
		device,/close
		set_plot,'x'
	endif
	
h2.warm_avg = mean(h2.temp_warm[where(h2.temp_warm ne 0.)])
h2.hot_avg = mean(h2.temp_hot[where(h2.temp_hot ne 0.)])

cd,'~/Astronomy/Research/Spitzer/control'
save, h2, filename='~/Astronomy/Research/Spitzer/control/h2_con.sav'

end
