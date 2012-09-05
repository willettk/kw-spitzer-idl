pro lockettplot_dusty, ps=ps, stop = stop, first = first
;+
; NAME:
;       
;	LOCKETTPLOT_DUSTY
;
; PURPOSE:
;
;	Plot dust temperature vs. optical depth in V to contrast our data w/predictions of maser strength
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
;	PS - create hard copy
;
; EXAMPLE:
;
;	IDL> lockettplot_dusty
;
; NOTES:
;
;	Follows plot in Lockett & Elitzur (2008), Figure 2
;
; REVISION HISTORY
;       Written by K. Willett              Jan 10  
;-

device, decomposed = 0

; OH data

a1667 = float(archdat('f1667'))
a1420 = float(archdat('f1420'))
atauoh = -1d * alog((a1667 + a1420)/a1420)

o1667 = float(ohmdat('f1667'))
o1420 = float(ohmdat('f1420'))
otauoh = -1d * alog((o1667 + o1420)/o1420)

ohm_tauoh = [transpose(atauoh), transpose(otauoh)]

c1667 = float(condat('f1667'))
c1420 = float(condat('f1420'))
con_tauoh = -1d * alog((c1667 + c1420)/c1420)

; DUSTY parameter fits

restore, '~/Astronomy/Research/Spitzer/dusty_pahfit_results.sav'	; [arch, mega]

ohmobj = [transpose(archdat('obj')), transpose(ohmdat('obj'))]
conobj = [transpose(condat('obj'))]

; Cull for objects with bad data in one or more parameters

ohm_goodind = where(finite(ohm_tauoh) eq 1 and ohm_tauoh ne 0 and ohm_tauvarr gt 0 and ohm_tdustarr gt 0)
con_goodind = where(finite(con_tauoh) eq 1 and con_tauvarr gt 0 and con_tdustarr gt 0)

ohm_tauoh = ohm_tauoh[ohm_goodind]
ohm_tauv = ohm_tauvarr[ohm_goodind]
ohm_tdust = ohm_tdustarr[ohm_goodind]
ohm_tags = ohmtags[ohm_goodind]
ohm_objs = ohmobj[ohm_goodind]

con_tauoh = con_tauoh[con_goodind]
con_tauv = con_tauvarr[con_goodind]
con_tdust = con_tdustarr[con_goodind]
con_tags = contags[con_goodind]
con_objs = conobj[con_goodind]

; Color indices

violet = fsc_color("Black")		; -4.0 < t_oh < -3.5
royalblue = fsc_color("Royal Blue")	; -3.5 < t_oh < -3.0
cyan = fsc_color("Cyan")		; -3.0 < t_oh < -2.5
seagreen = fsc_color("Sea Green")	; -2.5 < t_oh < -2.0
green = fsc_color("Green")		; -2.0 < t_oh < -1.5
yellow = fsc_color("Yellow")		; -1.5 < t_oh < -1.0
orange = fsc_color("Orange")		; -1.0 < t_oh < -0.5
defcolor = fsc_color("White")
white = fsc_color("White")
black = fsc_color("Black")
red = fsc_color("Red")

royalblueind = where(ohm_tauoh gt -3.5 and ohm_tauoh lt -3.0)
cyanind = where(ohm_tauoh gt -3.0 and ohm_tauoh lt -2.5)
seagreenind = where(ohm_tauoh gt -2.5 and ohm_tauoh lt -2.0)
greenind = where(ohm_tauoh gt -2.0 and ohm_tauoh lt -1.5)
yellowind = where(ohm_tauoh gt -1.5 and ohm_tauoh lt -1.0)
orangeind = where(ohm_tauoh gt -1.0 and ohm_tauoh lt -0.5)
defind = where(ohm_tauoh gt -0.5 and ohm_tauoh lt -0.0)

croyalblueind = where(con_tauoh gt -3.5 and con_tauoh lt -3.0)
ccyanind = where(con_tauoh gt -3.0 and con_tauoh lt -2.5)
cseagreenind = where(con_tauoh gt -2.5 and con_tauoh lt -2.0)
cgreenind = where(con_tauoh gt -2.0 and con_tauoh lt -1.5)
cyellowind = where(con_tauoh gt -1.5 and con_tauoh lt -1.0)
corangeind = where(con_tauoh gt -1.0 and con_tauoh lt -0.5)
cdefind = where(con_tauoh gt -0.5 and con_tauoh lt -0.0)

colors=[violet,royalblue,cyan,seagreen,green,yellow,orange, defcolor, defcolor]
barcolors = reverse(colors[0:n_elements(colors)-2])

; Load in contours from Lockett & Elitzur 2008

lockfile = '~/Astronomy/Research/Spitzer/ohm/lockett.txt'
readcol, lockfile, tdust, tauv, taumaser, format='i,i,f', skipline=1, /silent
tdust_arr = tdust(rem_dup(tdust))
tauv_arr = tauv(rem_dup(tauv))
taumaser_arr = transpose(reform(taumaser,10,9))

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/lockettplot_dusty.ps', /portrait, /color, xs=16, ys=16
	defcolor=fsc_color("Black")
	cs = 1.5
	barcs = 1.0
	th = 4
endif else begin
	defcolor=fsc_color("White")
	cs = 2
	barcs = 1.0
	th = 1
endelse


;loadct, 3, /silent
contour, taumaser_arr, tdust_arr, tauv_arr,$
	levels = fillarr(0.5,-4,0), $
	path_xy=xy, path_info=info, /path_data_coords

contour, taumaser_arr, tdust_arr, tauv_arr,$
	levels = fillarr(0.5,-4,0), $
	c_colors = colors, $
	/fill, $
	xtitle = 'T!Idust!N [K]', $
	ytitle = '!7s!3!IV!N', $
	xrange = [30,150], /xstyle, $
	yrange = [0,450], ystyle=9, $
	charsize = cs, $
	charthick = th, $
	thick = th, $
	xthick = th, $
	ythick = th, $
	position=[0.15, 0.15, 0.85, 0.8], $
	color=defcolor


; Screw it. I'm making my own colorbar. 

	for i=0,7 do polyfill, 0.15 + 0.7*[i, i, i+1, i+1, i] / 8., $
		[0.88, 0.95, 0.95, 0.88, 0.88], $
		/norm, color = barcolors[i]
	plots, [0.15, 0.15, 0.85, 0.85, 0.15], [0.88, 0.95, 0.95, 0.88, 0.88], /norm, color = defcolor, thick = th
	for i = 0,8 do begin
		plots, [0.15 + 0.7*i/8, 0.15+0.7*i/8], [0.88, 0.95], /norm, color=defcolor, thick = th
		xyouts, 0.13 + 0.7*i/8, 0.85, string(0. - 0.5*i,format='(f4.1)'), /norm, $
			color=defcolor, charsize=barcs, charthick = th
	endfor
	xyouts, 0.87, 0.90, '!7s!3(1667)', /normal, color=defcolor, charsize=barcs, charthick=th

if defind(0) ne -1 then oplot,ohm_tdust(defind), ohm_tauv(defind), psym=symcat(14), color=black, symsize=1.3
if royalblueind(0) ne -1 then oplot,ohm_tdust(royalblueind), ohm_tauv(royalblueind), psym=symcat(14), color=black, symsize=1.3
if cyanind(0) ne -1 then oplot,ohm_tdust(cyanind), ohm_tauv(cyanind), psym=symcat(14), color=black, symsize=1.3
if seagreenind(0) ne -1 then oplot,ohm_tdust(seagreenind), ohm_tauv(seagreenind), psym=symcat(14), color=black, symsize=1.3
if greenind(0) ne -1 then oplot,ohm_tdust(greenind), ohm_tauv(greenind), psym=symcat(14), color=black, symsize=1.3
if yellowind(0) ne -1 then oplot,ohm_tdust(yellowind), ohm_tauv(yellowind), psym=symcat(14), color=black, symsize=1.3
if orangeind(0) ne -1 then oplot,ohm_tdust(orangeind), ohm_tauv(orangeind), psym=symcat(14), color=black, symsize=1.3

if royalblueind(0) ne -1 then oplot,ohm_tdust(royalblueind), ohm_tauv(royalblueind), psym=symcat(14), color=royalblue, symsize=0.8
if cyanind(0) ne -1 then oplot,ohm_tdust(cyanind), ohm_tauv(cyanind), psym=symcat(14), color=cyan, symsize=1.0
if seagreenind(0) ne -1 then oplot,ohm_tdust(seagreenind), ohm_tauv(seagreenind), psym=symcat(14), color=seagreen, symsize=1.0
if greenind(0) ne -1 then oplot,ohm_tdust(greenind), ohm_tauv(greenind), psym=symcat(14), color=green, symsize=1.0
if yellowind(0) ne -1 then oplot,ohm_tdust(yellowind), ohm_tauv(yellowind), psym=symcat(14), color=yellow, symsize=1.0
if orangeind(0) ne -1 then oplot,ohm_tdust(orangeind), ohm_tauv(orangeind), psym=symcat(14), color=orange, symsize=1.0
if defind(0) ne -1 then oplot,ohm_tdust(defind), ohm_tauv(defind), psym=symcat(14), color=white, symsize=1.0

oplot, con_tdust, con_tauv, psym=symcat(7), color=black, thick=5, symsize=1.3

; Deal with overlapping points due to coarse grid

ohm_sum = ohm_tdust + ohm_tauv

for i=0, n_elements(ohm_tdust) - 1 do begin
	
	ind0 = where(ohm_tdust eq ohm_tdust[i] and ohm_tauv eq ohm_tauv[i], count)

;	if count gt 1 then print, ohm_tdust[i], ohm_tauv[i]

endfor


; Plausible physical limits for ULIRGs

; Yun and Carilli find max of 75K in their sample of 23 dusty starburst galaxies; for galaxies w/reliable IRAS fluxes, none above 100K
; Max. optical depth detected in IRS data has S_sil = -4.0 -> tau_v = 68

;hor, 68, linestyle=2
;ver, 75, linestyle=2

; Average error bar for the sample

oploterror, [127], [403], [20], [20], thick=3

axis, yaxis=1, yrange=[0,450] * 1.086, /save, color = defcolor, $
	ytitle='A!IV!N', $
	charsize = cs, charthick = th

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Testing on what contour a point lies

tauoh_predicted = fltarr(n_elements(ohm_tdust))

for k=0, n_elements(ohm_tdust) - 1 do begin

;	plot,indgen(10),/nodata,xr=[30,180],yr=[-20,320]
	
	tlevels = fillarr(0.5,-4,0)
	wlevels = intarr(n_elements(tlevels))
	
	for i=0, n_elements(tlevels)-1 do begin
		ncont = where(info.value eq tlevels[i])
		temp_wlevel = 0
		for j=0, n_elements(ncont) - 1 do begin
			rb_contour = xy(*,info[ncont[j]].offset+[indgen(info[ncont[j]].n),0])
			if tlevels[i] le -2.0 then begin
				rb_ind = where(rb_contour(0,*) ne 30. and rb_contour(1,*) ne 300. and rb_contour(1,*) ne 0.)
				rb_contour = rb_contour[*,rb_ind]
			endif
			rx = transpose(rb_contour[0,*]) & ry = transpose(rb_contour[1,*])
;			plots,rb_contour, color=colors[i]
;			xyouts,10 + i*15, 350, string(tlevels[i],format='(f4.1)'),/data,color=colors[i]
		
			; Test if point is within this contour
		
			x = ohm_tdust[k]
			y = ohm_tauv[k]
	
;			plots,[x],[y],psym=symcat(14),color=red
		
			temp_wlevel = temp_wlevel + inside(x,y,rx,ry)
	
		endfor
	
			wlevels[i] = temp_wlevel
	
	endfor
	
	wl = where(wlevels eq 1)
	if wl[0] ne -1 then begin
;		print,tlevels[where(wlevels eq 1)]
		if min(tlevels[wl]) le -2.0 then tpred = min(tlevels[wl]) else tpred = max(tlevels[wl])
	endif else tpred = -1.5

	tauoh_predicted[k] = tpred
endfor

;if not keyword_set(first) then begin
;
;	if keyword_set(ps) then begin
;		set_plot,'ps'
;		device, filename='~/Astronomy/Research/Spitzer/papers/locketthist.ps', /portrait, /color
;	endif
;	
;	plothist, (tauoh_predicted - ohm_tauoh), bin = 0.5, $
;		xtitle = textoidl('\tau^{OH}_{pred} - \tau^{OH}_{meas}'), $
;		ytitle = 'Count', $
;		charsize = cs, $
;		thick = th, $
;		xthick = th, $
;		ythick = th, $
;		charthick = th
;	
;	ver, mean(tauoh_predicted - ohm_tauoh), linestyle = 1, thick = th
;	
;	if keyword_set(ps) then begin
;		device,/close
;		set_plot,'x'
;	endif
;
;endif

;; Perform chi-squared goodness of fit test
;
;ohms_1420 = [o1420[ogoodind], a1420[agoodind]]
;ohms_1667 = [o1667[ogoodind], a1667[agoodind]]
;nohm = n_elements(ohms_1420)
;
;; Assume standard error of 0.1 mJy for the continuum flux (FIRST survey), 0.01 mJy for OHM peak
;
;sig_1420 = replicate(0.1, nohm)
;sig_1667 = replicate(0.01, nohm)
;
;dtau_d1420 = ohms_1667 / ((ohms_1667 + ohms_1420) * ohms_1420)
;dtau_d1667 = -1d / (ohms_1667 + ohms_1420)
;
;smean_1420 = mean(ohms_1420)
;smean_1667 = mean(ohms_1667)
;
;; Calculate sample covariance between f1420 and f1667 (Bevington, Eqn. 11.21)
;
;cov_1420_1667 = 1d / (nohm - 1) * total((ohms_1420 - smean_1420) * (ohms_1667 - smean_1667))
;
;tau_err     = sqrt(sig_1420^2 * (dtau_d1420)^2 + sig_1667^2 * (dtau_d1667)^2)
;tau_err_cov = sqrt(sig_1420^2 * (dtau_d1420)^2 + sig_1667^2 * (dtau_d1667)^2 + 2d * cov_1420_1667 * (dtau_d1420) * (dtau_d1667))
;
;tau_err = replicate(0.25,nohm)
;
;chisq = total((tauoh_predicted - ohm_tauoh)^2 / tau_err^2)
;reduced_chisq = chisq / (nohm - 2)
;
;print,''
;print, 'Reduced chi^2 = ',reduced_chisq
;print,''

if keyword_set(stop) then stop
end
