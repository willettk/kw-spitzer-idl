pro feature2_ohm, ps = ps, stop = stop, label=label

;+
; NAME:
;       
;	FEATURE2_OHM
;
; PURPOSE:
;
;	Plot feature-feature diagram for the silicate 10 and 18 um features in the CSO paper
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
;	LABEL - label data points
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
;       Based off FEATURE2.pro - Jan 10
;	Removed plots for p=1, 2 (not well fit by DUSTY) - Mar 2010
;-

ohmjunk = ohmdat('sil')
ohmerrjunk = ohmdat('silerr')
ohmobjjunk = ohmdat('obj')

archjunk = archdat('sil')
archerrjunk = archdat('silerr')
archobjjunk = archdat('obj')

sil10_ohm = ohmjunk[0,*]
sil18_ohm = ohmjunk[1,*]
sil10_ohm_err = abs(ohmerrjunk[0,*])
sil18_ohm_err = abs(ohmerrjunk[1,*])

sil10_arch = archjunk[0,*]
sil18_arch = archjunk[1,*]
sil10_arch_err = abs(archerrjunk[0,*])
sil18_arch_err = abs(archerrjunk[1,*])

conjunk = condat('sil')
conerrjunk = condat('silerr')
conobjjunk = condat('obj')

sil10_con = conjunk[0,*]
sil18_con = conjunk[1,*]
sil10_con_err = abs(conerrjunk[0,*])
sil18_con_err = abs(conerrjunk[1,*])

; Remove control033, mega034

ohmind = where(strtrim(ohmobjjunk,2) ne 'IRAS 23028+0725')
conind = where(strtrim(conobjjunk,2) ne 'IRAS 20460+1925')

sil10_ohm = sil10_ohm[ohmind]
sil18_ohm = sil18_ohm[ohmind]
sil10_ohm_err = sil10_ohm_err[ohmind]
sil18_ohm_err = sil18_ohm_err[ohmind]
ohmobj = ohmobjjunk[ohmind]

sil10_con = sil10_con[conind]
sil18_con = sil18_con[conind]
sil10_con_err = sil10_con_err[conind]
sil18_con_err = sil18_con_err[conind]
conobj = conobjjunk[conind]

archobj = archobjjunk

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/feature2_ohm.ps',/color,/portrait, xs=20, ys=16
	defcolor = fsc_color("Black")
	cs = 1.8
	lthick = 5
	cthick = 5
endif else begin
	defcolor=fsc_color("White")
	cs = 1.5
	lthick = 1
	cthick = 1
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
yellow = fsc_color("Goldenrod")
;green = fsc_color("Dark Green")
green = fsc_color("Green")

plot, indgen(10), $
	/nodata, $
	;xtitle='S!I9.7!N', $
	;ytitle='S!I18!N', $
	xtitle='Silicate strength at 10 !7l!3m', $
	ytitle='Silicate strength at 18 !7l!3m', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charthick = cthick, $
	charsize=cs, $
	color=defcolor, $
	xr=[-4,1.5], /xstyle, $
	yr=[-2.0,1.0], /ystyle

xyouts, /data, 0.0, -0.5, 'Clumpy', color=green, charsize=cs, charthick=cthick
xyouts, /data, -2.8, 0.2, 'Smooth', color=defcolor, charsize=cs, charthick=cthick

; DUSTY tracks from Sirocky 2008, measured using Dexter

xstart = 1.26 & ystart = 0.67

dexdir='/Applications/Dexter/'
readcol,dexdir+'f7.gif.yellow',skipline=1,x,y,/silent
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=yellow, thick = lthick

readcol,dexdir+'f7.gif.green1',skipline=1,x,y,/silent
plots, [x,xstart], [y,ystart], color=green, thick = lthick
readcol,dexdir+'f7.gif.green2',skipline=1,x,y,/silent
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=green, thick = lthick
readcol,dexdir+'f7.gif.green3',skipline=1,x,y,/silent
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=green, thick = lthick

readcol,dexdir+'f7.gif.blue2',skipline=1,x,y,/silent
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick

readcol,dexdir+'f7.gif.black_dash2',skipline=1,x,y,/silent
	black_dash_ind = where(x gt -4 and y gt -1)
	x = x[black_dash_ind] & y = y[black_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=2
readcol,dexdir+'f7.gif.black_dot2',skipline=1,x,y,/silent
	black_dash_ind = where(x gt -4 and y gt -1)
	x = x[black_dash_ind] & y = y[black_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=1
readcol,dexdir+'f7.gif.black_solid2',skipline=1,x,y,/silent
	black_dash_ind = where(x gt -4 and y gt -1)
	x = x[black_dash_ind] & y = y[black_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=0

readcol,dexdir+'f7.gif.red_dash2',skipline=1,x,y,/silent
	red_dash_ind = where(x gt -4 and y gt -1)
	x = x[red_dash_ind] & y = y[red_dash_ind]
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=2
readcol,dexdir+'f7.gif.red_dot2',skipline=1,x,y,/silent
	red_dot_ind = where(x gt -4 and y gt -1)
	x = x[red_dot_ind] & y = y[red_dot_ind]
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=1
readcol,dexdir+'f7.gif.red_solid2',skipline=1,x,y,/silent
	red_solid_ind = where(x gt -4 and y gt -1)
	x = x[red_solid_ind] & y = y[red_solid_ind]
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=0

; ULIRG data

;oplot, sil10_ohm, sil18_ohm, psym=symcat(14), color=red
;oplot, sil10_arch, sil18_arch, psym=symcat(14), color=red
;oplot, sil10_con, sil18_con, psym=symcat(14), color=red

oploterror, sil10_ohm, sil18_ohm, sil10_ohm_err, sil18_ohm_err, $
	psym=symcat(14), color=red, errthick = lthick, errcolor=red, /nohat
oploterror, sil10_arch, sil18_arch, sil10_arch_err, sil18_arch_err, $
	psym=symcat(14), color=red, errthick = lthick, errcolor=red, /nohat
oploterror, sil10_con, sil18_con, sil10_con_err, sil18_con_err, $
	psym=symcat(7), color=blue, errthick = lthick, errcolor=blue, /nohat, thick=lthick, symsize=1.2

; CSO data

;oploterror, sil10_cso, sil18_cso, sil10_cso_err, sil18_cso_err, $
;	psym=symcat(16), color=defcolor, errthick = 4, errcolor=defcolor, /nohat

; Starting point

oplot,[xstart],[ystart], psym=symcat(45), symsize=4.5
oplot,[xstart],[ystart], psym=symcat(45), symsize=4.25
oplot,[xstart],[ystart], psym=symcat(45), symsize=4
oplot,[xstart],[ystart], psym=symcat(45), symsize=3.75
oplot,[xstart],[ystart], psym=symcat(45), symsize=3.5

;legend,['Y=100', 'Y=200','Y=400','clumpy'],linestyle=[0,1,2,0], $
;	color=[defcolor,defcolor,defcolor,green], $ 
;	thick = lthick, charthick = cthick, /top, /left

;legend,['CSOs','ULIRGs'],psym=[16,14],color=[defcolor,red], $ 
;	thick = lthick, charthick = cthick, /bottom, /right

if keyword_set(label) then begin
	xyouts, sil10_ohm, sil18_ohm, ohmobj, charsize=1, color=fsc_color("Orange"), charthick=1
	xyouts, sil10_arch, sil18_arch, archobj, charsize=1, color=fsc_color("Orange"), charthick=1
	xyouts, sil10_con, sil18_con, conobj, charsize=1, color=fsc_color("Orange"), charthick=1
endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
