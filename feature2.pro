pro feature2, ps = ps, stop = stop, label=label

;+
; NAME:
;       
;	FEATURE2
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
;       Written by K. Willett                Jul 09
;	Added more accurate tracks using Dexter - Aug 09
;	Added LABEL keyword - Oct 09
;	Changed label of power-law index from p to q - Feb 10
;-

csojunk = csodat('sil')
csoerrjunk = csodat('silerr')
csoobjjunk = csodat('obj')

sil10_cso = csojunk[0,*]
sil18_cso = csojunk[1,*]
sil10_cso_err = abs(csoerrjunk[0,*])
sil18_cso_err = abs(csoerrjunk[1,*])

; Remove VII Zw 485, NGC 5793

csoind = [0,1,3,4,5,6,7,8]

sil10_cso = sil10_cso[csoind]
sil18_cso = sil18_cso[csoind]
sil10_cso_err = sil10_cso_err[csoind]
sil18_cso_err = sil18_cso_err[csoind]
csoobj=csoobjjunk[csoind]

; ohmjunk = ohmdat('sil')
; 
; sil10_ohm = ohmjunk[0,*]
; sil18_ohm = ohmjunk[1,*]
; 
; archjunk = archdat('sil')
; 
; sil10_arch = archjunk[0,*]
; sil18_arch = archjunk[1,*]
; 
; conjunk = condat('sil')
; 
; sil10_con = conjunk[0,*]
; sil18_con = conjunk[1,*]

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/feature2.ps',/color,/portrait
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
green = fsc_color("Dark Green")

plot, indgen(10), $
	/nodata, $
	xtitle='S!I9.7!N', $
	ytitle='S!I18!N', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charthick = cthick, $
	charsize=cs, $
	color=defcolor, $
	xr=[-3,1.5], /xstyle, $
	yr=[-1,1.0]

; DUSTY tracks from Sirocky 2008, measured using Dexter

xstart = 1.26 & ystart = 0.67

dexdir='/Applications/Dexter/'
readcol,dexdir+'f7.gif.yellow',skipline=1,x,y,/silent
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=yellow, thick = 6

readcol,dexdir+'f7.gif.green1',skipline=1,x,y,/silent
plots, [x,xstart], [y,ystart], color=green, thick = 6
readcol,dexdir+'f7.gif.green2',skipline=1,x,y,/silent
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=green, thick = 6
readcol,dexdir+'f7.gif.green3',skipline=1,x,y,/silent
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=green, thick = 6

readcol,dexdir+'f7.gif.blue',skipline=1,x,y,/silent
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=blue, thick = 6

readcol,dexdir+'f7.gif.black_dash',skipline=1,x,y,/silent
	black_dash_ind = where(x gt -3 and y gt -1)
	x = x[black_dash_ind] & y = y[black_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = 6, linestyle=2
readcol,dexdir+'f7.gif.black_dot',skipline=1,x,y,/silent
	black_dash_ind = where(x gt -3 and y gt -1)
	x = x[black_dash_ind] & y = y[black_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = 6, linestyle=1
readcol,dexdir+'f7.gif.black_solid',skipline=1,x,y,/silent
	black_dash_ind = where(x gt -3 and y gt -1)
	x = x[black_dash_ind] & y = y[black_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = 6, linestyle=0

readcol,dexdir+'f7.gif.red_dash',skipline=1,x,y,/silent
	red_dash_ind = where(x gt -3 and y gt -1)
	x = x[red_dash_ind] & y = y[red_dash_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=red, thick = 6, linestyle=2
readcol,dexdir+'f7.gif.red_dot',skipline=1,x,y,/silent
	red_dot_ind = where(x gt -3 and y gt -1)
	x = x[red_dot_ind] & y = y[red_dot_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=red, thick = 6, linestyle=1
readcol,dexdir+'f7.gif.red_solid',skipline=1,x,y,/silent
	red_solid_ind = where(x gt -3 and y gt -1)
	x = x[red_solid_ind] & y = y[red_solid_ind]
plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=red, thick = 6, linestyle=0

; ULIRG data

;oplot, sil10_ohm, sil18_ohm, psym=symcat(14), color=red
;oplot, sil10_arch, sil18_arch, psym=symcat(14), color=red
;oplot, sil10_con, sil18_con, psym=symcat(14), color=red

; CSO data

oploterror, sil10_cso, sil18_cso, sil10_cso_err, sil18_cso_err, $
	psym=symcat(16), color=defcolor, errthick = 4, errcolor=defcolor, /nohat

; Starting point

oplot,[xstart],[ystart], psym=symcat(46), symsize=3

legend,['Y=100', 'Y=200','Y=400','q=0','q=1','q=2','slab','clumps'],linestyle=[0,1,2,0,0,0,0,0], $
	color=[defcolor,defcolor,defcolor,defcolor,red,blue,yellow,green], $ 
	thick = lthick, charthick = cthick, /top, /left

;legend,['CSOs','ULIRGs'],psym=[16,14],color=[defcolor,red], $ 
;	thick = lthick, charthick = cthick, /bottom, /right

if keyword_set(label) then begin
	xyouts, sil10_cso, sil18_cso, csoobj, charsize=1, color=fsc_color("Orange"), charthick=2
endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
