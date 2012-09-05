;+
; NAME: 
;       OHSIL_POSTER
;
; PURPOSE:
;
;	Plot the peak OH luminosity vs. the 9.7 um silicate strength for OHMs and control sample galaxies
;
; CATEGORY:
;	ASTRONOMY
;
; KEYWORDS:
;
;	PS - create PS hard copy of plot
;
; EXAMPLE:
;
;	IDL> ohsil
;         
; MODIFICATION HISTORY:
;
;	Written by KW - Sep 07
;	Added archived OHMs - KW, Dec 07
;-
pro ohsil_poster, ps = ps, bw = bw


f = ohmdat('sil')
fsil = float(f(0,*))

c = condat('sil',sz=2)
csil = float(c(*,0))

a = archdat('sil')
asil = float(a(0,*))

oh = float(ohmdat('logoh'))
oh_con = float(condat('logoh'))
oh_arch = float(archdat('logoh'))


!p.multi=[0,1,1]

plotname='~/Desktop/jila/ohsil_jila.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color,/landscape
	cs = 1.7
	cthick = 3
	lthick = 3
	defcolor=fsc_color("Black")
endif else begin
	cs = 2
	cthick = 1
	lthick = 1
	defcolor=fsc_color("White")
endelse


red = fsc_color("Red")
blue = fsc_color("Blue")
green = fsc_color("Green")
if keyword_set(bw) then begin
	red = defcolor
	blue = defcolor
	green = defcolor
endif

bs=1
plot, oh, fsil, /nodata, $
	xtitle='log L!IOH!N [L'+sunsymbol()+']', $
	ytitle='9.7 !7l!3m silicate strength', $
	title='OHM strength vs. silicate depth', $
	xr=[0,4.5], /xstyle,$
	yr=[-5,0.4], /ystyle, $
	charsize=cs, thick=lthick, charthick=cthick

hor,-1.1,linestyle=2,thick=4

oplot, oh, fsil, color=red, psym = symcat(14),symsize=1.5
oplot, oh_arch, asil, color=red, psym = symcat(14),symsize=1.5
oplot, oh_con, csil, color=blue, psym = symcat(15),symsize=1.5
arrow, oh_con, csil, oh_con-0.2, csil, color=blue,/data,thick=3



if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif



stop
end
