pro ohsil, ps = ps, bw = bw, label=label
;+
; NAME: 
;       OHSIL 
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
;	Added completed data set - KW, Mar 08
;-

oname=ohmdat('tag')
o = ohmdat('sil')
osil = float(o(0,*))
badosil = where(oname eq 'mega034')
goodosil = setdifference(indgen(n_elements(osil)),badosil)
osil = osil(goodosil)

cname=condat('tag')
c = condat('sil')
csil = float(c(0,*))
badcsil = where(oname eq 'control033')
goodcsil = setdifference(indgen(n_elements(csil)),badcsil)
csil = csil(goodcsil)

aname=archdat('tag')
a = archdat('sil')
asil = float(a(0,*))

oh = float(ohmdat('logoh'))
oh_con = float(condat('logoh'))
oh_arch = float(archdat('logoh'))

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = '~/Astronomy/Research/Spitzer/plots/ohsil.ps', /color
	cs = 1.5
	cthick = 2
	lthick = 2
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

plot, oh, osil, /nodata, $
	xtitle='log L!IOH!N!E(max)!N [L'+sunsymbol()+']', $
	ytitle='9.7 !7l!3m silicate strength', $
	title='OHM strength vs. silicate depth', $
	xr=[0,4.5], /xstyle,$
	yr=[-5,0.4], /ystyle, $
	charsize=cs, $
	thick=lthick, $
	charthick=cthick

hor,-1.1,linestyle=2,thick=2

oplot, oh, osil, color=red, psym = symcat(14),symsize=1.5
oplot, oh_arch, asil, color=red, psym = symcat(16),symsize=1.5
oplot, oh_con, csil, color=blue, psym = symcat(15),symsize=1.5
arrow, oh_con, csil, oh_con-0.2, csil, color=blue,/data,thick=2

legend, /bottom,/right, ['OHMs (Darling)','OHMs (archived)','Control sample'],psym=[14,16,15],$
	color=[red,red,blue], charsize=cs, thick=lthick, charthick=cthick

if keyword_set(label) then begin
	xyouts, oh, osil-0.1, fname, color=red
	xyouts, oh_arch, asil-0.1, aname, color=red
	xyouts, oh_con, csil-0.1, cname, color=blue
endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif



end
