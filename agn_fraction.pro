pro agn_fraction, ps = ps, stop = stop, label=label, oivstop = oivstop
;+
; NAME:
;       
;	AGN_FRACTION
;
; PURPOSE:
;
;	Plot AGN fraction for CSOs, other galaxies using diagnostics from ARmus et al. (2007)
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
;	PS - hard copy
;
;	STOP - stop at end of program
;
; EXAMPLE:
;
;	IDL> agn_fraction, /ps
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Dec 09
;-

; PAH 6.2 um flux

pahjunk = csodat('pah62ew')
pah62cso = float(pahjunk[0,*]) 
pah62cso_err = float(pahjunk[1,*]) 
pah62cso[1] = [0.12]			; Correct 4C+31.04 for water ice absorption
csoobjjunk = csodat('obj')

; OIV and NeII fluxes

; Neon

restore,'~/Astronomy/Research/Spitzer/cso/linedata/neII.sav'
neII_cso = line

neII_tag     = neII_cso.tag
neII_flux    = neII_cso.flux
neII_fluxerr = neII_cso.fluxerr

restore,'~/Astronomy/Research/Spitzer/cso/linedata/oIV.sav'
oIV_cso = line

oIV_tag     = oIV_cso.tag
oIV_flux    = oIV_cso.flux
oIV_fluxerr = oIV_cso.fluxerr

; Pad for non-detections of atomic lines

neII_tag = [neII_tag[0:3],'cso005',neII_tag[4:7],'cso010']
neII_flux = [neII_flux[0:3],0,neII_flux[4:7],0]
neII_fluxerr = [neII_fluxerr[0:3],0,neII_fluxerr[4:7],0]

oIV_tag = [oIV_tag[0:3],'cso005',oIV_tag[4],'cso007',oIV_tag[5:6],'cso010']
oIV_flux = [oIV_flux[0:3],0,oIV_flux[4],0,oIV_flux[5:6],0]
oIV_fluxerr = [oIV_fluxerr[0:3],0,oIV_fluxerr[4],0,oIV_fluxerr[5:6],0]

; Remove objects without PAH, NeII and OIV

csoind = where(pah62cso gt 0 and neII_flux gt 0 and oIV_flux gt 0)
csoind = csoind[where(csoind ne 2 and csoind ne 9 and csoind ne 7)]			; Remove NGC 5793 and VII Zw 485

pah62 = pah62cso[csoind]
pah62_err = pah62cso_err[csoind]
neII = neII_flux[csoind]
oIV = oIV_flux[csoind]
csoobj=csoobjjunk[csoind]

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/agn_PAH_OIV.ps',/color,/portrait
	defcolor = fsc_color("Black")
	cs = 1.5
	lthick = 5
	cthick = 5
	hsize = 200
endif else begin
	defcolor=fsc_color("White")
	cs = 1.5
	lthick = 1
	cthick = 1
	hsize = 10
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
yellow = fsc_color("Goldenrod")
green = fsc_color("Dark Green")

plot, indgen(10)+1, $
	/nodata, $
	/xlog, $
	/ylog, $
	xtitle='PAH 6.2 EW [!7l!3m]', $
	ytitle='[OIV]/[NeII]', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charthick = cthick, $
	charsize=cs, $
	color=defcolor, $
	xr=[2d-3,2d0], /xstyle, $
	yr=[3d-3,2.5d1], /ystyle

; Data and tracks from Dexter on Armus et al. (2007), Fig. 6, bottom graph

dexdir='/Applications/Dexter/'
readcol,dexdir+'armus2007_fig6.gif.agn',skipline=1,x,y,/silent

redline = where(x lt 4d-3)
greenline = where(y lt 8d-3)

agn_ind = [5,6,7,9,11,12,14,15,17,18,19,21]
ulirg_ind = [8,13,16,20,22,24,26,27,29,28,31,33,41]
starburst_ind = [32,33,34,35,36,38,39,40,42,43,44,45]

agnlim_left = [5,7,6,9,11,18,14,17]
agnlim_down = [12]
ulirglim_down = [13,16,19,29,28,41,33]
starburstlim_down = [39,40]

;xyouts, x*1.1, y, strtrim(indgen(n_elements(x)),2), /data
oplot, x[agn_ind], y[agn_ind], color=red, psym = symcat(6, thick=4), symsize=1.5
oplot, x[ulirg_ind], y[ulirg_ind], color=blue, psym = symcat(4,  thick=4), symsize=1.5
oplot, x[starburst_ind], y[starburst_ind], color=green, psym = symcat(9, thick=4), symsize=1.5

arrow, x[agnlim_left], y[agnlim_left], x[agnlim_left] * 0.75, y[agnlim_left], color=red, /data, hsize = hsize
arrow, x[agnlim_down], y[agnlim_down], x[agnlim_down], y[agnlim_down] * 0.65, color=red, /data, hsize = hsize
arrow, x[ulirglim_down], y[ulirglim_down], x[ulirglim_down], y[ulirglim_down] * 0.65, color=blue, /data, hsize = hsize
arrow, x[starburstlim_down], y[starburstlim_down], x[starburstlim_down], y[starburstlim_down] * 0.65, color=green, /data, hsize = hsize

oplot, replicate(mean(x[redline]),n_elements(redline)),   y[redline], color=red, psym=-symcat(3), thick=2
oplot, replicate(mean(x[redline]),n_elements(redline)-1), y[redline[1:4]], color=red, psym=-symcat(28), thick=2
oplot, x[greenline], replicate(mean(y[greenline]),n_elements(greenline)), color=green, psym=-symcat(3), thick=2
oplot, x[greenline[1:4]], replicate(mean(y[greenline]),n_elements(greenline)-1), color=green, psym=-symcat(26), thick=2

xyouts, 2.5d-3, 3d-2, 'AGN fraction', /data, color=red, charthick = cthick, charsize=1.2
xyouts, 3d-3, 7d-3, 'Starburst fraction', /data, color=green, charthick = cthick, charsize = 1.2

xyouts, x[redline[4]] * 0.65, y[redline[4]],'100%', color=red, charthick = cthick
xyouts, x[redline[3]] * 0.65, y[redline[3]], '50%', color=red, charthick = cthick
xyouts, x[redline[2]] * 0.65, y[redline[2]], '25%', color=red, charthick = cthick
xyouts, x[redline[1]] * 0.65, y[redline[1]], '10%', color=red, charthick = cthick

xyouts, x[greenline[4]], y[greenline[4]] * 0.7,'100%', color=green, charthick = cthick
xyouts, x[greenline[3]], y[greenline[3]] * 0.7, '50%', color=green, charthick = cthick
xyouts, x[greenline[2]], y[greenline[2]] * 0.7, '25%', color=green, charthick = cthick
xyouts, x[greenline[1]], y[greenline[1]] * 0.7, '10%', color=green, charthick = cthick

; CSO data

oplot, pah62, oIV / neII, psym=symcat(46), color=defcolor, symsize=1.5

; 4C37.11 - limit on PAH flux

pah62limarr = csodat('pah62ew_lim')

oplot, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], psym=symcat(46), color=defcolor, symsize=1.5
arrow, pah62limarr[7], oIV_flux[7]/neII_flux[7], pah62limarr[7] * 0.7, oIV_flux[7]/neII_flux[7], /data, hsize=hsize, color=defcolor

; 1946+70 - limit on OIV flux

oIVlim_cso007 = linelim('cso007','oIV',/lh,/noplot) * 1d-21

oplot, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], psym=symcat(46), color=defcolor, symsize=1.5
arrow, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], $
	[pah62cso[6]], [oIVlim_cso007/neII_flux[6]] * 0.6, /data, hsize=hsize, color=defcolor

legend,['AGN','ULIRG','Starburst','CSO'],psym=[6,4,9,46], $
	charsize = 1.2, $
	color=[red, blue, green, defcolor], $ 
	symthick = 4, $
	thick = lthick, charthick = cthick, /top, /right

if keyword_set(label) then begin
	xyouts, pah62, oIV / neII, csoobj, charsize=1, color=defcolor, charthick = cthick
	xyouts, [pah62limarr[7]] * 1.1, [oIV_flux[7]/neII_flux[7]], '4C37.11', color=defcolor, charsize=1, charthick = cthick
	xyouts, [pah62cso[6]] * 1.1, [oIVlim_cso007/neII_flux[6]], '1946+70', color=defcolor, charsize=1, charthick = cthick
endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(oivstop) then stop

; Figure 8 from Armus+07 (based on classification diagram in Laurent+00)

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/agn_triangle.ps',/color,/portrait
	defcolor = fsc_color("Black")
	cs = 1.5
	lthick = 5
	cthick = 5
	hsize = 100
endif else begin
	defcolor=fsc_color("White")
	cs = 1.5
	lthick = 1
	cthick = 1
	hsize = 10
endelse


readcol,dexdir+'armus2007_fig8.gif.agn',skipline=1,x2,y2,/silent

; CSO data

agnind2 = [0,1,3,4,5,6,7,8] + 1

f5  = fltarr(n_elements(agnind2))
f15 = fltarr(n_elements(agnind2))
f62 = fltarr(n_elements(agnind2))
f62cont = fltarr(n_elements(agnind2))
csofluxobj = strarr(n_elements(agnind2))

for i=0,n_elements(agnind2) - 1 do begin

	filename = 'cso'+string(agnind2[i],format='(i03)')
	restore, '~/Astronomy/Research/Spitzer/cso/data/structures/'+filename+'.sav'
	
	stdpivot_62 = [5.15, 5.55, 5.95, 6.55, 7.1]

	case agnind2[i] of
		0: pivots_62 = [5.55, 5.95, 6.55, 7.1]		; cso001
		1: pivots_62 = [5.2,7.2,7.7]			; cso002
		2: pivots_62 = stdpivot_62			; cso003
		3: pivots_62 = stdpivot_62			; cso004
		4: pivots_62 = stdpivot_62			; cso005
		5: pivots_62 = stdpivot_62			; cso006
		6: pivots_62 = stdpivot_62			; cso007
		7: pivots_62 = [5.15,5.95,6.55,7.2]		; cso008
		8: pivots_62 = stdpivot_62			; cso009
		9: pivots_62 = [5.15,5.95,6.55,7.5]		; cso010
	endcase
		
	cso_spl = [4,2,3,4,4,1,1,5,1,1]
	pahspl_nu, filename, spl=cso_spl[agnind2[i]], cont62 = cont62, pivots_62 = pivots_62, /quiet, /noplot, ps=ps

	wave = sed.wave_lr
	flux = sed.flux_lr
	
	f5[i]  = flux(closetomed(wave,flux,5.5))
	f15[i] = flux(closetomed(wave,flux,15.0))
	f62[i] = max(flux[closeto(wave,5.95):closeto(wave,6.55)])
	f62cont[i] = cont62

	csofluxobj[i] = sed.obj

endfor

plot, indgen(10)+1, $
	/nodata, $
	/xlog, $
	/ylog, $
	xtitle='PAH 6.2/f(5.5)', $
	ytitle='f(15)/f(5.5)', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charthick = cthick, $
	charsize=cs, $
	color=defcolor, $
	xr=[5d-3,2d1], /xstyle, $
	yr=[3d-1,5d1], /ystyle

line90 = [18,21,22,25,28,30,32,33,35,36,39,40]
line75 = [27,31,34,42,43,46,50,51,53,54,57]
line50 = [38,44,47,52,60,62,63,66,67,68,69]
triangle = [4,10,16,27,37,48,55,61,65,72,73,79,82,84,69,57,40,4]

agn2 = [0,1,3,5,7,8,9,12,13,15,17,19,20]
ulirg2 = [2,11,14,23,24,29,45,49,58,59,71,64,74]
starburst2 = [81,88,89,87,80,78,77,76,85,86,90,91]

agn2lim = [7,5,9,0,12,3,17]
ulirg2lim = [2,14]

plots, x2[line90],y2[line90],color=defcolor,linestyle=2, thick = lthick
plots, x2[line75],y2[line75],color=defcolor,linestyle=2, thick = lthick
plots, x2[line50],y2[line50],color=defcolor,linestyle=2, thick = lthick
plots, x2[triangle],y2[triangle],color=defcolor,linestyle=1, thick = lthick

oplot, x2[agn2], y2[agn2], color=red, psym=symcat(6, thick=4), symsize=1.5
oplot, x2[ulirg2], y2[ulirg2], color=blue, psym=symcat(4, thick=4), symsize=1.5
oplot, x2[starburst2], y2[starburst2], color=green, psym=symcat(9, thick=4), symsize=1.5

arrow, x2[agn2lim],       y2[agn2lim],      x2[agn2lim]*0.75,       y2[agn2lim],       color=red, /data, hsize = hsize
arrow, x2[ulirg2lim],     y2[ulirg2lim],    x2[ulirg2lim]*0.75,     y2[ulirg2lim],     color=blue, /data, hsize = hsize

;xyouts, x2*1.1, y2, strtrim(indgen(n_elements(x2)),2), /data

xyouts, x2[triangle[4]] * 0.7, y2[triangle[4]], '50%', color=defcolor, charthick = cthick
xyouts, x2[triangle[3]] * 0.7, y2[triangle[3]], '75%', color=defcolor, charthick = cthick
xyouts, x2[triangle[2]] * 0.7, y2[triangle[2]] * 1.1, '90%', color=defcolor, charthick = cthick

xyouts, x2[triangle[16]], y2[triangle[16]] * 0.75, '90%', color=defcolor, charthick = cthick
xyouts, x2[triangle[15]], y2[triangle[15]] * 0.75, '75%', color=defcolor, charthick = cthick
xyouts, x2[triangle[14]], y2[triangle[14]] * 0.75, '50%', color=defcolor, charthick = cthick

xyouts, 5d-1, 4d1, 'HII', /data, color=defcolor, charthick = cthick
xyouts, 2.5d0,4d-1,'PDR', /data, color=defcolor, charthick = cthick
xyouts, 1d-2,5d-1,'AGN', /data, color=defcolor, charthick = cthick

; CSO data

cso2lim = [3,5,6]
xval2 = (f62 - f62cont)/f5
yval2 = f15/f5

oplot, xval2, yval2, psym=symcat(46), color=defcolor, symsize=1.5
arrow, xval2[cso2lim], yval2[cso2lim], xval2[cso2lim] * 0.75, yval2[cso2lim], color=defcolor, /data, hsize = hsize

if keyword_set(label) then begin
	xyouts, xval2[0] * 1.05, yval2[0] * 0.8, csofluxobj[0], color=defcolor, charthick = cthick
	xyouts, xval2[1] * 1.05, yval2[1], csofluxobj[1], color=defcolor, charthick = cthick
	xyouts, xval2[2] * 1.05, yval2[2], csofluxobj[2], color=defcolor, charthick = cthick
	xyouts, xval2[3] * 0.12, yval2[3] * 1.05, csofluxobj[3], color=defcolor, charthick = cthick
	xyouts, xval2[4] * 1.05, yval2[4], csofluxobj[4], color=defcolor, charthick = cthick
	xyouts, xval2[5] * 1.05, yval2[5] * 0.85, csofluxobj[5], color=defcolor, charthick = cthick
	xyouts, xval2[6] * 1.05, yval2[6], csofluxobj[6], color=defcolor, charthick = cthick
	xyouts, xval2[7] * 1.05, yval2[7], csofluxobj[7], color=defcolor, charthick = cthick
endif

legend,['AGN','ULIRG','Starburst','CSO'],psym=[6,4,9,46], $
	charsize = 1.2, $
	color=[red, blue, green, defcolor], $ 
	symthick = 4, $
	thick = lthick, charthick = cthick, /top, /right


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
