pro agn_fraction_trends, ps = ps, stop = stop, label=label
;+
; NAME:
;       
;	AGN_FRACTION
;
; PURPOSE:
;
;	Plot AGN fraction for CSOs, looking at trends across sample
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
;	IDL> agn_fraction_trends, /ps
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

!p.multi = [0,1,3]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/agn_fraction_trends.ps',/color,/portrait, xoff=1,yoff=1, xs=18, ys=24
	defcolor = fsc_color("Black")
	cs = 1
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

xyouts, 2.5d-3, 3d-2, 'AGN fraction', /data, color=red, charthick = cthick
xyouts, 3d-3, 7d-3, 'Starburst fraction', /data, color=green, charthick = cthick

xyouts, x[redline[4]] * 0.65, y[redline[4]],'100%', color=red, charthick = cthick
xyouts, x[redline[3]] * 0.65, y[redline[3]], '50%', color=red, charthick = cthick
xyouts, x[redline[2]] * 0.65, y[redline[2]], '25%', color=red, charthick = cthick
xyouts, x[redline[1]] * 0.65, y[redline[1]], '10%', color=red, charthick = cthick

xyouts, x[greenline[4]], y[greenline[4]] * 0.7,'100%', color=green, charthick = cthick
xyouts, x[greenline[3]], y[greenline[3]] * 0.7, '50%', color=green, charthick = cthick
xyouts, x[greenline[2]], y[greenline[2]] * 0.7, '25%', color=green, charthick = cthick
xyouts, x[greenline[1]], y[greenline[1]] * 0.7, '10%', color=green, charthick = cthick

; CSO data

oplot, pah62, oIV / neII, psym=symcat(3), color=defcolor, symsize=1.5

; 4C37.11 - limit on PAH flux

pah62limarr = csodat('pah62ew_lim')

oplot, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], psym=symcat(3), color=defcolor, symsize=1.5
arrow, pah62limarr[7], oIV_flux[7]/neII_flux[7], pah62limarr[7] * 0.7, oIV_flux[7]/neII_flux[7], /data, hsize=hsize, color=defcolor

; 1946+70 - limit on OIV flux

oIVlim_cso007 = linelim('cso007','oIV',/lh,/noplot) * 1d-21

oplot, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], psym=symcat(3), color=defcolor, symsize=1.5
arrow, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], $
	[pah62cso[6]], [oIVlim_cso007/neII_flux[6]] * 0.6, /data, hsize=hsize, color=defcolor

; Kinematic jet age

	trend = alog10([34,550,92,1700,1d5])
	trend7 = alog10([502])
	trend6 = alog10([4000])

	for i=0, n_elements(trend) - 1 do $
		oplot, [pah62[i]], [oIV[i] / neII[i]], psym=symcat(16), color=defcolor, symsize=trend[i]
	oplot, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], psym=symcat(16), color=defcolor, symsize=trend7
	oplot, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], psym=symcat(16), color=defcolor, symsize=trend6

	xyouts, /data, 0.3, 5.0, 'Kinematic jet age', charsize=cs, color=defcolor, charthick=cthick

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

xyouts, 2.5d-3, 3d-2, 'AGN fraction', /data, color=red, charthick = cthick
xyouts, 3d-3, 7d-3, 'Starburst fraction', /data, color=green, charthick = cthick

xyouts, x[redline[4]] * 0.65, y[redline[4]],'100%', color=red, charthick = cthick
xyouts, x[redline[3]] * 0.65, y[redline[3]], '50%', color=red, charthick = cthick
xyouts, x[redline[2]] * 0.65, y[redline[2]], '25%', color=red, charthick = cthick
xyouts, x[redline[1]] * 0.65, y[redline[1]], '10%', color=red, charthick = cthick

xyouts, x[greenline[4]], y[greenline[4]] * 0.7,'100%', color=green, charthick = cthick
xyouts, x[greenline[3]], y[greenline[3]] * 0.7, '50%', color=green, charthick = cthick
xyouts, x[greenline[2]], y[greenline[2]] * 0.7, '25%', color=green, charthick = cthick
xyouts, x[greenline[1]], y[greenline[1]] * 0.7, '10%', color=green, charthick = cthick

; CSO data

oplot, pah62, oIV / neII, psym=symcat(3), color=defcolor, symsize=1.5

; 4C37.11 - limit on PAH flux

pah62limarr = csodat('pah62ew_lim')

oplot, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], psym=symcat(3), color=defcolor, symsize=1.5
arrow, pah62limarr[7], oIV_flux[7]/neII_flux[7], pah62limarr[7] * 0.7, oIV_flux[7]/neII_flux[7], /data, hsize=hsize, color=defcolor

; 1946+70 - limit on OIV flux

oIVlim_cso007 = linelim('cso007','oIV',/lh,/noplot) * 1d-21

oplot, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], psym=symcat(3), color=defcolor, symsize=1.5
arrow, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], $
	[pah62cso[6]], [oIVlim_cso007/neII_flux[6]] * 0.6, /data, hsize=hsize, color=defcolor

; Silicate depth

	trend = [-0.7,-0.8,0.3,-0.5,-0.4] * (-1) + 0.5
	trend7 = [-0.5] * (-1) + 0.5
	trend6 = [-0.6] * (-1) + 0.5

	for i=0, n_elements(trend) - 1 do $
		oplot, [pah62[i]], [oIV[i] / neII[i]], psym=symcat(16), color=defcolor, symsize=trend[i]
	oplot, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], psym=symcat(16), color=defcolor, symsize=trend7
	oplot, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], psym=symcat(16), color=defcolor, symsize=trend6

	xyouts, /data, 0.3, 5.0, 'Silicate depth', charsize=cs, color=defcolor, charthick=cthick

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

xyouts, 2.5d-3, 3d-2, 'AGN fraction', /data, color=red, charthick = cthick
xyouts, 3d-3, 7d-3, 'Starburst fraction', /data, color=green, charthick = cthick

xyouts, x[redline[4]] * 0.65, y[redline[4]],'100%', color=red, charthick = cthick
xyouts, x[redline[3]] * 0.65, y[redline[3]], '50%', color=red, charthick = cthick
xyouts, x[redline[2]] * 0.65, y[redline[2]], '25%', color=red, charthick = cthick
xyouts, x[redline[1]] * 0.65, y[redline[1]], '10%', color=red, charthick = cthick

xyouts, x[greenline[4]], y[greenline[4]] * 0.7,'100%', color=green, charthick = cthick
xyouts, x[greenline[3]], y[greenline[3]] * 0.7, '50%', color=green, charthick = cthick
xyouts, x[greenline[2]], y[greenline[2]] * 0.7, '25%', color=green, charthick = cthick
xyouts, x[greenline[1]], y[greenline[1]] * 0.7, '10%', color=green, charthick = cthick

; CSO data

oplot, pah62, oIV / neII, psym=symcat(3), color=defcolor, symsize=1.5

; 4C37.11 - limit on PAH flux

pah62limarr = csodat('pah62ew_lim')

oplot, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], psym=symcat(3), color=defcolor, symsize=1.5
arrow, pah62limarr[7], oIV_flux[7]/neII_flux[7], pah62limarr[7] * 0.7, oIV_flux[7]/neII_flux[7], /data, hsize=hsize, color=defcolor

; 1946+70 - limit on OIV flux

oIVlim_cso007 = linelim('cso007','oIV',/lh,/noplot) * 1d-21

oplot, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], psym=symcat(3), color=defcolor, symsize=1.5
arrow, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], $
	[pah62cso[6]], [oIVlim_cso007/neII_flux[6]] * 0.6, /data, hsize=hsize, color=defcolor

; Optical type

	trend = ['LINER','WLRG','BLRG','NLRG','LINER']
	trend7 = ['NLRG']
	trend6 = ['?']

	for i=0, n_elements(trend) - 1 do $
		xyouts, [pah62[i]], [oIV[i] / neII[i]], trend[i], /data, charsize=cs, color=defcolor
	xyouts, [pah62limarr[7]], [oIV_flux[7]/neII_flux[7]], trend7, /data, charsize=cs, color=defcolor
	xyouts, [pah62cso[6]], [oIVlim_cso007/neII_flux[6]], trend6, /data, charsize=cs, color=defcolor

	xyouts, /data, 0.3, 5.0, 'Optical type', charsize=cs, color=defcolor, charthick=cthick

;legend,['AGN','ULIRG','Starburst','CSO'],psym=[6,4,9,46], $
;	color=[red, blue, green, defcolor], $ 
;	symthick = 4, $
;	thick = lthick, charthick = cthick, /top, /right

;if keyword_set(label) then begin
;	xyouts, pah62, oIV / neII, csoobj, charsize=1, color=defcolor, charthick = cthick
;	xyouts, [pah62limarr[7]] * 1.1, [oIV_flux[7]/neII_flux[7]], '4C37.11', color=defcolor, charsize=1, charthick = cthick
;	xyouts, [pah62cso[6]] * 1.1, [oIVlim_cso007/neII_flux[6]], '1946+70', color=defcolor, charsize=1, charthick = cthick
;endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
