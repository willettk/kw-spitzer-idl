;+
; NAME:
;       
;	CSO_POSTER
;
; PURPOSE:
;
;	Create plots for poster presentation at Spitzer legacy conference
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
;	IDL> 
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Oct 09
;-

!p.font=0

tags = csodat('tag',/rasort)
tags = tags[1,*]
tags = tags[where(tags ne 'cso010')]
ncso = n_elements(tags)

set_plot,'ps'
device,filename = '~/Astronomy/Research/Spitzer/cso/conference/lrspectra.ps', $
	/color, /portrait, xs=24, ys=24, xoff=1, yoff=1, /encap
cs = 1
ls = 2
	
erase 

multiplot, [3,3], $
	/square, $
	;mxtitle = 'Wavelength [um]', $
	;mytitle = 'Flux density [Jy]', $
	mxtitsize = 2, $
	mytitsize = 2, $
	mtitsize = 1.5, $
	mcharthick = 4, $
	mxcharthick = 4, $
	mycharthick = 4

yr = [$
	[5d-4,5d-2], $		; cso002
	[3d-4,1d-1], $		; cso008
	[4d-3,1d-1], $		; cso001
	[8d-3,5d-0], $		; cso006
	[3d-2,1d-0], $		; cso004
	[7d-3,3d-1], $		; cso005
	[1d-2,5d-0], $		; cso003
	[3d-3,5d-1], $ 		; cso009
	[3d-4,4d-2]]		; cso007

x1 = 0.11
x2 = 0.4
x3 = 0.7

y1 = 0.9
y2 = 0.6
y3 = 0.31

xobj = [x1,x2,x3,x1,x2,x3,x1,x2,x3]
yobj = [y1,y1,y1,y2,y2,y2,y3,y3,y3]

for i = 0, ncso - 1 do begin

	fname = tags[i]

	tag, fname, dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
	
	flux = sed.flux_lr
	wave = sed.wave_lr
	order = sed.order_lr
	
	noneg = where(sed.flux_lr gt 0)
	flux = flux[noneg]
	wave = wave[noneg]
	order = order[noneg]
	
	; Identified IR lines for overplotting
	
	templines = ir_lines(/lr)
	lines = templines(*,0) & line_id = templines(*,1)
	
	; Plot data
	
	defcolor = fsc_color('Black')
	lthick = 7
	cthick = 7
	cs = 1.5

	xr = [fix(min(wave)),fix(max(wave))+3]
	;yr = [1d-3,max(flux)]
	
	sl_bonus = where(order eq 3 and wave lt 10.)
	ll_bonus = where(order eq 3 and wave gt 10.)
	
	sl2_index = where((order eq 2) and (wave lt 10.))
	sl1_index = where((order eq 1) and (wave lt 15.))
	ll2_index = where((order eq 2) and (wave gt 10.))
	ll1_index = where((order eq 1) and (wave gt 15.))
	
	plot, wave, flux, $
		/xlog, /ylog, $
	;	xticks = 13, $
	;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
		xrange = xr, /xstyle, $
		yrange = yr[*,i], /ystyle, $
;		title = sed.obj, $
		charsize = cs, $
		color = defcolor, $
		thick = lthick, $
		xthick = lthick, $
		ythick = lthick, $
		charthick = cthick, $
		/nodata

	xyouts, xobj[i], yobj[i], sed.obj, /normal, charthick=4
	
	oplot, wave,flux, psym = 10, thick = 4

		case i of
			0: begin
				pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
			   	ylabel = 0.02
			   end
			1: begin
				pahlines = [            7.7,    11.3,12.7                    ]
			   	ylabel = 0.03
			   end
			2: begin
				pahlines =             [7.7,8.6,11.3,12.7,     16.4,17.1]
			   	ylabel = 0.05
			   end
;			3: pahlines = [0]
			3: begin
				pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2]
			   	ylabel = 1.1
			   end
			4: begin
				pahlines =         [6.2,7.7,8.6,11.3,12.7,14.2,16.4]
			   	ylabel = 0.56
			   end
			5: begin
				pahlines = [0]
			   	ylabel = 0.1
			   end
			6: begin
				pahlines = [5.3,5.7,6.2,7.7,8.6,11.3,12.7,14.2,16.4,     17.4]
			   	ylabel = 1
			   end
			7: begin
				pahlines = [    5.7,6.2,7.7,8.6,11.3,12.7,     16.4,     17.4]
			   	ylabel = 0.12
			   end
			8: begin
				pahlines = [            7.7,8.6,11.3,12.7,14.2               ]
			   	ylabel = 0.015
			   end
		endcase

	; Fine-structure emission

	labelsize = 0.7

	if min(where('cso00'+strtrim([1,2,3,4,6,8,9],2) eq fname)) ne -1 then begin
		lw = 25.9 & line_id = '[OIV]' & st = 7d-2
		xyouts, 25.9, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif

	if min(where('cso00'+strtrim([1,2,3,4,6,7,8,9],2) eq fname)) ne -1 then begin
		lw = 15.6 & line_id = '[NeIII]' & st = 3d-2
		xyouts, 15.8, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif

	if min(where('cso00'+strtrim([1,2,3,4,6,7,8,9],2) eq fname)) ne -1 then begin
		lw = 12.8 & line_id = '[NeII]' & st = 6d-2
		xyouts, 13.0, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif
		
	if min(where('cso00'+strtrim([1,2,3,4,6,8,9],2) eq fname)) ne -1 then begin
		lw = 18.7 & line_id = '[SIII]' & st = 3d-2
		xyouts, lw, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif

	if min(where('cso00'+strtrim([3,6],2) eq fname)) ne -1 then begin
		lw = 10.5 & line_id = '[SIV]' & st = 9d-3
		xyouts, lw, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif

	if min(where('cso00'+strtrim([1,2,3,4,6,7,8,9],2) eq fname)) ne -1 then begin
		lw = 17.0 & line_id = 'H!I2!N' & st = 3.5d-2
		xyouts, lw, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif

	if min(where('cso00'+strtrim([1,2,3,7,6,8,9],2) eq fname)) ne -1 then begin
		lw = 9.7 & line_id = 'H!I2!N' & st = 2d-2
		xyouts, lw, ylabel, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	endif


		if fname ne 'cso005' then begin

			bheight = 10d^(!y.crange[0]+0.13*(!y.crange[1] - !y.crange[0]))
			blevel  = 10d^(!y.crange[0]+0.10*(!y.crange[1] - !y.crange[0]))
			blabel  = 10d^(!y.crange[0]+0.04*(!y.crange[1] - !y.crange[0]))
		
			plots, [pahlines(0),pahlines[n_elements(pahlines)-1]], [blevel,blevel], thick = lthick
			for k=0,n_elements(pahlines)-1 do plots, [pahlines(k),pahlines(k)], [blevel,bheight], thick = lthick
			xyouts,sqrt(pahlines(0)*pahlines[n_elements(pahlines)-1]),blabel,'PAH',charsize=0.9, charthick = lthick

		endif
	
		; Label the silicate features
	
		sillevel = 10d^(!y.crange[0]+0.7*(!y.crange(1) - !y.crange(0)))
		sillevel2 = 10d^(!y.crange[0]+0.75*(!y.crange(1) - !y.crange(0)))
		sillabel = 10d^(!y.crange[0]+0.85*(!y.crange(1) - !y.crange(0)))
;		plots, [8,12], [sillevel2,sillevel2], linestyle=1, thick=lthick
;		plots, [17,20], [sillevel2,sillevel2], linestyle=1, thick=lthick
;		xyouts, 12, sillabel, 'Silicate', charsize = 0.7, charthick = lthick
	

	multiplot

endfor

multiplot,/reset

device, /close
set_plot,'x'

;!p.font=-1

; Feature-feature Sirocky diagram

csojunk = csodat('sil')
csoerrjunk = csodat('silerr')
csoobjjunk = csodat('obj')

sil10_cso = csojunk[0,*]
sil18_cso = csojunk[1,*]
sil10_cso_err = abs(csoerrjunk[0,*])
sil18_cso_err = abs(csoerrjunk[1,*])

; Remove VII Zw 485

sil10_cso = sil10_cso[0:8]
sil18_cso = sil18_cso[0:8]
sil10_cso_err = sil10_cso_err[0:8]
sil18_cso_err = sil18_cso_err[0:8]
csoobj=csoobjjunk[0:8]

set_plot,'ps'
device,filename='~/Astronomy/Research/Spitzer/cso/conference/feature2.ps',$
	/color,/portrait, xs=24, ys=24, yoff=1, xoff=1, /encap
defcolor = fsc_color("Black")
cs = 1.8
lthick = 10
cthick = 5

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

; CSO data

oploterror, sil10_cso, sil18_cso, sil10_cso_err, sil18_cso_err, $
	psym=symcat(16), color=defcolor, errthick = 4, errcolor=defcolor, /nohat

; Starting point

oplot,[xstart],[ystart], psym=symcat(46), symsize=3

legend,['Y=100', 'Y=200','Y=400','p=0','p=1','p=2','slab','clumps'],linestyle=[0,1,2,0,0,0,0,0], $
	color=[defcolor,defcolor,defcolor,defcolor,red,blue,yellow,green], $ 
	thick = lthick, charthick = cthick, /top, /left, $
	charsize=2

device,/close
set_plot,'x'


; Fork (Spoon) diagram 

restore,'~/Astronomy/Research/Spitzer/icelist62.sav'

; PAH data (spline fit should be used, since Henrik's data was measured using the same method)

sil_cso = csodat('sil') & sil_cso = sil_cso[0,*]
pah62ew_cso = csodat('pah62ew') & pah62ew_cso = pah62ew_cso[0,*]
pah62ew_csoice = csodat('pah62ew_ice') & pah62ew_csoice = pah62ew_csoice[0,*]
pah62ewlim_cso = csodat('pah62ew_lim') & pah62ewlim_cso = pah62ewlim_cso[0,*]
tag_cso = transpose(csodat('tag'))

match, tag_cso, icelist62, cind, oiceind, count = ocount
if ocount gt 0 then pah62ew_cso[cind] = pah62ew_csoice[cind]

; Remove spurious data for VII Zw 485 (cso010)

badcso = where(tag_cso eq 'cso010')
goodcso = setdifference(indgen(n_elements(tag_cso)),badcso)
sil_cso = sil_cso[goodcso]
pah62ew_cso = pah62ew_cso[goodcso]
pah62ewlim_cso = pah62ewlim_cso[goodcso]
tag_cso = tag_cso[goodcso]

; MEGA data

sil_ohm = ohmdat('sil') & sil_ohm = sil_ohm[0,*]
pah62ew_ohm = ohmdat('pah62ew') & pah62ew_ohm = pah62ew_ohm[0,*]
pah62ewlim_ohm = ohmdat('pah62ew_lim') & pah62ewlim_ohm = pah62ewlim_ohm[0,*]
pah62ew_ohmice = ohmdat('pah62ew_ice') & pah62ew_ohmice = pah62ew_ohmice[0,*]
tag_ohm = transpose(ohmdat('tag'))

match, tag_ohm, icelist62, oind, oiceind, count = ocount
if ocount gt 0 then pah62ew_ohm[oind] = pah62ew_ohmice[oind]

; Remove mega034

badohm = where(tag_ohm eq 'mega034')
goodohm = setdifference(indgen(n_elements(tag_ohm)),badohm)
sil_ohm = sil_ohm[goodohm]
pah62ew_ohm = pah62ew_ohm[goodohm]
pah62ewlim_ohm = pah62ewlim_ohm[goodohm]
tag_ohm = tag_ohm[goodohm]

; PAHFIT data

pahfit62ew_cso = csodat('pahfit62ew')     & pahfit62ew_cso = pahfit62ew_cso[0,*]
pahfit62ew_csoice = csodat('pahfit62ew_ice') & pahfit62ew_csoice = pahfit62ew_csoice[0,*]

if ocount gt 0 then pahfit62ew_cso[oind] = pahfit62ew_csoice[oind]
pahfit62ew_cso = pahfit62ew_cso[goodcso]

; Read in Henrik's data

restore, '~/Astronomy/Research/Spitzer/henrik/fork.diagram.parameters.xdr'

; Remove overlapping objects in both samples

obj_cso = csodat('obj') & obj_cso = obj_cso[goodcso]

allobj = [obj_cso]
henobj = 'IRAS '+target
henobj[where(henobj eq 'IRAS NGC7469')] = 'IRAS 23007+0836'
henobj[where(henobj eq 'IRAS Mrk1014')] = 'IRAS 01572+0009'
henobj[where(henobj eq 'IRAS Arp220')]  = 'IRAS 15327+2340'
henobj[where(henobj eq 'IRAS Mrk273')]  = 'IRAS 13428+5608'
henobj[where(henobj eq 'IRAS Mrk231')]  = 'IRAS 12540+5708'
henobj[where(henobj eq 'IRAS UGC5101')] = 'IRAS 09320+6134'

match, allobj, henobj, myind, henind, count = count
goodhen = setdifference(indgen(n_elements(target)),henind)
target = target[goodhen]
galtype = galtype[goodhen]
iceew62 = iceew62[goodhen]
silstrength = silstrength[goodhen]

; Objects on which to plot limits

limlist = ['cso005','cso006','cso007','cso008']
;limlist = ['cso005','cso006','cso007','cso008','cso010']

match, tag_cso, limlist, csolim, csolim2
nolimits = setdifference(indgen(n_elements(tag_cso)),csolim)

olimlist = ['mega007']

match, tag_ohm, olimlist, ilim, ilim2
onolimits = setdifference(indgen(n_elements(tag_ohm)),ilim)

; Sort Henrik's data by galaxy type

wherea = where(galtype eq 'A')				; AGN
whereo = where(galtype eq 'O')				; Obscured
whereu = where(galtype eq 'U' or galtype eq 'H')	; ULIRGs and HyLIRGs
wheres = where(galtype eq 'S')				; Starbursts
wheren = where(galtype eq 'N')				; Normal

pahew_a = iceew62[wherea] & sil_a = silstrength[wherea]
pahew_o = iceew62[whereo] & sil_o = silstrength[whereo]
pahew_u = iceew62[whereu] & sil_u = silstrength[whereu]
pahew_s = iceew62[wheres] & sil_s = silstrength[wheres]
pahew_n = iceew62[wheren] & sil_n = silstrength[wheren]

	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/conference/fork.ps', $
		/color, /portrait, xs=24, ys=24, xoff=1, yoff=1, /encap
	cs = 1.5
	ls = 4
	defcolor = fsc_color("Black")

red = fsc_color("Red")
blue = fsc_color("Blue")
lightblue = fsc_color("Cyan")
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")
green = fsc_color("Dark Green")
purple = fsc_color("Purple")
color6 = fsc_color("Rosy Brown")

plot, pah62ew_cso, sil_cso, $
	/nodata, $
	xtitle = '6.2 um PAH EW [um]', $
	ytitle = 'Silicate strength', $
	charsize = cs, $
	thick = ls, $
	xthick = 6, $
	ythick = 6, $
	charthick = ls, $
	/xlog, /xstyle, $
	ystyle = 1, $
	xr = [1d-3, 3], $
	yr = [1, -5], $
	color = defcolor

; Insert Henrik's divisions

;ver, 0.07
;ver, 0.45
;hor, -0.8
;hor, -2.4
hor, 0, linestyle = 1, thick = ls

; Henrik's data

oplot, pahew_a, sil_a, $
	psym = symcat(18), thick = ls, color = orange

oplot, pahew_o, sil_o, $
	psym = symcat(7), thick = ls, color = green

oplot, pahew_s, sil_s, $
	psym = symcat(46), thick = ls, color = purple

oplot, pahew_u, sil_u, $
	psym = symcat(2), thick = ls, color = defcolor

oplot, pahew_n, sil_n, $
	psym = symcat(15), thick = ls, color = blue

; OHM data

oplot, pah62ew_ohm[onolimits], sil_ohm[onolimits], $
	psym = symcat(2), thick = ls, color = defcolor

; My data

oplot, pah62ew_cso[nolimits], sil_cso[nolimits], $
	psym = symcat(16), thick = ls, color = red

oplot, pah62ewlim_cso[csolim], sil_cso[csolim], $
	psym = symcat(16), thick = ls, color = red

arrow, pah62ewlim_cso[csolim], sil_cso[csolim], $
	pah62ewlim_cso[csolim] * 0.7, sil_cso[csolim], $
	hsize = 290, $
	thick = 5, color = red, /data

; Difference between spline and PAHFIT results

mean62 = mean(double(pah62ew_cso[nolimits]))
diffvec = abs((mean(double(pah62ew_cso[nolimits])) - mean(double(pahfit62ew_cso[nolimits]))))

arrow, mean62, 0.5, mean62 + diffvec, 0.5, color=defcolor, thick=4, /data

if keyword_set(label) then xyouts, pah62ew_cso, sil_cso, obj_cso, color=red, charthick = ls, charsize = 1.0, /data 

legend, /top, /right, $
	['CSOs', 'AGN', 'Obscured AGN', 'Normal galaxies', 'Starbursts', 'ULIRG/HyLIRGs'], $
	psym = [16,18,7,15,46,2], $
	color = [red, orange, green, blue, purple, defcolor], $
	thick = ls, charthick = ls, charsize = 1.5

device,/close
set_plot,'x'

; Neon ratio (ionization levels)

; Read in data from Farrah Tbl. 3

readcol,'~/Astronomy/Research/Spitzer/OHM/farrah_tbl3.txt',$
	gal, arIII, h2s3, sIV, h2s2, neII, neV14, neIII, h2s1, sIII, neV24, oIV, h2s0, sIII33, siII34, $
	format = 'a,a,a,a,a,a,a,a,a,a,a,a,a,a,a',/silent

all = [[gal], [arIII], [h2s3], [sIV], [h2s2], [neII], [neV14], [neIII], [h2s1], [sIII], $
	[neV24], [oIV], [h2s0], [sIII33], [siII34]]

sz = size(all)
lim=intarr(sz(1),sz(2))

for j = 0, sz(2)-1 do begin
	for i = 0, sz(1)-1 do begin
		cell = all(i,j)
		if strmid(cell,0,1) eq '<'  then begin
			cell = strmid(cell,1,strlen(cell)-1)
			lim(i,j) = 1
		endif
		if strmid(cell,0,2) eq '\l' then cell = '0.'
		if strmid(cell,strlen(cell)-2,2) eq '::' then cell = strmid(cell,0,strlen(cell)-2)
		if strmid(cell,strlen(cell)-1,1) eq ':'  then cell = strmid(cell,0,strlen(cell)-1)
		all(i,j) = cell
	endfor
endfor


sIII = float(all(*,9))
sIV  = float(all(*,3))
neII = float(all(*,5))
neIII= float(all(*,7))
sulfur = where(sIII ne 0. and sIV ne 0.)
neon   = where(neII ne 0. and neIII ne 0.)

s3lim = where(lim(*,9) eq 1)
s4lim = where(lim(*,3) eq 1)

slim = setunion(s3lim,s4lim)
notlimits = setdifference(indgen(53),slim)

; Load in data for our OHMs 

	restore,'~/Astronomy/Research/Spitzer/linedata/neII.sav'
	neII_ohm = line
	restore,'~/Astronomy/Research/Spitzer/linedata/neIII.sav'
	neIII_ohm = line

	; Isolate targets where both neon lines are detected
	
	match, neii_ohm.tag, neiii_ohm.tag, a, b
	
	neII_hr = neii_ohm.flux[a]
	neIII_hr = neiii_ohm.flux[b]
	
	neII_tag = neii_ohm.tag[a]
	neIII_tag = neiii_ohm.tag[b]

	; Darling OHMs
	
	neII_hr_ohm = neII_hr[where(strmid(neII_tag,0,4) eq 'mega')]
	neIII_hr_ohm = neIII_hr[where(strmid(neIII_tag,0,4) eq 'mega')]
	
	neII_tag_ohm = neII_tag[where(strmid(neII_tag,0,4) eq 'mega')]
	neIII_tag_ohm = neIII_tag[where(strmid(neIII_tag,0,4) eq 'mega')]
	
	; Archived OHMs
	
	neII_hr_arch = neII_hr[where(strmid(neII_tag,0,4) eq 'arch')]
	neIII_hr_arch = neIII_hr[where(strmid(neIII_tag,0,4) eq 'arch')]
	
	neII_tag_arch = neII_tag[where(strmid(neII_tag,0,4) eq 'arch')]
	neIII_tag_arch = neIII_tag[where(strmid(neIII_tag,0,4) eq 'arch')]
	
	; Control sample
	
	neII_hr_control = neII_hr[where(strmid(neII_tag,0,7) eq 'control')]
	neIII_hr_control = neIII_hr[where(strmid(neIII_tag,0,7) eq 'control')]
	
	neII_tag_control = neII_tag[where(strmid(neII_tag,0,7) eq 'control')]
	neIII_tag_control = neIII_tag[where(strmid(neIII_tag,0,7) eq 'control')]
	
; Add data for AGN and starburst galaxies

; Sturm et al 2002 (AGN, ISO)

sturm_names = ['NGC 1068','NGC 1365','NGC7469']
sturm_neII = [70.0, 40.9, 22.6] * 1d-1			; Sturm measures fluxes in 10^-20 W/cm^2
sturm_neIII = [160.0, 7.7, 2.2] * 1d-1			; doesn't really matter, since they're ratios, but should be consistent

; Verma et al 2003 (starbursts, ISO)

verma_names = ['NGC 3256','NGC 3690A','NGC 3690B/C']
verma_neII = [89.2, 31.9, 27.7] * 1d-1			; Verma measures fluxes in 10^-20 W/cm^2
verma_neIII = [14.2, 9.3, 20.0] * 1d-1

; Tommasin et al 2008 (AGN, IRS)

tommasin_names = ['IRAS 00198-7926','ESO 541-IG012','Mrk 1034 NED02',$
	'IRAS F15091-2107','IRAS F22017+0319','NGC 7496']
tommasin_neII = [6.19, 1.87, 34.99, 11.52, 5.95, 48.08]
tommasin_neIII = [14.03, 2.02, 3.61, 16.29, 14.07, 6.67]

red = fsc_color("Red")
blue = fsc_color("Blue")
purple = fsc_color("purple")
green = fsc_color("Forest Green")
orange = fsc_color("Orange")
purple = fsc_color("purple")

ne_farrah = neIII / neII

ne_ohm = neIII_hr_ohm / neII_hr_ohm

ne_arch = neIII_hr_arch / neII_hr_arch

ne_control = neIII_hr_control / neII_hr_control

ne_agn = [tommasin_neIII, sturm_neIII] / [tommasin_neII, sturm_neII]

ne_starburst = verma_neIII / verma_neII

; Plot data for CSO paper

; Load in CSO data

	restore,'~/Astronomy/Research/Spitzer/cso/linedata/neII.sav'
	neII_cso = line
	restore,'~/Astronomy/Research/Spitzer/cso/linedata/neIII.sav'
	neIII_cso = line

	; Isolate targets where both neon lines are detected
	
	match, neii_cso.tag, neiii_cso.tag, a, b
	
	neII_hr = neii_cso.flux[a]
	neIII_hr = neiii_cso.flux[b]
	
	neII_tag = neii_cso.tag[a]
	neIII_tag = neiii_cso.tag[b]

	neII_hr_cso = neII_hr[where(strmid(neII_tag,0,3) eq 'cso')]
	neIII_hr_cso = neIII_hr[where(strmid(neIII_tag,0,3) eq 'cso')]
	
	neII_tag_cso = neII_tag[where(strmid(neII_tag,0,3) eq 'cso')]
	neIII_tag_cso = neIII_tag[where(strmid(neIII_tag,0,3) eq 'cso')]
	
; For the ancillary data, now include ALL objects with line ratios (CSOs span a large range in L_IR)

; Sturm et al 2002 (AGN, ISO)

sturm_names = ['Cen A','Circinus','NGC 1068','NGC 1365','NGC 4151','NGC 5506','NGC 7469','NGC 7582']
sturm_neII = [22.1,90.0,70.0,40.9,11.8,5.9,22.6,14.8] * 1d-1	; Sturm measures fluxes in 10^-20 W/cm^2
sturm_neIII = [14.1,33.5,160.0,7.7,20.7,5.8,2.2,6.7] * 1d-1	

; Verma et al 2003 (starbursts, ISO)

verma_names = ['II Zw 40','M 82','NGC 3256','NGC 3690A','NGC 3690B/C','NGC 4038','NGC 5236','NGC 5253','NGC 7552']
verma_neII = [1.7,714.0,89.2,31.9,27.7,7.7,133.9,7.9,68.0] * 1d-1		; Verma measures fluxes in 10^-20 W/cm^2
verma_neIII = [17.0,126.0,14.2,9.3,20.0,6.5,6.8,28.9,5.1] * 1d-1

; Tommasin et al 2008 (AGN, IRS)

tommasin_names = ['IRAS 00198-7926','ESO 012-G021','ESO 541-IG012','NGC 424','NGC 513','Mrk 1034 NED02',$
	'ESO 545-G013','NGC 931','NGC 1125','ESO 033-G002','Mrk 6','Mrk 9','NGC 3516','NGC 3660',$
	'TOLOLO 1238-364','NGC 4748','NGC 4968','Mrk 817','IRAS F15091-2107','ESO 141-G055','NGC 6890',$
	'IRAS F22017+0319','NGC 7496']
tommasin_neII = [6.19,11.95,1.87,8.70,12.76,34.99,10.05,5.47,16.37,2.13,28.00,3.23,8.07,6.51,$
	45.15,7.37,24.90,3.83,11.52,2.24,11.32,5.95,48.08]
tommasin_neIII = [14.03,6.42,2.02,18.45,4.43,3.61,10.50,15.41,15.55,9.22,49.34,1.90,17.72,1.49,27.00,$
	15.93,33.80,4.58,16.29,5.62,6.57,14.07,6.67]

ne_cso = neIII_hr_cso / neII_hr_cso

ne_agn = [tommasin_neIII, sturm_neIII] / [tommasin_neII, sturm_neII]

ne_starburst = verma_neIII / verma_neII

; Add limits for the three CSOs without SIV detections

csolim = ['cso001','cso002','cso004']
cso_neII_lim = dblarr(n_elements(csolim))
cso_neIII_lim = dblarr(n_elements(csolim))

for i = 0, n_elements(csolim)-1 do begin
	junk2 = getlineflux(csolim[i],'neIII')
	cso_neIII_lim[i] = junk2[0]
	junk3 = getlineflux(csolim[i],'neII')
	cso_neII_lim[i] = junk3[0]
endfor

ne_csolim = cso_neIII_lim / cso_neII_lim

plotname='~/Astronomy/Research/Spitzer/cso/conference/irex_cso.ps'
;!p.multi=[0,1,4]

;erase

x0 = 0.15
x1 = 0.8
x2 = 0.95

y0 = 0.15
y1 = 0.8
y2 = 0.95

!x.style = 1
!y.style = 1

	set_plot,'ps'
	device, filename = plotname, $
		/color, /portrait, xs = 18, ys = 20, xoff = 1, yoff = 1, /encap
	cthick = 4
	lthick = 4
	labelsize = 1.5
	cs = 1
	arrowcolor=fsc_color("Grey")
	defcolor=fsc_color("Black")

erase 

multiplot, [1,4], $
	mxtitle = 'log ([Ne III]/[Ne II])', $
	mytitle = 'Count', $
	mxtitsize = 1.5, $
	mytitsize = 1.5, $
	mtitsize = 1.5, $
	mcharthick = 4, $
	mxcharthick = 4, $
	mycharthick = 4

bs = 0.4

; ULIRGs - Farrah, Darling OHMs and control

plothist, alog10([ne_farrah[notlimits], ne_ohm, ne_arch, ne_control]), /halfbin, $
;	position = [0.1, 0.73, 0.95, 0.94], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,40], $
	xticks = 5, xtickv = fillarr(0.5,-2.0,1.5), xtickname = replicate(' ',6), $
	yticks = 5, yminor = 0, ytickv = fillarr(10,0,40), ytickname = ['0','10','20','30','40'], $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

multiplot

xyouts, 0.15, 0.87, 'ULIRGs', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10([ne_farrah[notlimits], ne_ohm, ne_arch, ne_control])), linestyle = 1, thick = lthick

; AGN

plothist, alog10(ne_agn), /halfbin, $
;	position = [0.1, 0.52, 0.95, 0.73], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,20], $
	xticks = 5, xtickv = fillarr(0.5,-2.0,1.5), xtickname = replicate(' ',6), $
	yticks = 5, yminor = 0, ytickv = fillarr(5,0,15), ytickname = ['0','5','10','15'], $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

multiplot

xyouts, 0.15, 0.64, 'AGN', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10(ne_agn)), linestyle = 1, thick = lthick

; Starbursts

plothist, alog10(ne_starburst), /halfbin, $
	;position = [0.1, 0.31, 0.95, 0.52], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,5], $
	xticks = 5, xtickv = fillarr(0.5,-2.0,1.5), xtickname = replicate(' ',6), $
	yticks = 6, yminor = 0, ytickv = fillarr(2,0,6), ytickname = ['0','2','4',' '], $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

multiplot

xyouts, 0.15, 0.45, 'Starbursts', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10(ne_starburst)), linestyle = 1, thick = lthick

; CSOs

plothist, alog10(ne_cso), /halfbin, $
	;position = [0.1, 0.1, 0.95, 0.31], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,6], $
	yticks = 6, yminor = 0, ytickv = fillarr(2,0,6), ytickname = ['0','2','4',' '], $
;	xtitle = 'log ([NeIII]/[NeII])', $
;	ytitle = 'Count', $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick


xyouts, 0.15, 0.24, 'CSOs', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10(ne_cso)), linestyle = 1, thick = lthick
	
multiplot
multiplot, /default

	device,/close
	set_plot,'x'

;stop

end
