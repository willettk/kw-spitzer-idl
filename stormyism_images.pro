
;+
; NAME:
;       
;	STORMYISM_IMAGES.pro
;
; PURPOSE:
;
;	Make poster-quality images for Stormy ISM conference in Pasadena
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
;       Written by K. Willett                Sep 10
;-

; SPECAVG

;######################
;#######  LR  #########
;######################

; Load names of files from data structures

onames = ohmdat('tag')
badohm = where(onames eq 'mega034')
goodohm = setdifference(indgen(n_elements(onames)),badohm)
onames = onames(goodohm)
no = n_elements(onames)

anames = transpose(archdat('tag'))
na = n_elements(anames)

cnames = condat('tag')
badcon = where(cnames eq 'control033')
goodcon = setdifference(indgen(n_elements(cnames)),badcon)
cnames = cnames(goodcon)
nc = n_elements(cnames)

; Wavelength grid over which to plot data

wavesep = 0.087
grid = fillarr(wavesep,5.0,30.0)		; Avg. separation between wave bins is 0.087 um

allspec = fltarr(n_elements(grid),no+na)
conspec = fltarr(n_elements(grid),nc)

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Dark Grey")
white = fsc_color("White")
black = fsc_color("Black")

; Loop over OHMs
	 
allohms = [onames,anames]

for i = 0, n_elements(allohms) - 1 do begin

	; Restore data

	tag, allohms[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,13.2))
	newy = newy * 0.1 / newy15

	allspec[*,i] = newy
endfor

; Plot the median template from all spectra

meanohm = median(allspec,dim=2)
stdohm = meanohm
for i=0,n_elements(stdohm)-1 do stdohm[i] = stddev(allspec[i,*])

; Control sample

for i = 0, nc - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/control/data/structures/'+cnames[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,13.2))

	newy = newy * 0.1 / newy15

	conspec[*,i] = newy
endfor

meancon = median(conspec,dim=2)
stdcon = meancon
for i=0, n_elements(stdcon)-1 do stdcon[i] = stddev(conspec[i,*])

!p.multi=[0,1,1]

	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/stormyism/specavg.ps', /color, /portrait, xs=18, ys=12
	cs = 1.8
	lthick = 7
	cthick = 7
	medthick = 5
	defcolor=black

plot,indgen(35), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [-2.5,0.5], /ystyle, $
	color=defcolor, $
	/xlog, $
	;/ylog, $
	xticks = 5, $
	xtickv = [5,10,15,20,25,30], $
	yticks = 3, $
	ytickv = [-2.0,-1.0,0.0], $
	;xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	xtitle = '!7k!3!Irest!N [!7l!3m]', $
	ytitle = 'log S!I!7m!3!N (normalized)', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, charthick = cthick, charsize = cs

ohmcolor=fsc_color("Red")
concolor=fsc_color("Blue")

;oploterror,newx,meanohm,stdohm,color=defcolor,thick=lthick,   psym=10, errthick=0.5, errcolor=grey, /nohat
;oploterror,newx+wavesep/2,meancon,stdcon,color=concolor,thick=medthick, psym=10, errthick=0.5, errcolor=fsc_color("Salmon"), /nohat

oplot,newx,alog10(meanohm),color=ohmcolor,thick=lthick, psym=10 
oplot,newx+wavesep/2,alog10(meancon),color=concolor,thick=lthick, psym=10 

;legend, /top, /left, ['OHMs ('+strtrim(n_elements(allohms),2)+')', 'non-masing ('+strtrim(nc,2)+')'], $
;	color=[defcolor,concolor], $
;	linestyle=[0,0], $
;	thick = [lthick,medthick], $
;	charsize = 1, charthick=cthick

xyouts, 15, -1.5, 'OHMs', charsize=1.8, /data, color=ohmcolor, charthick=5
xyouts, 6, -0.5, 'Non-masing', charsize=1.8, /data, color=concolor, charthick=5

ver, 13.2, linestyle=2, thick=5

	device,/close
	set_plot,'x'

; SIL_A3020

asil = archdat('sil')
osil = ohmdat('sil')
csil = condat('sil')

ohmsil = float([transpose(asil[0,*]),transpose(osil[0,*])])
consil = float([transpose(csil[0,*])])

aspindex = archdat('spindex')
ospindex = ohmdat('spindex')
cspindex = condat('spindex')

ohmspindex = float([transpose(aspindex[1,*]),transpose(ospindex[1,*])])
conspindex = float([transpose(cspindex[1,*])])

aloh = archdat('logoh')
oloh = ohmdat('logoh')
cloh = condat('logoh')

ohmloh = float([transpose(aloh[0,*]),transpose(oloh[0,*])])
conloh = float([transpose(cloh[0,*])])

	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/stormyism/sil_a3020.ps', /color, /portrait, xs=15, ys=12
	thick = 7
	cthick = 7
	symthick = 4
	csize = 2
	legsize = 1.3

plot, ohmsil, ohmspindex, $
	/nodata, $
	color=defcolor, $
	xr = [-4,1], $
	yr = [0,8], $
	xtitle='S!I9.7!N', $
	ytitle='30-20 !7l!3m slope', $
	thick = thick, $
	xthick = thick, $
	ythick = thick, $
	charthick = cthick, $
	charsize = csize

for i = 0, n_elements(ohmsil)-1 do $
	oplot, [ohmsil[i]], [ohmspindex[i]], psym = symcat(9,thick=symthick), symsize = ohmloh[i] - 1, color=ohmcolor
oplot, consil, conspindex, psym = symcat(7), thick=thick, color=concolor

;legend, /top, /right, $
;	psym=[9,7], ['OHMs','non-masing'], $
;	charsize=legsize, thick=thick, charthick=cthick, symthick=symthick

xyouts, -1.3, 6.8, 'OHMs', charsize=1.8, /data, color=ohmcolor, charthick=5
xyouts, -3.0, 1, 'Non-masing', charsize=1.8, /data, color=concolor, charthick=5

; Draw rough locus of separation between masing and non-masing populations; somewhat arbitrary
;	regarding exact location

xlocus = [-2, -0.5]
ylocus = [2.9, 4.3]

;plots, xlocus, ylocus, linestyle=2, thick=cthick

	device,/close
	set_plot,'x'

; FORK

; SPOON_FORKII

restore,'~/Astronomy/Research/Spitzer/icelist62.sav'

; MEGA

sil_ohm = ohmdat('sil') & sil_ohm = sil_ohm[0,*]
pah62ew_ohm = ohmdat('pah62ew') & pah62ew_ohm = pah62ew_ohm[0,*]
pah62ew_ohmice = ohmdat('pah62ew_ice') & pah62ew_ohmice = pah62ew_ohmice[0,*]
pah62ewlim_ohm = ohmdat('pah62ew_lim') & pah62ewlim_ohm = pah62ewlim_ohm[0,*]
tag_ohm = transpose(ohmdat('tag'))

match, tag_ohm, icelist62, oind, count = ocount
if ocount gt 0 then pah62ew_ohm[oind] = pah62ew_ohmice[oind]

; ARCH

sil_arch = archdat('sil') & sil_arch = sil_arch[0,*]
pah62ew_arch = archdat('pah62ew') & pah62ew_arch = pah62ew_arch[0,*]
pah62ew_archice = archdat('pah62ew_ice') & pah62ew_archice = pah62ew_archice[0,*]
pah62ewlim_arch = archdat('pah62ew_lim') & pah62ewlim_arch = pah62ewlim_arch[0,*]
tag_arch = transpose(archdat('tag'))

match, tag_arch, icelist62, aind, count = acount
if acount gt 0 then pah62ew_arch[aind] = pah62ew_archice[aind]

; CONTROL

sil_con = condat('sil') & sil_con = sil_con[0,*]
pah62ew_con = condat('pah62ew') & pah62ew_con = pah62ew_con[0,*]
pah62ew_conice = condat('pah62ew_ice') & pah62ew_conice = pah62ew_conice[0,*]
pah62ewlim_con = condat('pah62ew_lim') & pah62ewlim_con = pah62ewlim_con[0,*]
tag_con = transpose(condat('tag'))

match, tag_con, icelist62, cind, count = ccount
if ccount gt 0 then pah62ew_con[cind] = pah62ew_conice[cind]

; Remove spurious data for mega034, control033

badohm = where(tag_ohm eq 'mega034')
goodohm = setdifference(indgen(n_elements(tag_ohm)),badohm)
sil_ohm = sil_ohm[goodohm]
pah62ew_ohm = pah62ew_ohm[goodohm]
pah62ewlim_ohm = pah62ewlim_ohm[goodohm]
tag_ohm = tag_ohm[goodohm]

badcon = where(tag_con eq 'control033')
goodcon = setdifference(indgen(n_elements(tag_con)),badcon)
sil_con = sil_con[goodcon]
pah62ew_con = pah62ew_con[goodcon]
pah62ewlim_con = pah62ewlim_con[goodcon]
tag_con = tag_con[goodcon]

; Read in Henrik's data

restore, '~/Astronomy/Research/Spitzer/henrik/fork.diagram.parameters.xdr'

; Remove overlapping objects in both samples

obj_ohm = ohmdat('obj') & obj_ohm = obj_ohm[goodohm]
obj_con = condat('obj') & obj_con = obj_con[goodcon]
obj_arch= transpose(archdat('obj'))

allobj = [obj_ohm,obj_arch,obj_con]
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

megalim = ['mega007']
archlim = ['arch023']
conlim  = ['control008']

match, tag_ohm, megalim, olim, olim2
nomegalim = setdifference(indgen(n_elements(tag_ohm)),olim)
match, tag_arch, archlim, alim, alim2
noarchlim = setdifference(indgen(n_elements(tag_arch)),alim)
match, tag_con, conlim, clim, clim2
noconlim = setdifference(indgen(n_elements(tag_con)),clim)

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

; Begin plotting

	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/stormyism/spoon_forkII.ps', /color, /portrait
	cs = 2.0
	ls = 7
	defcolor = fsc_color("Black")

red = fsc_color("Red")
blue = fsc_color("Blue")
lightblue = fsc_color("Cyan")
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")
green = fsc_color("Dark Green")
purple = fsc_color("Purple")
color6 = fsc_color("Rosy Brown")
grey = fsc_color("Dark Grey")

!p.multi=[0,1,1]

erase

x0 = 0.15
x1 = 0.8
x2 = 0.95

y0 = 0.15
y1 = 0.8
y2 = 0.95

!x.style = 1
!y.style = 1

plot, pah62ew_ohm, sil_ohm, $
	/nodata, $
	xtitle = '6.2 !7l!3m PAH EW [!7l!3m]', $
	ytitle = 'S!I9.7!N', $
	charsize = cs, $
	thick = ls, $
	xthick = ls, $
	ythick = ls, $
	charthick = ls, $
	/xlog, /xstyle, $
	ystyle = 1, $
	xr = [1d-3, 4], $
	yr = [1, -5], $
	yticks = 6, yminor = 0, ytickv = reverse(fillarr(1,-5,1)), ytickname = ['1','0','-1','-2','-3','-4',' '], $
;	position = [x0,y0,x1,y1], $
	color = defcolor

; Insert Henrik's divisions

;ver, 0.07
;ver, 0.45
;hor, -0.8
;hor, -2.4
hor, 0, linestyle = 1, thick = ls

; Henrik's data

oplot, pahew_a, sil_a, $
	psym = symcat(18), thick = ls, color = grey

oplot, pahew_o, sil_o, $
	psym = symcat(7), thick = ls, color = grey

oplot, pahew_s, sil_s, $
	psym = symcat(46), thick = ls, color = grey

oplot, pahew_u, sil_u, $
	psym = symcat(2), thick = ls, color = grey

;oplot, pahew_n, sil_n, $
;	psym = symcat(4), thick = ls, color = color6

; My data

oplot, pah62ew_ohm[nomegalim], sil_ohm[nomegalim], $
	psym = symcat(16), thick = ls, color = red
oplot, pah62ew_arch[noarchlim], sil_arch[noarchlim], $
	psym = symcat(16), thick = ls, color = red
oplot, pah62ew_con[noconlim], sil_con[noconlim], $
	psym = symcat(15), thick = ls, color = blue

oplot, pah62ewlim_ohm[olim], sil_ohm[olim], $
	psym = symcat(16), thick = ls, color = red
arrow, pah62ewlim_ohm[olim], sil_ohm[olim], $
	pah62ewlim_ohm[olim] * 0.7, sil_ohm[olim], $
	thick = ls, color = red, /data

oplot, pah62ewlim_arch[alim], sil_arch[alim], $
	psym = symcat(16), thick = ls, color = red
arrow, pah62ewlim_arch[alim], sil_arch[alim], $
	pah62ewlim_arch[alim] * 0.7, sil_arch[alim], $
	thick = ls, color = red, /data

oplot, pah62ewlim_con[clim], sil_con[clim], $
	psym = symcat(15), thick = ls, color = blue
arrow, pah62ewlim_con[clim], sil_con[clim], $
	pah62ewlim_con[clim] * 0.7, sil_con[clim], $
	thick = ls, color = blue, /data

legend, /top, /right, ['OHMs', 'Non-masing', 'AGN', 'Obscured', 'Starbursts', 'ULIRG/HyLIRGs'], $
	psym = [16,15,18,7,46,2], $
	color = [red, blue, grey, grey, grey, grey], $
	thick = ls, charthick = ls

;; Overplot histograms of silicate strength and PAH EW on the sides
;
;; PAH EW on top
;
;plot, indgen(10), /nodata, $
;	/xstyle, /ystyle, $
;	xr = [-3, alog10(3)], $
;	yr = [0,25], $
;	xticks = 1, $
;	xtickv = [-3.0,0.0], $
;	xtickname = replicate(' ',2), $
;	yticks = 2, $
;	ytickv = fillarr(10,0,20), $
;	ytickname = ['0','10','20'], $
;	thick = ls, $
;	charthick = ls, $
;	xthick = ls, $
;	ythick = ls, $
;	charsize = cs, $
;	position = [x0,y1,x1,y2],$
;	/noerase
;
;ohmew = alog10([pah62ew_ohm, transpose(pah62ew_arch)])
;conew = alog10([pah62ew_con])
;
;ohmew = ohmew[setdifference(indgen(n_elements(ohmew)),where(finite(ohmew,/nan)))]
;conew = conew[setdifference(indgen(n_elements(conew)),where(finite(conew,/nan)))]
;
;plothist, ohmew, datacolor = red, /overplot, bin = 0.3, thick = ls
;plothist, conew, datacolor = blue, /overplot,bin = 0.3, thick = ls
;
;xyouts, -2.5, 15, 'PAH EW', charsize = 1.0, /data, charthick = ls
;
;; Silicate on side
;
;plot, indgen(10), /nodata, $
;	/xstyle, /ystyle, $
;	xr = [0,25], $
;	yr = [1,-5], $
;	yticks = 1, $
;	ytickv = [1,-5], $
;	ytickname = replicate(' ',2), $
;	xticks = 2, $
;	xtickv = fillarr(10,0,20), $
;	xtickname = ['0','10','20'], $
;	thick = ls, $
;	charthick = ls, $
;	xthick = ls, $
;	ythick = ls, $
;	charsize = cs, $
;	position = [x1,y0,x2,y1],$
;	/noerase
;
;xyouts, 8, -4.5, 'Silicate', charsize = 1.0, /data, charthick = ls
;
;ohmsil = [sil_ohm, transpose(sil_arch)]
;consil = sil_con
;
;bs = 0.5
;plothist, ohmsil, obin, ocount, color = red, /noplot, bin = bs, thick = 3
;plothist, consil, cbin, ccount, color = blue, /noplot,bin = bs, thick = 3
;
;verbin = bs / 2. 
;
;; Vertical histogram 
;
;for i = 0, n_elements(obin) - 1 do begin
;	plots, [ocount[i], ocount[i]], [obin[i] - verbin, obin[i] + verbin],  color = red, thick = ls
;	if i ne n_elements(obin)-1 then plots, [ocount[i],ocount[i+1]], [obin[i] + verbin, obin[i] + verbin], color = red, thick = 3 else $
;		plots, [ocount[i],0], [obin[i] + verbin, obin[i] + verbin], color = red, thick = ls 
;	if i ne 0 then plots, [ocount[i],ocount[i-1]], [obin[i] - verbin, obin[i] - verbin], color = red, thick = 3 else $
;		plots, [ocount[i],0], [obin[i] - verbin, obin[i] - verbin], color = red, thick = ls 
;endfor
;
;for i = 0, n_elements(cbin) - 1 do begin
;	plots, [ccount[i], ccount[i]], [cbin[i] - verbin, cbin[i] + verbin],  color = blue, thick = ls
;	if i ne n_elements(cbin)-1 then plots, [ccount[i],ccount[i+1]], [cbin[i] + verbin, cbin[i] + verbin], color = blue, thick = 3 else $
;		plots, [ccount[i],0], [cbin[i] + verbin, cbin[i] + verbin], color = blue, thick = ls 
;	if i ne 0 then plots, [ccount[i],ccount[i-1]], [cbin[i] - verbin, cbin[i] - verbin], color = blue, thick = 3 else $
;		plots, [ccount[i],0], [cbin[i] - verbin, cbin[i] - verbin], color = blue, thick = ls 
;endfor

	device,/close
	set_plot,'x'

;; FEATURE-FEATURE
;
;ohmjunk = ohmdat('sil')
;ohmerrjunk = ohmdat('silerr')
;ohmobjjunk = ohmdat('obj')
;
;archjunk = archdat('sil')
;archerrjunk = archdat('silerr')
;archobjjunk = archdat('obj')
;
;sil10_ohm = ohmjunk[0,*]
;sil18_ohm = ohmjunk[1,*]
;sil10_ohm_err = abs(ohmerrjunk[0,*])
;sil18_ohm_err = abs(ohmerrjunk[1,*])
;
;sil10_arch = archjunk[0,*]
;sil18_arch = archjunk[1,*]
;sil10_arch_err = abs(archerrjunk[0,*])
;sil18_arch_err = abs(archerrjunk[1,*])
;
;conjunk = condat('sil')
;conerrjunk = condat('silerr')
;conobjjunk = condat('obj')
;
;sil10_con = conjunk[0,*]
;sil18_con = conjunk[1,*]
;sil10_con_err = abs(conerrjunk[0,*])
;sil18_con_err = abs(conerrjunk[1,*])
;
;; Remove control033, mega034
;
;ohmind = where(strtrim(ohmobjjunk,2) ne 'IRAS 23028+0725')
;conind = where(strtrim(conobjjunk,2) ne 'IRAS 20460+1925')
;
;sil10_ohm = sil10_ohm[ohmind]
;sil18_ohm = sil18_ohm[ohmind]
;sil10_ohm_err = sil10_ohm_err[ohmind]
;sil18_ohm_err = sil18_ohm_err[ohmind]
;ohmobj = ohmobjjunk[ohmind]
;
;sil10_con = sil10_con[conind]
;sil18_con = sil18_con[conind]
;sil10_con_err = sil10_con_err[conind]
;sil18_con_err = sil18_con_err[conind]
;conobj = conobjjunk[conind]
;
;archobj = archobjjunk
;
;!p.multi = [0,1,1]
;
;	set_plot,'ps'
;	device,filename='~/Astronomy/Research/Spitzer/stormyism/feature2_ohm.ps',/color,/portrait
;	defcolor = fsc_color("Black")
;	cs = 2.0
;	lthick = 7
;	cthick = 7
;
;red = fsc_color("Red")
;blue = fsc_color("Blue")
;yellow = fsc_color("Goldenrod")
;green = fsc_color("Dark Green")
;
;plot, indgen(10), $
;	/nodata, $
;	xtitle='S!I9.7!N', $
;	ytitle='S!I18!N', $
;	xthick = lthick, $
;	ythick = lthick, $
;	thick = lthick, $
;	charthick = cthick, $
;	charsize=cs, $
;	color=defcolor, $
;	xr=[-4,1.5], /xstyle, $
;	yr=[-2.0,1.0], /ystyle
;
;; DUSTY tracks from Sirocky 2008, measured using Dexter
;
;xstart = 1.26 & ystart = 0.67
;
;dexdir='/Applications/Dexter/'
;readcol,dexdir+'f7.gif.yellow',skipline=1,x,y,/silent
;;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=yellow, thick = lthick
;
;readcol,dexdir+'f7.gif.green1',skipline=1,x,y,/silent
;plots, [x,xstart], [y,ystart], color=green, thick = lthick
;readcol,dexdir+'f7.gif.green2',skipline=1,x,y,/silent
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=green, thick = lthick
;readcol,dexdir+'f7.gif.green3',skipline=1,x,y,/silent
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=green, thick = lthick
;
;readcol,dexdir+'f7.gif.blue2',skipline=1,x,y,/silent
;;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick
;
;readcol,dexdir+'f7.gif.black_dash2',skipline=1,x,y,/silent
;	black_dash_ind = where(x gt -4 and y gt -1)
;	x = x[black_dash_ind] & y = y[black_dash_ind]
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=2
;readcol,dexdir+'f7.gif.black_dot2',skipline=1,x,y,/silent
;	black_dash_ind = where(x gt -4 and y gt -1)
;	x = x[black_dash_ind] & y = y[black_dash_ind]
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=1
;readcol,dexdir+'f7.gif.black_solid2',skipline=1,x,y,/silent
;	black_dash_ind = where(x gt -4 and y gt -1)
;	x = x[black_dash_ind] & y = y[black_dash_ind]
;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=0
;
;readcol,dexdir+'f7.gif.red_dash2',skipline=1,x,y,/silent
;	red_dash_ind = where(x gt -4 and y gt -1)
;	x = x[red_dash_ind] & y = y[red_dash_ind]
;;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=2
;readcol,dexdir+'f7.gif.red_dot2',skipline=1,x,y,/silent
;	red_dot_ind = where(x gt -4 and y gt -1)
;	x = x[red_dot_ind] & y = y[red_dot_ind]
;;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=1
;readcol,dexdir+'f7.gif.red_solid2',skipline=1,x,y,/silent
;	red_solid_ind = where(x gt -4 and y gt -1)
;	x = x[red_solid_ind] & y = y[red_solid_ind]
;;plots, [x[sort(y)],xstart], [y[sort(y)],ystart], color=defcolor, thick = lthick, linestyle=0
;
;; ULIRG data
;
;oploterror, sil10_ohm, sil18_ohm, sil10_ohm_err, sil18_ohm_err, $
;	psym=symcat(14), color=red, errthick = lthick, errcolor=red, /nohat
;oploterror, sil10_arch, sil18_arch, sil10_arch_err, sil18_arch_err, $
;	psym=symcat(14), color=red, errthick = lthick, errcolor=red, /nohat
;oploterror, sil10_con, sil18_con, sil10_con_err, sil18_con_err, $
;	psym=symcat(7), color=blue, errthick = lthick, errcolor=blue, /nohat, thick=lthick, symsize=1.2
;
;; Starting point
;
;oplot,[xstart],[ystart], psym=symcat(46), symsize=4
;
;legend,['Y=100', 'Y=200','Y=400','clumpy'],linestyle=[0,1,2,0], $
;	color=[defcolor,defcolor,defcolor,green], $ 
;	thick = lthick, charthick = cthick, /top, /left
;
;	device,/close
;	set_plot,'x'
;
;; LE08 models
;
;; OH data
;
;o1667 = float(ohmdat('f1667'))
;o1420 = float(ohmdat('f1420'))
;otauoh = -1d * alog((o1667 + o1420)/o1420)
;
;a1667 = float(archdat('f1667'))
;a1420 = float(archdat('f1420'))
;atauoh = -1d * alog((a1667 + a1420)/a1420)
;
;c1667 = float(condat('f1667'))
;c1420 = float(condat('f1420'))
;ctauoh = -1d * alog((c1667 + c1420)/c1420)
;
;; Silicate data
;
;osil = float(ohmdat('sil'))
;osil = osil(0,*)
;otauv = -1d * osil * 18.5 / 1.085
;
;asil = float(archdat('sil'))
;asil = asil(0,*)
;atauv = -1d * asil * 18.5 / 1.085
;
;csil = float(condat('sil'))
;csil = csil(0,*)
;ctauv = -1d * csil * 18.5 / 1.085
;
;; Dust temperature data
;
;odtemp = float(ohmdat('dtemp'))
;odtemp = odtemp(0,*)
;
;adtemp = float(archdat('dtemp'))
;adtemp = adtemp(0,*)
;
;cdtemp = float(condat('dtemp'))
;cdtemp = cdtemp(0,*)
;
;; Object IDs
;
;otag = ohmdat('tag') & oobj = ohmdat('obj')
;atag = archdat('tag') & aobj = archdat('obj')
;ctag = condat('tag') & cobj = condat('obj')
;
;; Cull for objects with bad data in one or more parameters
;
;ogoodind = where(finite(otauoh) eq 1 and otauv gt 0 and odtemp gt 0)
;agoodind = where(finite(atauoh) eq 1 and atauv gt 0 and adtemp gt 0)
;cgoodind = where(ctauv gt 0 and cdtemp gt 0)
;
;alltauoh = [otauoh(ogoodind),atauoh(agoodind)]
;alltauv = [otauv(ogoodind),atauv(agoodind)]
;alldtemp = [odtemp(ogoodind),adtemp(agoodind)]
;alltags = [otag(ogoodind),atag(agoodind)]
;allobjs = [oobj(ogoodind),aobj(agoodind)]
;
;contauoh = ctauoh[cgoodind]
;contauv = ctauv[cgoodind]
;condtemp = cdtemp[cgoodind]
;contags = ctag[cgoodind]
;conobjs = cobj[cgoodind]
;
;; Color indices
;
;violet = fsc_color("Black")		; -4.0 < t_oh < -3.5
;royalblue = fsc_color("Royal Blue")	; -3.5 < t_oh < -3.0
;cyan = fsc_color("Cyan")		; -3.0 < t_oh < -2.5
;seagreen = fsc_color("Sea Green")	; -2.5 < t_oh < -2.0
;green = fsc_color("Green")		; -2.0 < t_oh < -1.5
;yellow = fsc_color("Yellow")		; -1.5 < t_oh < -1.0
;orange = fsc_color("Orange")		; -1.0 < t_oh < -0.5
;defcolor = fsc_color("White")
;white = fsc_color("White")
;black = fsc_color("Black")
;red = fsc_color("Red")
;
;royalblueind = where(alltauoh gt -3.5 and alltauoh lt -3.0)
;cyanind = where(alltauoh gt -3.0 and alltauoh lt -2.5)
;seagreenind = where(alltauoh gt -2.5 and alltauoh lt -2.0)
;greenind = where(alltauoh gt -2.0 and alltauoh lt -1.5)
;yellowind = where(alltauoh gt -1.5 and alltauoh lt -1.0)
;orangeind = where(alltauoh gt -1.0 and alltauoh lt -0.5)
;defind = where(alltauoh gt -0.5 and alltauoh lt -0.0)
;
;croyalblueind = where(contauoh gt -3.5 and contauoh lt -3.0)
;ccyanind = where(contauoh gt -3.0 and contauoh lt -2.5)
;cseagreenind = where(contauoh gt -2.5 and contauoh lt -2.0)
;cgreenind = where(contauoh gt -2.0 and contauoh lt -1.5)
;cyellowind = where(contauoh gt -1.5 and contauoh lt -1.0)
;corangeind = where(contauoh gt -1.0 and contauoh lt -0.5)
;cdefind = where(contauoh gt -0.5 and contauoh lt -0.0)
;
;colors=[violet,royalblue,cyan,seagreen,green,yellow,orange, defcolor, defcolor]
;barcolors = reverse(colors[0:n_elements(colors)-2])
;
;; Load in contours from Lockett & Elitzur 2008
;
;lockfile = '~/Astronomy/Research/Spitzer/ohm/lockett.txt'
;readcol, lockfile, tdust, tauv, taumaser, format='i,i,f', skipline=1, /silent
;tdust_arr = tdust(rem_dup(tdust))
;tauv_arr = tauv(rem_dup(tauv))
;taumaser_arr = transpose(reform(taumaser,10,9))
;
;!p.multi=[0,1,2]
;
;	set_plot,'ps'
;	device, filename='~/Astronomy/Research/Spitzer/stormyism/lockettplot.ps', /portrait, /color, xs=14, ys=26, xoff=0, yoff=0
;	defcolor=fsc_color("Black")
;	cs = 1.5
;	barcs = 1.0
;	th = 4
;
;contour, taumaser_arr, tdust_arr, tauv_arr,$
;	levels = fillarr(0.5,-4,0), $
;	path_xy=xy, path_info=info, /path_data_coords
;
;contour, taumaser_arr, tdust_arr, tauv_arr,$
;	levels = fillarr(0.5,-4,0), $
;	c_colors = colors, $
;	/fill, $
;;	xtitle = 'T!Idust!N [K]', $
;	ytitle = '!7s!3!IV!N', $
;;	xrange = [30,180], /xstyle, $
;;	yrange = [0,150], ystyle=9, $
;	xrange = [30,120], /xstyle, $
;	yrange = [0,60], ystyle=9, $
;	yticks=6, $
;	ytickv = [0,10,20,30,40,50,60], $
;	charsize = cs, $
;	charthick = th, $
;	thick = th, $
;	xthick = th, $
;	ythick = th, $
;	xtickname=replicate(' ',5), $
;	position=[0.15, 0.45, 0.85, 0.80], $
;	color=defcolor
;
;; Colorbar
;
;	for i=0,7 do polyfill, 0.15 + 0.7*[i, i, i+1, i+1, i] / 8., $
;		[0.88, 0.95, 0.95, 0.88, 0.88], $
;		/norm, color = barcolors[i]
;	plots, [0.15, 0.15, 0.85, 0.85, 0.15], [0.88, 0.95, 0.95, 0.88, 0.88], /norm, color = defcolor, thick = th
;	for i = 0,8 do begin
;		plots, [0.15 + 0.7*i/8, 0.15+0.7*i/8], [0.88, 0.95], /norm, color=defcolor, thick = th
;		xyouts, 0.13 + 0.7*i/8, 0.85, string(0. - 0.5*i,format='(f4.1)'), /norm, $
;			color=defcolor, charsize=barcs, charthick = th
;	endfor
;	xyouts, 0.87, 0.90, '!7s!3(1667)', /normal, color=defcolor, charsize=barcs, charthick=th
;
;if defind(0) ne -1 then oplot,alldtemp(defind), alltauv(defind), psym=symcat(14), color=black, symsize = 1.3
;if royalblueind(0) ne -1 then oplot,alldtemp(royalblueind), alltauv(royalblueind), psym=symcat(14), color=black, symsize = 1.3
;if cyanind(0) ne -1 then oplot,alldtemp(cyanind), alltauv(cyanind), psym=symcat(14), color=black, symsize = 1.3
;if seagreenind(0) ne -1 then oplot,alldtemp(seagreenind), alltauv(seagreenind), psym=symcat(14), color=black, symsize = 1.3
;if greenind(0) ne -1 then oplot,alldtemp(greenind), alltauv(greenind), psym=symcat(14), color=black, symsize = 1.3
;if yellowind(0) ne -1 then oplot,alldtemp(yellowind), alltauv(yellowind), psym=symcat(14), color=black, symsize = 1.3
;if orangeind(0) ne -1 then oplot,alldtemp(orangeind), alltauv(orangeind), psym=symcat(14), color=black, symsize = 1.3
;
;if royalblueind(0) ne -1 then oplot,alldtemp(royalblueind), alltauv(royalblueind), psym=symcat(14), color=royalblue, symsize=1.0
;if cyanind(0) ne -1 then oplot,alldtemp(cyanind), alltauv(cyanind), psym=symcat(14), color=cyan, symsize=1.0
;if seagreenind(0) ne -1 then oplot,alldtemp(seagreenind), alltauv(seagreenind), psym=symcat(14), color=seagreen, symsize=1.0
;if greenind(0) ne -1 then oplot,alldtemp(greenind), alltauv(greenind), psym=symcat(14), color=green, symsize=1.0
;if yellowind(0) ne -1 then oplot,alldtemp(yellowind), alltauv(yellowind), psym=symcat(14), color=yellow, symsize=1.0
;if orangeind(0) ne -1 then oplot,alldtemp(orangeind), alltauv(orangeind), psym=symcat(14), color=orange, symsize=1.0
;if defind(0) ne -1 then oplot,alldtemp(defind), alltauv(defind), psym=symcat(14), color=white, symsize=1.0
;
;oplot, condtemp, contauv, psym=symcat(7), color=black, thick=5, symsize=1.3
;
;; Add average error bar
;
;oploterror, [98], [50], [20], [5], thick=3
;
;axis, yaxis=1, yrange=[0,150] * 1.086, /save, color = defcolor, $
;	ytitle='A!IV!N', $
;	charsize = cs, charthick = th
;
;; Bottom plot - LE08 with DUSTY models
;
;; OH data
;
;a1667 = float(archdat('f1667'))
;a1420 = float(archdat('f1420'))
;atauoh = -1d * alog((a1667 + a1420)/a1420)
;
;o1667 = float(ohmdat('f1667'))
;o1420 = float(ohmdat('f1420'))
;otauoh = -1d * alog((o1667 + o1420)/o1420)
;
;ohm_tauoh = [transpose(atauoh), transpose(otauoh)]
;
;c1667 = float(condat('f1667'))
;c1420 = float(condat('f1420'))
;con_tauoh = -1d * alog((c1667 + c1420)/c1420)
;
;; DUSTY parameter fits
;
;restore, '~/Astronomy/Research/Spitzer/dusty_pahfit_results.sav'	; [arch, mega]
;
;ohmobj = [transpose(archdat('obj')), transpose(ohmdat('obj'))]
;conobj = [transpose(condat('obj'))]
;
;; Cull for objects with bad data in one or more parameters
;
;ohm_goodind = where(finite(ohm_tauoh) eq 1 and ohm_tauoh ne 0 and ohm_tauvarr gt 0 and ohm_tdustarr gt 0)
;con_goodind = where(finite(con_tauoh) eq 1 and con_tauvarr gt 0 and con_tdustarr gt 0)
;
;ohm_tauoh = ohm_tauoh[ohm_goodind]
;ohm_tauv = ohm_tauvarr[ohm_goodind]
;ohm_tdust = ohm_tdustarr[ohm_goodind]
;ohm_tags = ohmtags[ohm_goodind]
;ohm_objs = ohmobj[ohm_goodind]
;
;con_tauoh = con_tauoh[con_goodind]
;con_tauv = con_tauvarr[con_goodind]
;con_tdust = con_tdustarr[con_goodind]
;con_tags = contags[con_goodind]
;con_objs = conobj[con_goodind]
;
;; Color indices
;
;violet = fsc_color("Black")		; -4.0 < t_oh < -3.5
;royalblue = fsc_color("Royal Blue")	; -3.5 < t_oh < -3.0
;cyan = fsc_color("Cyan")		; -3.0 < t_oh < -2.5
;seagreen = fsc_color("Sea Green")	; -2.5 < t_oh < -2.0
;green = fsc_color("Green")		; -2.0 < t_oh < -1.5
;yellow = fsc_color("Yellow")		; -1.5 < t_oh < -1.0
;orange = fsc_color("Orange")		; -1.0 < t_oh < -0.5
;defcolor = fsc_color("Black")
;white = fsc_color("White")
;black = fsc_color("Black")
;red = fsc_color("Red")
;
;royalblueind = where(ohm_tauoh gt -3.5 and ohm_tauoh lt -3.0)
;cyanind = where(ohm_tauoh gt -3.0 and ohm_tauoh lt -2.5)
;seagreenind = where(ohm_tauoh gt -2.5 and ohm_tauoh lt -2.0)
;greenind = where(ohm_tauoh gt -2.0 and ohm_tauoh lt -1.5)
;yellowind = where(ohm_tauoh gt -1.5 and ohm_tauoh lt -1.0)
;orangeind = where(ohm_tauoh gt -1.0 and ohm_tauoh lt -0.5)
;defind = where(ohm_tauoh gt -0.5 and ohm_tauoh lt -0.0)
;
;croyalblueind = where(con_tauoh gt -3.5 and con_tauoh lt -3.0)
;ccyanind = where(con_tauoh gt -3.0 and con_tauoh lt -2.5)
;cseagreenind = where(con_tauoh gt -2.5 and con_tauoh lt -2.0)
;cgreenind = where(con_tauoh gt -2.0 and con_tauoh lt -1.5)
;cyellowind = where(con_tauoh gt -1.5 and con_tauoh lt -1.0)
;corangeind = where(con_tauoh gt -1.0 and con_tauoh lt -0.5)
;cdefind = where(con_tauoh gt -0.5 and con_tauoh lt -0.0)
;
;colors=[violet,royalblue,cyan,seagreen,green,yellow,orange, white, white]
;barcolors = reverse(colors[0:n_elements(colors)-2])
;
;; Load in contours from Lockett & Elitzur 2008
;
;lockfile = '~/Astronomy/Research/Spitzer/ohm/lockett.txt'
;readcol, lockfile, tdust, tauv, taumaser, format='i,i,f', skipline=1, /silent
;tdust_arr = tdust(rem_dup(tdust))
;tauv_arr = tauv(rem_dup(tauv))
;taumaser_arr = transpose(reform(taumaser,10,9))
;
;contour, taumaser_arr, tdust_arr, tauv_arr,$
;	levels = fillarr(0.5,-4,0), $
;	path_xy=xy, path_info=info, /path_data_coords
;
;contour, taumaser_arr, tdust_arr, tauv_arr,$
;	levels = fillarr(0.5,-4,0), $
;	c_colors = colors, $
;	/fill, $
;	xtitle = 'T!Idust!N [K]', $
;	ytitle = '!7s!3!IV!N', $
;	xrange = [30,150], /xstyle, $
;	yrange = [0,450], ystyle=9, $
;	charsize = cs, $
;	charthick = th, $
;	thick = th, $
;	xthick = th, $
;	ythick = th, $
;	position=[0.15, 0.10, 0.85, 0.45], $
;	color=defcolor
;
;
;if defind(0) ne -1 then oplot,ohm_tdust(defind), ohm_tauv(defind), psym=symcat(14), color=black, symsize=1.3
;if royalblueind(0) ne -1 then oplot,ohm_tdust(royalblueind), ohm_tauv(royalblueind), psym=symcat(14), color=black, symsize=1.3
;if cyanind(0) ne -1 then oplot,ohm_tdust(cyanind), ohm_tauv(cyanind), psym=symcat(14), color=black, symsize=1.3
;if seagreenind(0) ne -1 then oplot,ohm_tdust(seagreenind), ohm_tauv(seagreenind), psym=symcat(14), color=black, symsize=1.3
;if greenind(0) ne -1 then oplot,ohm_tdust(greenind), ohm_tauv(greenind), psym=symcat(14), color=black, symsize=1.3
;if yellowind(0) ne -1 then oplot,ohm_tdust(yellowind), ohm_tauv(yellowind), psym=symcat(14), color=black, symsize=1.3
;if orangeind(0) ne -1 then oplot,ohm_tdust(orangeind), ohm_tauv(orangeind), psym=symcat(14), color=black, symsize=1.3
;
;if royalblueind(0) ne -1 then oplot,ohm_tdust(royalblueind), ohm_tauv(royalblueind), psym=symcat(14), color=royalblue, symsize=0.8
;if cyanind(0) ne -1 then oplot,ohm_tdust(cyanind), ohm_tauv(cyanind), psym=symcat(14), color=cyan, symsize=1.0
;if seagreenind(0) ne -1 then oplot,ohm_tdust(seagreenind), ohm_tauv(seagreenind), psym=symcat(14), color=seagreen, symsize=1.0
;if greenind(0) ne -1 then oplot,ohm_tdust(greenind), ohm_tauv(greenind), psym=symcat(14), color=green, symsize=1.0
;if yellowind(0) ne -1 then oplot,ohm_tdust(yellowind), ohm_tauv(yellowind), psym=symcat(14), color=yellow, symsize=1.0
;if orangeind(0) ne -1 then oplot,ohm_tdust(orangeind), ohm_tauv(orangeind), psym=symcat(14), color=orange, symsize=1.0
;if defind(0) ne -1 then oplot,ohm_tdust(defind), ohm_tauv(defind), psym=symcat(14), color=white, symsize=1.0
;
;oplot, con_tdust, con_tauv, psym=symcat(7), color=black, thick=5, symsize=1.3
;
;; Average error bar for the sample
;
;oploterror, [127], [403], [20], [20], thick=3
;
;axis, yaxis=1, yrange=[0,450] * 1.086, /save, color = defcolor, $
;	ytitle='A!IV!N', $
;	charsize = cs, charthick = th
;
;	device,/close
;	set_plot,'x'

end
