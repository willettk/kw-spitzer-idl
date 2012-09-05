pro spoon_forkII, ps = ps, label = label, stop = stop
;+
; NAME:
;       
;	SPOON_FORKII
;
; PURPOSE:
;
;	Plot PAH EW vs. silicate optical depth for entire sample (Paper II)
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
;	PS - 		hard copy of plot
;
;	LABEL - 	label with tags (eg, 'mega001') for each object in my data
;
;	STOP - 		stop program at end
;
; EXAMPLE:
;
;	IDL> spoon_forkII, /ps
;
; NOTES:
;
;	This diagram follows the classification scheme of Spoon et al., ApJL, 2007 (astro-ph 0611918).
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;-

; Read in data

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

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/spoon_forkII.ps', /color, /portrait
	cs = 1.5
	ls = 4
	defcolor = fsc_color("Black")
endif else begin
	cs = 2
	ls = 1
	defcolor = fsc_color("White")
endelse

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
	position = [x0,y0,x1,y1], $
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

if keyword_set(label) then begin
	xyouts, pah62ew_ohm,  sil_ohm-0.1,  tag_ohm, color=red
	xyouts, pah62ew_con,  sil_con-0.1,  tag_con, color=blue
	xyouts, pah62ew_arch, sil_arch-0.1, tag_arch, color=red
	xyouts, pah62ewlim_con[clim],  sil_con[clim]-0.1,  tag_con[clim], color=blue
	xyouts, pah62ewlim_arch[alim], sil_arch[alim]-0.1, tag_arch[alim], color=red
	xyouts, pah62ewlim_arch[alim], sil_arch[alim]-0.1, tag_arch[alim], color=red
	xyouts, pah62ewlim_arch[olim], sil_arch[olim]-0.1, tag_arch[olim], color=red
endif

legend, /top, /right, ['OHMs', 'Non-masing', 'AGN', 'Obscured', 'Starbursts', 'ULIRG/HyLIRGs'], $
	psym = [16,15,18,7,46,2], $
	color = [red, blue, grey, grey, grey, grey], $
	thick = ls, charthick = ls

; Overplot histograms of silicate strength and PAH EW on the sides

; PAH EW on top

plot, indgen(10), /nodata, $
	/xstyle, /ystyle, $
	xr = [-3, alog10(3)], $
	yr = [0,25], $
	xticks = 1, $
	xtickv = [-3.0,0.0], $
	xtickname = replicate(' ',2), $
	yticks = 2, $
	ytickv = fillarr(10,0,20), $
	ytickname = ['0','10','20'], $
	thick = ls, $
	charthick = ls, $
	xthick = ls, $
	ythick = ls, $
	charsize = cs, $
	position = [x0,y1,x1,y2],$
	/noerase

ohmew = alog10([pah62ew_ohm, transpose(pah62ew_arch)])
conew = alog10([pah62ew_con])

ohmew = ohmew[setdifference(indgen(n_elements(ohmew)),where(finite(ohmew,/nan)))]
conew = conew[setdifference(indgen(n_elements(conew)),where(finite(conew,/nan)))]

plothist, ohmew, datacolor = red, /overplot, bin = 0.3, thick = ls
plothist, conew, datacolor = blue, /overplot,bin = 0.3, thick = ls

xyouts, -2.5, 15, 'PAH EW', charsize = 1.0, /data, charthick = ls

; Silicate on side

plot, indgen(10), /nodata, $
	/xstyle, /ystyle, $
	xr = [0,25], $
	yr = [1,-5], $
	yticks = 1, $
	ytickv = [1,-5], $
	ytickname = replicate(' ',2), $
	xticks = 2, $
	xtickv = fillarr(10,0,20), $
	xtickname = ['0','10','20'], $
	thick = ls, $
	charthick = ls, $
	xthick = ls, $
	ythick = ls, $
	charsize = cs, $
	position = [x1,y0,x2,y1],$
	/noerase

xyouts, 8, -4.5, 'Silicate', charsize = 1.0, /data, charthick = ls

ohmsil = [sil_ohm, transpose(sil_arch)]
consil = sil_con

bs = 0.5
plothist, ohmsil, obin, ocount, color = red, /noplot, bin = bs, thick = 3
plothist, consil, cbin, ccount, color = blue, /noplot,bin = bs, thick = 3

verbin = bs / 2. 

; Vertical histogram 

for i = 0, n_elements(obin) - 1 do begin
	plots, [ocount[i], ocount[i]], [obin[i] - verbin, obin[i] + verbin],  color = red, thick = ls
	if i ne n_elements(obin)-1 then plots, [ocount[i],ocount[i+1]], [obin[i] + verbin, obin[i] + verbin], color = red, thick = 3 else $
		plots, [ocount[i],0], [obin[i] + verbin, obin[i] + verbin], color = red, thick = ls 
	if i ne 0 then plots, [ocount[i],ocount[i-1]], [obin[i] - verbin, obin[i] - verbin], color = red, thick = 3 else $
		plots, [ocount[i],0], [obin[i] - verbin, obin[i] - verbin], color = red, thick = ls 
endfor

for i = 0, n_elements(cbin) - 1 do begin
	plots, [ccount[i], ccount[i]], [cbin[i] - verbin, cbin[i] + verbin],  color = blue, thick = ls
	if i ne n_elements(cbin)-1 then plots, [ccount[i],ccount[i+1]], [cbin[i] + verbin, cbin[i] + verbin], color = blue, thick = 3 else $
		plots, [ccount[i],0], [cbin[i] + verbin, cbin[i] + verbin], color = blue, thick = ls 
	if i ne 0 then plots, [ccount[i],ccount[i-1]], [cbin[i] - verbin, cbin[i] - verbin], color = blue, thick = 3 else $
		plots, [ccount[i],0], [cbin[i] - verbin, cbin[i] - verbin], color = blue, thick = ls 
endfor

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Histogram luminosities

;if keyword_set(ps) then begin
;	set_plot,'ps'
;	device,filename='~/Astronomy/Research/Spitzer/papers/pah62lum_hist.ps', /color, /portrait
;	cs = 1.5
;	ls = 3
;	defcolor = fsc_color("Black")
;endif else begin
;	cs = 2
;	ls = 1
;	defcolor = fsc_color("White")
;endelse
;
;!p.multi=[0,1,1]
;
;plot, indgen(10), /nodata, /xstyle, /ystyle, $
;	xr = [6,11], $
;	yr = [0,20], $
;	xtitle = 'log (6.2 !7l!3m L!IPAH!N/L'+sunsymbol()+')', $
;	ytitle = 'Count', $
;	thick = ls, $
;	charthick = ls, $
;	charsize = cs
;
;ohmlum = [pah62lum_ohm,transpose(pah62lum_arch)]
;badohmlum = where(finite(ohmlum,/nan))
;ohmlum = ohmlum[setdifference(indgen(n_elements(ohmlum)),badohmlum)]
;conlum = pah62lum_con
;badconlum = where(finite(conlum,/nan))
;conlum = conlum[setdifference(indgen(n_elements(conlum)),badconlum)]
;
;plothist, ohmlum, color = red, /overplot, bin = 0.4, thick = 3;, /fill, fcolor = red, /fline, forientation = 45
;plothist, conlum, color = blue, /overplot, bin = 0.4, thick = 3;, /fill, fcolor = blue, /fline, forientation = -45, thick = 2
;
;if keyword_set(ps) then begin
;	device,/close
;	set_plot,'x'
;endif
;
if keyword_set(stop) then stop

end
