pro spoon_fork_cso, ps = ps, label = label, stop = stop
;+
; NAME:
;       
;	SPOON_FORK_CSO
;
; PURPOSE:
;
;	Plot PAH EW vs. silicate optical depth for entire sample (Paper I)
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
;	IDL> spoon_fork_cso, /ps
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

; PAH data (spline fit should be used, since Henrik's data was measured using the same method)

sil_cso = csodat('sil') & sil_cso = sil_cso[0,*]
pah62ew_cso = csodat('pah62ew') & pah62ew_cso = pah62ew_cso[0,*]
pah62ew_csoice = csodat('pah62ew_ice') & pah62ew_csoice = pah62ew_csoice[0,*]
pah62ewlim_cso = csodat('pah62ew_lim') & pah62ewlim_cso = pah62ewlim_cso[0,*]
tag_cso = transpose(csodat('tag'))

match, tag_cso, icelist62, cind, oiceind, count = ocount
if ocount gt 0 then pah62ew_cso[cind] = pah62ew_csoice[cind]

; Remove spurious data for VII Zw 485 (cso010)

goodcso = where(tag_cso ne 'cso010' and tag_cso ne 'cso003')
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

goodohm = where(tag_ohm ne 'mega034')
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

; Begin plotting

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/spoon_fork_cso.ps', /color, /portrait
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

!p.multi=[0,1,1]

plot, pah62ew_cso, sil_cso, $
	/nodata, $
	xtitle = '6.2 !7l!3m PAH EW [!7l!3m]', $
	ytitle = 'Silicate strength', $
	charsize = cs, $
	thick = ls, $
	xthick = ls, $
	ythick = ls, $
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
	thick = ls, color = red, /data

; Difference between spline and PAHFIT results

mean62 = mean(double(pah62ew_cso[nolimits]))
diffvec = abs((mean(double(pah62ew_cso[nolimits])) - mean(double(pahfit62ew_cso[nolimits]))))

arrow, mean62, 0.5, mean62 + diffvec, 0.5, color=defcolor, thick=4, /data

if keyword_set(label) then begin
;	xyouts, 2d-3, -1,  'PKS 1413+135', color=red, charthick = ls, charsize = 1.0
;	xyouts, 6.5d-2,-0.5,'4C+31.04', color=red, charthick = ls, charsize = 1.0
;	xyouts, 5d-1, -1.9, 'NGC 5793', color=red, charthick = ls, charsize = 1.0
;	plots, [5d-1, pah62ew_cso[nolimits[2]]], [-1.85, sil_cso[nolimits[2]]], color=red, thick = ls
;	xyouts, pah62ew_cso[nolimits[3]],  sil_cso[nolimits[3]]-0.1,  obj_cso[nolimits[3]], color=red, charthick = ls, charsize = 1.0
;	xyouts, 2d-2, -1.1,    '1146+59',   color=red, charthick = ls, charsize = 1.0
;	plots, [2d-2, pah62ew_cso[nolimits[0]]], [-0.9, sil_cso[nolimits[0]]], color=red, thick = ls
;	xyouts, 9d-3,0.6,'4C+12.50',   color=red, charthick = ls, charsize = 1.0
;	plots, [1.5d-2, pah62ewlim_cso[5]], [0.4, sil_cso[5]], color=red, thick = ls
	xyouts, pah62ew_cso, sil_cso, obj_cso, color=red, charthick = ls, charsize = 1.0, /data
endif

legend, /top, /right, $
	['CSOs', 'AGN', 'Obscured AGN', 'Normal', 'Starbursts', 'ULIRG/HyLIRGs'], $
	psym = [16,18,7,15,46,2], $
	color = [red, orange, green, blue, purple, defcolor], $
	thick = ls, charthick = ls, charsize = 1.0

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
