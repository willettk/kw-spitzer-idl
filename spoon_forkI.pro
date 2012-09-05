pro spoon_forkI, ps = ps, label = label, stop = stop
;+
; NAME:
;       
;	SPOON_FORKI
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
;	IDL> spoon_forkI, /ps
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

; Spline-fit data (should be used, since Henrik's data was measured using the same method)

sil_ohm = ohmdat('sil') & sil_ohm = sil_ohm[0,*]
pah62ew_ohm = ohmdat('pah62ew') & pah62ew_ohm = pah62ew_ohm[0,*]
pah62ew_ice = ohmdat('pah62ew_ice') & pah62ew_ice = pah62ew_ice[0,*]
pah62ewlim_ohm = ohmdat('pah62ew_lim') & pah62ewlim_ohm = pah62ewlim_ohm[0,*]
tag_ohm = transpose(ohmdat('tag'))

match, tag_ohm, icelist62, oind, oiceind, count = ocount
if ocount gt 0 then pah62ew_ohm[oind] = pah62ew_ice[oind]

; Remove spurious data for mega034

badohm = where(tag_ohm eq 'mega034')
goodohm = setdifference(indgen(n_elements(tag_ohm)),badohm)
sil_ohm = sil_ohm[goodohm]
pah62ew_ohm = pah62ew_ohm[goodohm]
pah62ewlim_ohm = pah62ewlim_ohm[goodohm]
tag_ohm = tag_ohm[goodohm]

; PAHFIT data

pahfit62ew_ohm = ohmdat('pahfit62ew')     & pahfit62ew_ohm = pahfit62ew_ohm[0,*]
pahfit62ew_ice = ohmdat('pahfit62ew_ice') & pahfit62ew_ice = pahfit62ew_ice[0,*]

if ocount gt 0 then pahfit62ew_ohm[oind] = pahfit62ew_ice[oind]
pahfit62ew_ohm = pahfit62ew_ohm[goodohm]

; Read in Henrik's data

restore, '~/Astronomy/Research/Spitzer/henrik/fork.diagram.parameters.xdr'

; Remove overlapping objects in both samples

obj_ohm = ohmdat('obj') & obj_ohm = obj_ohm[goodohm]

allobj = [obj_ohm]
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

limlist = ['mega007']

match, tag_ohm, limlist, ilim, ilim2
nolimits = setdifference(indgen(n_elements(tag_ohm)),ilim)

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
	device,filename='~/Astronomy/Research/Spitzer/papers/spoon_forkI.ps', /color, /portrait
	cs = 1.5
	ls = 3
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

plot, pah62ew_ohm, sil_ohm, $
	/nodata, $
	xtitle = '6.2 !7l!3m PAH EW [!7l!3m]', $
	ytitle = 'Silicate strength', $
	charsize = cs, $
	thick = ls, $
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

;oplot, pahew_n, sil_n, $
;	psym = symcat(4), thick = ls, color = color6

; My data

oplot, pah62ew_ohm[nolimits], sil_ohm[nolimits], $
	psym = symcat(16), thick = ls, color = red

oplot, pah62ewlim_ohm[ilim], sil_ohm[ilim], $
	psym = symcat(16), thick = ls, color = red
arrow, pah62ewlim_ohm[ilim], sil_ohm[ilim], $
	pah62ewlim_ohm[ilim] * 0.7, sil_ohm[ilim], $
	thick = ls, color = red, /data

; Difference between spline and PAHFIT results

mean62 = mean(double(pah62ew_ohm[nolimits]))
diffvec = abs((mean(double(pah62ew_ohm[nolimits])) - mean(double(pahfit62ew_ohm[nolimits]))))

arrow, mean62, 0.5, mean62 + diffvec, 0.5, color=defcolor, thick=4, /data

if keyword_set(label) then begin
	xyouts, pah62ew_ohm,  sil_ohm-0.1,  tag_ohm, color=red
endif

legend, /top, /right, $
	['OHMs', 'AGN', 'Obscured', 'Starbursts', 'ULIRG/HyLIRGs'], $
	psym = [16,18,7,46,2], $
	color = [red, orange, green, purple, defcolor], $
	thick = ls, charthick = ls

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
