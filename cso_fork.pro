pro cso_fork, ps = ps, label = label, stop = stop
;+
; NAME: 
;       CSO_FORK 
;
; PURPOSE:
; 	Plot silicate 9.7 optical depth vs. PAH 6.2 um EW
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
;
; INPUTS:
;	
; OUTPUTS:
;
; KEYWORDS:
;
; REQUIRES:
;
; EXAMPLE:
;
;	IDL> cso_fork, /ps
; 
; NOTES:
;
;	This diagram follows the classification scheme of Spoon et al., ApJL, 2006 (astro-ph 0611918).
;
;	Different from SPOON.pro - both the OHMs and the control sample are plotted as ULIRGs, while
;	the CSOs are plotted as the red points
;
; MODIFICATION HISTORY:
;
;	Written by KW, Oct 07
;-

restore,'~/Astronomy/Research/Spitzer/icelist62.sav'
mpc2cm = 3.086d24		; Mpc to cm conversion
lsun = 3.862d33			; L_sun in erg/s

; CSO data

restore, '~/Astronomy/Research/Spitzer/cso/data/idl_sav/sil_pah_cso.sav'
sil_cso = sil
pah62ew_cso = pah62ew
pah62flux_cso = pah62flux
dl_cso = double(csodat('dl')) * mpc2cm
pah62lum_cso = alog10(4d * !dpi * dl_cso^2 * pah62flux_cso * 1d7 / lsun)
tag_cso = csofiles_silpah

match, tag_cso, icelist62, oind, oiceind, count = ocount
if ocount gt 0 then pah62ew_cso[oind] = pah62ew_ice[oind]

; Use ice-corrected EW

match, tag_cso, icelist62, oind, oiceind, count = ocount
if ocount gt 0 then pah62ew_cso[oind] = pah62ew_ice[oind]

; MEGA

restore, '~/Astronomy/Research/Spitzer/ohm/data/idl_sav/sil_pah_ohm.sav'
sil_ohm = sil
pah62ew_ohm = pah62ew
pah62flux_ohm = pah62flux
pah62ewlim_ohm = pah62ew_lim
dl_ohm = double(ohmdat('dl')) * mpc2cm
pah62lum_ohm = alog10(4d * !dpi * dl_ohm^2 * pah62flux_ohm * 1d7 / lsun)
tag_ohm = ohmfiles_silpah

match, tag_ohm, icelist62, oind, oiceind, count = ocount
if ocount gt 0 then pah62ew_ohm[oind] = pah62ew_ice[oind]

badohm = where(tag_ohm eq 'mega034')
goodohm = setdifference(indgen(n_elements(tag_ohm)),badohm)
sil_ohm = sil_ohm[goodohm]
pah62ew_ohm = pah62ew_ohm[goodohm]
pah62flux_ohm = pah62flux_ohm[goodohm]
pah62ewlim_ohm = pah62ewlim_ohm[goodohm]
pah62lum_ohm = pah62lum_ohm[goodohm]
tag_ohm = tag_ohm[goodohm]

; ARCH

restore, '~/Astronomy/Research/Spitzer/archived/data/idl_sav/sil_pah_arch.sav'
sil_arch = sil
pah62ew_arch = pah62ew
pah62flux_arch = pah62flux
pah62ewlim_arch = pah62ew_lim
dl_arch = double(archdat('dl')) * mpc2cm
pah62lum_arch = alog10(4d * !dpi * dl_arch^2 * pah62flux_arch * 1d7 / lsun)
tag_arch = archfiles_silpah

match, tag_arch, icelist62, aind, aiceind, count = acount
if acount gt 0 then pah62ew_arch[aind] = pah62ew_ice[aind]

; Read in Henrik's data

restore, '~/Astronomy/Research/Spitzer/henrik/fork.diagram.parameters.xdr'

; Remove overlapping objects in both samples

obj_ohm = ohmdat('obj') & obj_ohm = obj_ohm[goodohm]
obj_arch= transpose(archdat('obj'))

allobj = [obj_ohm,obj_arch]
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

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/cso_fork.ps', /portrait, /color
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
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")
green = fsc_color("Dark Green")
purple = fsc_color("Purple")

; Plot the fork diagram in log space

!p.multi=[0,1,1]

plot, pah62ew_cso, sil_cso, $
	psym = 4, $
	xtitle = '6.2 !7l!3m PAH EW [!7l!3m]', $
	ytitle = 'Silicate strength', $
	charsize = cs, $
	thick = ls, $
	xthick = ls, $
	ythick = ls, $
	charthick = ls, $
	/xlog, /xstyle, xr = [1d-3, 2], $
	/ystyle, yr = [1,-4.7], $
	/nodata, $
	color = defcolor

; Insert Henrik's divisions

;ver, 0.07
;ver, 0.45
;hor, -0.8
;hor, -2.4
hor, 0, linestyle = 1, thick=2

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
	psym = symcat(14), thick = ls, color = blue

; OHM data

oplot, pah62ew_ohm, sil_ohm, $
	psym = symcat(2), thick = ls, color = defcolor
oplot, pah62ew_arch, sil_arch, $
	psym = symcat(2), thick = ls, color = defcolor

; CSO data

have_pah = [1,2,3,5]
oplot, pah62ew_cso(have_pah), sil_cso(have_pah), $
	psym = symcat(15), thick = ls, color=red

pks_limit = 0.013
pks_sil = sil_cso[4]

eleven_limit = 0.009
eleven_sil = sil_cso[0]

oplot, [pks_limit],[pks_sil], $					; Overplot limits for PKS 1413+135 and 1146+59
	psym = symcat(15), thick = ls, color=red

oplot, [eleven_limit],[eleven_sil], $					
	psym = symcat(15), thick = ls, color=red

arrow, pks_limit,pks_sil,pks_limit * 0.6 ,pks_sil, thick = ls, color = red, /data
arrow, eleven_limit,eleven_sil,eleven_limit * 0.6 ,eleven_sil, thick = ls, color = red, /data

if keyword_set(label) then begin
	labelsize = 1

	xyouts, 0.032, 0.25, 'OQ 208', color=red, charthick=ls, charsize = labelsize
	xyouts, 0.0015, -0.40, '4C 12.50', color=red, charthick=ls, charsize = labelsize
	xyouts, 0.06, -0.90, '4C 31.04', color=red, charthick=ls, charsize = labelsize
	xyouts, 0.60, -2.0, 'NGC 5793', color=red, charthick=ls, charsize = labelsize
	plots,[0.55,0.70],[-1.50,-1.90],thick=3,color=red
	xyouts, eleven_limit * 0.3, eleven_sil - 0.3, '1146+59', color=red, charthick=ls, charsize = labelsize
	xyouts, pks_limit * 0.9,pks_sil-0.2, 'PKS 1413+135', color=red, charthick=ls, charsize = labelsize
endif

legend, /top, /right, ['CSO','AGN', 'Obscured AGN', 'Normal', 'Starburst', 'ULIRG/HyLIRGs'], $
	psym = [15,18,7,14,46,2], $
	color = [red,orange, green, blue, purple, defcolor], $
	thick = ls, charthick = ls

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop
end
