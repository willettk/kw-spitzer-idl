pro sfr_plotI, ps = ps, stop = stop
;+
; NAME:
;       
;	SFR_PLOTI
;
; PURPOSE:
;
;	Print SFR diagrams for Ne and PAH emission in Paper I
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
;	IDL> sfr_plotI, /ps
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Sep 08
;-

mpc2cm = 3.086d24 
lsun   = 3.862d33 

; Restore Ne data

restore,'~/Astronomy/Research/Spitzer/linedata/neII.sav'
neII_ohm = line
restore,'~/Astronomy/Research/Spitzer/linedata/neIII.sav'
neIII_ohm = line

neII_ind = where(strmid(neII_ohm.tag,0,4) eq 'mega')
neIII_ind = where(strmid(neIII_ohm.tag,0,4) eq 'mega')

neII_tag = neII_ohm.tag[neII_ind]
tag_list = transpose(ohmdat('tag'))
dl_list = ohmdat('dl')
match, tag_list, neII_tag, a2, bb
dl2 = dl_list[a2]

neII_flux = neII_ohm.flux[neII_ind]
neII_lum = alog10(neII_flux * (dl2 * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

neIII_tag = neIII_ohm.tag[neIII_ind]
tag_list = transpose(ohmdat('tag'))
dl_list = ohmdat('dl')
match, tag_list, neIII_tag, a3, bb
dl3 = dl_list[a3]

neIII_flux = neIII_ohm.flux[neIII_ind]
neIII_lum = alog10(neIII_flux * (dl3 * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

; Objects in which both lines are measured

match, neII_tag, neIII_tag, i2, i3
match, neII_tag[i2], tag_list, ii2, ii3
dlboth = dl_list[ii3]
neboth_flux = neII_ohm.flux[i2] + neIII_ohm.flux[i3]
neboth_lum = alog10(neboth_flux * (dlboth * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

lfir = ohmdat('lfir')
lfir2 = lfir[a2]
lfir3 = lfir[a3]
lfir_both = lfir[ii3]

; Linear fit to Ne data

xarr = fillarr(1d-2,0,20)

expr='p[0]+x*p[1]'
start = [0,1]

fit2 = mpfitexpr(expr,lfir2, neII_lum,0.1*neII_lum,start,/quiet, perror = err2)
fitboth = mpfitexpr(expr,lfir_both, neboth_lum,0.1*neboth_lum,start,/quiet, perror = errboth)

; Plot results

red = fsc_color("Red")
green = fsc_color("Forest Green")

plotname='~/Astronomy/Research/Spitzer/papers/sfr_plotI.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color, /portrait
	cthick = 3
	lthick = 3
	defcolor=fsc_color("Black")
endif else begin
	cthick = 1
	lthick = 1
	defcolor=fsc_color("White")
endelse

plot, lfir, neII_lum, $
	/nodata, $
	xtitle = 'log [L!IFIR!N/L'+sunsymbol()+']', $
	ytitle = 'log [L!INe!N/L'+sunsymbol()+']', $
	charsize = 1.5, $
	xr = [11.4, 12.5], $
	yr = [7.5, 10.0], $
	color=defcolor, $
	thick = lthick, $
	charthick = cthick

oplot, lfir2, neII_lum, $
	psym = symcat(16), $
	color = red

oplot, lfir_both, neboth_lum, $
	psym = symcat(7), $
	color = green, thick = lthick

oplot, xarr, fit2[0] + xarr*fit2[1], color=red, thick = lthick
oplot, xarr, fitboth[0] + xarr*fitboth[1], color=green, thick = lthick

oplot, xarr, -3.44 + xarr * 1.01, color = red, linestyle = 1, thick = lthick
oplot, xarr, -2.78 + xarr * 0.98, color = green, linestyle = 1, thick = lthick

legend, /top, /left, $
	['NeII','NeII + NeIII', 'HK07 NeII','HK07 NeII + NeIII'], $
	linestyle = [0,0,1,1], $
	color = [red, green, red, green], $
	thick = lthick, charthick = cthick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Plot the mean star formation rate using both Ne and PAH methods

; Neon

nelumcgs = lsun * (10d)^neboth_lum
sfr_neon = nelumcgs * 4.73d-41

; PAH (spline)

p6 = ohmdat('pah62lum')
p11 = ohmdat('pah11lum')

ind = setdifference(indgen(26),[5,25])		; Remove mega007, mega034 
p6 = p6[ind] & p11 = p11[ind]

p6lum = lsun * (10d)^p6
p11lum = lsun * (10d)^p11
sfr_pah = (p6lum + p11lum) * 1.18d-41

; PAH (PAHFIT)

pf6 = ohmdat('pahfit62lum')
pf11 = ohmdat('pahfit11lum')

ind = setdifference(indgen(26),[5,25])		; Remove mega007, mega034 
pf6 = pf6[ind] & pf11 = pf11[ind]

pf6lum = lsun * (10d)^pf6
pf11lum = lsun * (10d)^pf11
sfr_pahfit = (pf6lum + pf11lum) * 1.18d-41

print,'SFR PAH spline = ',mean(sfr_pah), ' +_',stddev(sfr_pah)
print,'SFR neon =       ',mean(sfr_neon), ' +_',stddev(sfr_neon)
print,'SFR PAH PAHFIT = ',mean(sfr_pahfit), ' +_',stddev(sfr_pahfit)
print,'NeII fit ',fit2, err2
print,'Ne fit ',fitboth, errboth

if keyword_set(stop) then stop
end
