pro sfr_plotII, ps = ps, stop = stop
;+
; NAME:
;       
;	SFR_PLOTII
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
;	IDL> sfr_plotII, /ps
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Sep 08
;	Changed L_FIR to L_IR (used in Ho+2007)  - Feb 10
;	B&W plots - Feb 10
;-

mpc2cm = 3.086d24 
lsun   = 3.862d33 

; Restore Ne data

restore,'~/Astronomy/Research/Spitzer/linedata/neII.sav'
neII = line
restore,'~/Astronomy/Research/Spitzer/linedata/neIII.sav'
neIII = line

neII_ohm_ind = where(strmid(neII.tag,0,4) eq 'mega')
neIII_ohm_ind = where(strmid(neIII.tag,0,4) eq 'mega')
neII_arch_ind = where(strmid(neII.tag,0,4) eq 'arch')
neIII_arch_ind = where(strmid(neIII.tag,0,4) eq 'arch')
neII_con_ind = where(strmid(neII.tag,0,3) eq 'con')
neIII_con_ind = where(strmid(neIII.tag,0,3) eq 'con')

;;; OHM data ;;;;

neII_ohm_tag = neII.tag[[neII_ohm_ind, neII_arch_ind]]
taglist_ohm = [transpose(ohmdat('tag')),transpose(archdat('tag'))]
dllist_ohm = [transpose(ohmdat('dl')),transpose(archdat('dl'))]
match, taglist_ohm, neII_ohm_tag, a2, bb
dl2_ohm = dllist_ohm[a2]

neIII_ohm_tag = neIII.tag[[neIII_ohm_ind, neIII_arch_ind]]
match, taglist_ohm, neIII_ohm_tag, a3, bb
dl3_ohm = dllist_ohm[a3]

neII_con_tag = neII.tag[neII_con_ind]
taglist_con = [transpose(condat('tag'))]
dllist_con = [transpose(condat('dl'))]
match, taglist_con, neII_con_tag, ca2, cbb
dl2_con = dllist_con[ca2]

neIII_con_tag = neIII.tag[neIII_con_ind]
match, taglist_con, neIII_con_tag, ca3, cbb
dl3_con = dllist_con[ca3]

neII_ohm_flux = neII.flux[[neII_ohm_ind,neII_arch_ind]]
neII_ohm_fluxerr = neII.fluxerr[[neII_ohm_ind,neII_arch_ind]]
neII_ohm_lum = alog10(neII_ohm_flux * (dl2_ohm * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

neIII_ohm_flux = neIII.flux[[neIII_ohm_ind,neIII_arch_ind]]
neIII_ohm_fluxerr = neIII.fluxerr[[neIII_ohm_ind,neIII_arch_ind]]
neIII_ohm_lum = alog10(neIII_ohm_flux * (dl3_ohm * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

neII_con_flux = neII.flux[neII_con_ind]
neII_con_fluxerr = neII.fluxerr[neII_con_ind]
neII_con_lum = alog10(neII_con_flux * (dl2_con * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

neIII_con_flux = neIII.flux[neIII_con_ind]
neIII_con_fluxerr = neIII.fluxerr[neIII_con_ind]
neIII_con_lum = alog10(neIII_con_flux * (dl3_con * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

; Objects in which both lines are measured

match, neII_ohm_tag, neIII_ohm_tag, i2, i3
match, neII_ohm_tag[i2], taglist_ohm, ii2, ii3
dl_ohm_both = dllist_ohm[ii3]
neboth_ohm_flux = neII.flux[i2] + neIII.flux[i3]
neboth_ohm_lum = alog10(neboth_ohm_flux * (dl_ohm_both * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

lir_ohm = [transpose(ohmdat('lir')),transpose(archdat('lir'))]
lir2_ohm = lir_ohm[a2]
lir3_ohm = lir_ohm[a3]
lirboth_ohm = lir_ohm[ii3]

match, neII_con_tag, neIII_con_tag, ci2, ci3
match, neII_con_tag[ci2], taglist_con, cii2, cii3
dl_con_both = dllist_con[cii3]
neboth_con_flux = neII.flux[ci2] + neIII.flux[ci3]
neboth_con_lum = alog10(neboth_con_flux * (dl_con_both * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)

lir_con = [transpose(condat('lir')),transpose(archdat('lir'))]
lir2_con = lir_con[ca2]
lir3_con = lir_con[ca3]
lirboth_con = lir_con[cii3]

; Linear fit to Ne data

xarr = fillarr(1d-2,0,20)

expr='p[0]+x*p[1]'
start = [0,1]

fit2_ohm = mpfitexpr(expr,lir2_ohm, neII_ohm_lum,0.1*neII_ohm_lum,start,/quiet, perror = err2_ohm)
fitboth_ohm = mpfitexpr(expr,lirboth_ohm, neboth_ohm_lum,0.1*neboth_ohm_lum,start,/quiet, perror = errboth_ohm)

fit2_con = mpfitexpr(expr,lir2_con, neII_con_lum,0.1*neII_con_lum,start,/quiet, perror = err2_con)
fitboth_con = mpfitexpr(expr,lirboth_con, neboth_con_lum,0.1*neboth_con_lum,start,/quiet, perror = errboth_con)

; Plot results

red = fsc_color("Red")
blue = fsc_color("Blue")

plotname='~/Astronomy/Research/Spitzer/papers/sfr_plotII.ps'

if keyword_set(ps) then begin
;	!p.font = 1
	set_plot,'ps'
	device, filename = plotname, /color, /portrait;, set_font='Times', /tt_font
	cthick = 5
	lthick = 5
	defcolor=fsc_color("Black")
endif else begin
	cthick = 1
	lthick = 1
	defcolor=fsc_color("White")
endelse

!p.multi=[0,1,1]

plot, lir_ohm, neII_ohm_lum, $
	/nodata, $
	xtitle = 'log [L!IIR!N/L'+sunsymbol()+']', $
	ytitle = 'log [L!INe!N/L'+sunsymbol()+']', $
	charsize = 1.5, $
	xr = [11.0, 13.0], $
	yr = [6.0, 10.0], $
	color=defcolor, $
	thick = lthick, $
	xthick = lthick, $
	ythick = lthick, $
	charthick = cthick

oplot, lirboth_ohm, neboth_ohm_lum, $
	psym = symcat(14), $
	color = defcolor, thick = lthick

oplot, lirboth_con, neboth_con_lum, $
	psym = symcat(7), $
	color = defcolor, thick = lthick

oplot, xarr, -2.78 + xarr * 0.98, color = defcolor, linestyle = 1, thick = lthick
oplot, xarr, fitboth_ohm[0] + xarr*fitboth_ohm[1], color=defcolor, thick = lthick, linestyle=0
oplot, xarr, fitboth_con[0] + xarr*fitboth_con[1], color=defcolor, thick = lthick, linestyle=2

legend, /top, /left, $
	['OHMs', 'Non-masing','HK07'], $
	linestyle = [0,2,1], $
	color = [defcolor, defcolor, defcolor], $
	thick = lthick, charthick = cthick

if keyword_set(ps) then begin
;	!p.font = -1
	device,/close
	set_plot,'x'
endif

; Plot the mean star formation rate using both Ne and PAH methods

; Neon

nelumcgs_ohm = lsun * (10d)^neboth_ohm_lum
sfr_ohm_neon = nelumcgs_ohm * 4.73d-41

nelumcgs_con = lsun * (10d)^neboth_con_lum
sfr_con_neon = nelumcgs_con * 4.73d-41

; PAH (spline)

p6 =  [transpose(ohmdat('pah62lum')),transpose(archdat('pah62lum'))]
p11 = [transpose(ohmdat('pah11lum')),transpose(archdat('pah11lum'))]

match, taglist_ohm, ['mega007','mega034','arch023'], ta, tb
ind = setdifference(indgen(n_elements(p6)),ta)		
p6 = p6[ind] & p11 = p11[ind]

p6lum = lsun * (10d)^p6
p11lum = lsun * (10d)^p11
sfr_pah_ohm = (p6lum + p11lum) * 1.18d-41

p6_con =  [transpose(condat('pah62lum'))]
p11_con = [transpose(condat('pah11lum'))]

match, taglist_con, ['control008','control033'], ca, cb
ind_con = setdifference(indgen(n_elements(p6_con)),ca)		
p6_con = p6_con[ind_con] & p11_con = p11_con[ind_con]

p6lum_con = lsun * (10d)^p6_con
p11lum_con = lsun * (10d)^p11_con
sfr_pah_con = (p6lum_con + p11lum_con) * 1.18d-41

; PAH (PAHFIT)

pf6 =  [transpose(ohmdat('pahfit62lum')),transpose(archdat('pahfit62lum'))]
pf11 = [transpose(ohmdat('pahfit11lum')),transpose(archdat('pahfit11lum'))]

pf6 = pf6[ind] & pf11 = pf11[ind]

pf6lum = lsun * (10d)^pf6
pf11lum = lsun * (10d)^pf11
sfr_pahfit_ohm = (pf6lum + pf11lum) * 1.18d-41

pf6_con =  [transpose(condat('pahfit62lum'))]
pf11_con = [transpose(condat('pahfit11lum'))]

pf6_con = pf6_con[ind_con] & pf11_con = pf11_con[ind_con]

pf6lum_con = lsun * (10d)^pf6_con
pf11lum_con = lsun * (10d)^pf11_con
sfr_pahfit_con = (pf6lum_con + pf11lum_con) * 1.18d-41

print,''
print,'OHM SFR neon =       ',mean(sfr_ohm_neon), ' +_',stddev(sfr_ohm_neon)/sqrt(n_elements(sfr_ohm_neon))
print,'OHM SFR PAH spline = ',mean(sfr_pah_ohm), ' +_',stddev(sfr_pah_ohm)/sqrt(n_elements(sfr_pah_ohm))
print,'OHM SFR PAH PAHFIT = ',mean(sfr_pahfit_ohm), ' +_',stddev(sfr_pahfit_ohm)/sqrt(n_elements(sfr_pahfit_ohm))
print,'OHM NeII fit int:    ',fit2_ohm[0], ' +- ', err2_ohm[0]
print,'OHM NeII fit slope:  ',fit2_ohm[1], ' +- ', err2_ohm[1]
print,'OHM NeII + NeIII fit int:    ',fitboth_ohm[0], ' +- ', errboth_ohm[0]
print,'OHM NeII + NeIII fit slope:  ',fitboth_ohm[1], ' +- ', errboth_ohm[1]
print,''
print,'Control SFR neon =       ',mean(sfr_con_neon), ' +_',stddev(sfr_con_neon)/sqrt(n_elements(sfr_con_neon))
print,'Control SFR PAH spline = ',mean(sfr_pah_con), ' +_',stddev(sfr_pah_con)/sqrt(n_elements(sfr_pah_con))
print,'Control SFR PAH PAHFIT = ',mean(sfr_pahfit_con), ' +_',stddev(sfr_pahfit_con)/sqrt(n_elements(sfr_pahfit_con))
print,'Control NeII fit int:    ',fit2_con[0], ' +- ', err2_con[0]
print,'Control NeII fit slope:  ',fit2_con[1], ' +- ', err2_con[1]
print,'Control NeII + NeIII fit int:    ',fitboth_con[0], ' +- ', errboth_con[0]
print,'Control NeII + NeIII fit slope:  ',fitboth_con[1], ' +- ', errboth_con[1]
print,''

if keyword_set(stop) then stop
end
