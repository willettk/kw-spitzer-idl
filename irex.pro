pro irex, ps = ps, label = label, stop = stop, ohmstop = ohmstop
;+
; NAME:
;       
;	IREX
;
; PURPOSE:
;
;	Plot the excitation level from IRS data using neon and sulfur ratios
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
;	IDL> irex
;
; NOTES:
;
;	Follows Fig. 11 in Farrah et al. (2007)
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;-

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
	restore,'~/Astronomy/Research/Spitzer/linedata/sIII.sav'
	sIII_ohm = line
	restore,'~/Astronomy/Research/Spitzer/linedata/sIV.sav'
	sIV_ohm = line
	
	; Isolate targets where both neon lines are detected
	
	match, neii_ohm.tag, neiii_ohm.tag, a, b
	match, siii_ohm.tag, siv_ohm.tag, c, d
	
	neii_a = neii_ohm.tag[a] & neiii_b = neiii_ohm.tag[b]
	siii_c = siii_ohm.tag[c] & siv_d = siv_ohm.tag[d]
	
	match, neii_a, siii_c, aa, cc
	match, neiii_b, siv_d, bb, dd
	
	match, neii_ohm.tag, neii_a[aa], aaa, junk
	match, neiii_ohm.tag, neiii_b[bb], bbb, junk
	match, siii_ohm.tag, siii_c[cc], ccc, junk
	match, siv_ohm.tag, siv_d[dd], ddd, junk
	
	neII_hr = neii_ohm.flux[aaa]
	neIII_hr = neiii_ohm.flux[bbb]
	sIII_hr = siii_ohm.flux[ccc]
	sIV_hr = siv_ohm.flux[ddd]
	
	neII_tag = neii_ohm.tag[aaa]
	neIII_tag = neiii_ohm.tag[bbb]
	sIII_tag = siii_ohm.tag[ccc]
	sIV_tag = siv_ohm.tag[ddd]

	; Darling OHMs
	
	neII_hr_ohm = neII_hr[where(strmid(neII_tag,0,4) eq 'mega')]
	neIII_hr_ohm = neIII_hr[where(strmid(neIII_tag,0,4) eq 'mega')]
	sIII_hr_ohm = sIII_hr[where(strmid(sIII_tag,0,4) eq 'mega')]
	sIV_hr_ohm = sIV_hr[where(strmid(sIV_tag,0,4) eq 'mega')]
	
	neII_tag_ohm = neII_tag[where(strmid(neII_tag,0,4) eq 'mega')]
	neIII_tag_ohm = neIII_tag[where(strmid(neIII_tag,0,4) eq 'mega')]
	sIII_tag_ohm = sIII_tag[where(strmid(sIII_tag,0,4) eq 'mega')]
	sIV_tag_ohm = sIV_tag[where(strmid(sIV_tag,0,4) eq 'mega')]
	
	; Archived OHMs
	
	neII_hr_arch = neII_hr[where(strmid(neII_tag,0,4) eq 'arch')]
	neIII_hr_arch = neIII_hr[where(strmid(neIII_tag,0,4) eq 'arch')]
	sIII_hr_arch = sIII_hr[where(strmid(sIII_tag,0,4) eq 'arch')]
	sIV_hr_arch = sIV_hr[where(strmid(sIV_tag,0,4) eq 'arch')]
	
	neII_tag_arch = neII_tag[where(strmid(neII_tag,0,4) eq 'arch')]
	neIII_tag_arch = neIII_tag[where(strmid(neIII_tag,0,4) eq 'arch')]
	sIII_tag_arch = sIII_tag[where(strmid(sIII_tag,0,4) eq 'arch')]
	sIV_tag_arch = sIV_tag[where(strmid(sIV_tag,0,4) eq 'arch')]
	
	; Control sample
	
	neII_hr_control = neII_hr[where(strmid(neII_tag,0,7) eq 'control')]
	neIII_hr_control = neIII_hr[where(strmid(neIII_tag,0,7) eq 'control')]
	sIII_hr_control = sIII_hr[where(strmid(sIII_tag,0,7) eq 'control')]
	sIV_hr_control = sIV_hr[where(strmid(sIV_tag,0,7) eq 'control')]
	
	neII_tag_control = neII_tag[where(strmid(neII_tag,0,7) eq 'control')]
	neIII_tag_control = neIII_tag[where(strmid(neIII_tag,0,7) eq 'control')]
	sIII_tag_control = sIII_tag[where(strmid(sIII_tag,0,7) eq 'control')]
	sIV_tag_control = sIV_tag[where(strmid(sIV_tag,0,7) eq 'control')]
	
; Add data for AGN and starburst galaxies

; Sturm et al 2002 (AGN, ISO)

sturm_names = ['NGC 1068','NGC 1365','NGC7469']
sturm_neII = [70.0, 40.9, 22.6] * 1d-1			; Sturm measures fluxes in 10^-20 W/cm^2
sturm_neIII = [160.0, 7.7, 2.2] * 1d-1			; doesn't really matter, since they're ratios, but should be consistent
sturm_sIII = [40.0, 13.5, 9.2] * 1d-1
sturm_sIV = [58.0, 2.6, 0.9] * 1d-1

; Verma et al 2003 (starbursts, ISO)

verma_names = ['NGC 3256','NGC 3690A','NGC 3690B/C']
verma_neII = [89.2, 31.9, 27.7] * 1d-1			; Verma measures fluxes in 10^-20 W/cm^2
verma_neIII = [14.2, 9.3, 20.0] * 1d-1
verma_sIII = [32.5, 8.2, 18.4] * 1d-1
verma_sIV = [0.9, 1.0, 3.5] * 1d-1

; Tommasin et al 2008 (AGN, IRS)

tommasin_names = ['IRAS 00198-7926','ESO 541-IG012','Mrk 1034 NED02',$
	'IRAS F15091-2107','IRAS F22017+0319','NGC 7496']
tommasin_neII = [6.19, 1.87, 34.99, 11.52, 5.95, 48.08]
tommasin_neIII = [14.03, 2.02, 3.61, 16.29, 14.07, 6.67]
tommasin_sIII = [5.41, 1.65, 9.07, 9.95, 6.10, 23.48]			; Used SH measurement where available
tommasin_sIV = [8.10, 2.03, 1.89, 7.29, 10.31, 1.3]

; Slopes from Dale et al 2006

; SFR

x1 = 1.0
y1 = 4.2

; Seyfert

x2 = 1.0
y2 = 2.0

sfrslope = 0.75
seyslope = 0.75 ; originally 0.71

sfrb = alog10(y1) - sfrslope*alog10(x1)
seyb = alog10(y2) - seyslope*alog10(x2)

xarr = fillarr(1d-3,1d-2,10)

; Fit my hires OHM data

expr='p[0]+x*p[1]'
start = [sfrb,sfrslope]

fit = mpfitexpr(expr,alog10(sIV_hr/sIII_hr),alog10(neIII_hr/neII_hr),0.1*alog10(neIII_hr/neII_hr),start,/quiet)

; Plot data for paper I

cs = 1

plotname='~/Astronomy/Research/Spitzer/papers/excitation.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color, /portrait
	cthick = 3
	lthick = 3
	arrowcolor=fsc_color("Grey")
	defcolor=fsc_color("Black")
endif else begin
	cthick = 1
	lthick = 1
	arrowcolor=fsc_color("purple")
	defcolor=fsc_color("White")
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
purple = fsc_color("purple")
green = fsc_color("Forest Green")
orange = fsc_color("Orange")
purple = fsc_color("Purple")
grey = fsc_color("Dark Grey")

s_farrah = sIV / sIII
ne_farrah = neIII / neII

s_ohm = sIV_hr_ohm / sIII_hr_ohm
ne_ohm = neIII_hr_ohm / neII_hr_ohm

s_arch = sIV_hr_arch / sIII_hr_arch
ne_arch = neIII_hr_arch / neII_hr_arch

s_control = sIV_hr_control / sIII_hr_control
ne_control = neIII_hr_control / neII_hr_control

s_agn = [tommasin_sIV, sturm_sIV] / [tommasin_sIII, sturm_sIII]
ne_agn = [tommasin_neIII, sturm_neIII] / [tommasin_neII, sturm_neII]

s_tom = [tommasin_sIV] / [tommasin_sIII]
ne_tom = [tommasin_neIII] / [tommasin_neII]

s_stu = [sturm_sIV] / [sturm_sIII]
ne_stu = [sturm_neIII] / [sturm_neII]

s_starburst = verma_sIV / verma_sIII
ne_starburst = verma_neIII / verma_neII

!p.multi=[0,1,1]
plot, sIV/sIII, neIII/neII, $
	/nodata, $
	xtitle = '[SIV]/[SIII]', $
	ytitle = '[NeIII]/[NeII]', $
	xr=[0.015,2.5], /xlog, /xstyle, $
	yr=[0.075,4], /ylog, /ystyle, $
	charsize = cs, $
	charthick = cthick, $
	thick = lthick

; Farrah

oplot, s_farrah(notlimits), ne_farrah(notlimits), psym=symcat(16), color=defcolor
plotsym,6,2,thick=lthick
;oplot, s_farrah[slim], ne_farrah[slim], psym=8,color=defcolor

; Darling

oplot,s_ohm, ne_ohm, psym=symcat(14), color=red

; AGN

oplot, s_agn, ne_agn, psym=symcat(14), color=orange

; Starbursts

oplot, s_starburst, ne_starburst, psym=symcat(14), color=green

oplot,xarr,10^(sfrslope*alog10(xarr)+sfrb),linestyle=1,thick=lthick
oplot,xarr,10^(seyslope*alog10(xarr)+seyb),linestyle=2,thick=lthick
;oplot,xarr,10^(fit(1)*alog10(xarr)+fit(0)),linestyle=0,thick=lthick	; Likely too few points to plot a meaningful fit

legend, /top,/left, ['Star-forming','Seyfert'], linestyle=[1,2], thick=lthick, charthick=cthick, charsize = cs
legend, /bottom,/right, ['OHM','AGN','Starburst','ULIRG'], psym=[14,14,14,16], thick=lthick, charthick=cthick, $
	color=[red, orange, green, defcolor], charsize = cs

; Labels

if keyword_set(label) then begin
	xyouts,s_ohm * 1.1, ne_ohm, neII_tag_ohm, color=red
	xyouts,s_agn * 1.1, ne_agn, [tommasin_names, sturm_names], color=orange
	xyouts,s_starburst * 1.1, ne_starburst, verma_names, color=green
	xyouts,s_farrah[notlimits] * 1.1, ne_farrah[notlimits], gal[notlimits], color=defcolor
endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Plot data for paper II

plotname='~/Astronomy/Research/Spitzer/papers/excitationII.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color, /portrait
	cthick = 5
	lthick = 5
	cs=1.5
	arrowcolor=fsc_color("Grey")
	defcolor=fsc_color("Black")
endif else begin
	cthick = 1
	lthick = 1
	cs=1
	arrowcolor=fsc_color("purple")
	defcolor=fsc_color("White")
endelse

plot, sIV/sIII, neIII/neII, $
	/nodata, $
	xtitle = '[SIV]/[SIII]', $
	ytitle = '[NeIII]/[NeII]', $
	xr=[0.01,3], /xlog, /xstyle, $
	yr=[0.05,6], /ylog, /ystyle, $
	charsize = cs, $
	charthick = cthick, $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick

; Darling and archived OHMs

oplot, s_ohm, ne_ohm, psym=symcat(14), color=red
oplot, s_arch, ne_arch, psym=symcat(14), color=red

; Control sample

oplot, s_control, ne_control, psym=symcat(15), color=blue

; AGN

oplot, s_agn, ne_agn, psym=symcat(18), color=grey
;oplot, s_tom, ne_tom, psym=symcat(16), color=fsc_color("Hot Pink")
;oplot, s_stu, ne_stu, psym=symcat(16), color=fsc_color("Yellow")

; Starbursts

oplot, s_starburst, ne_starburst, psym=symcat(17), color=grey

; Farrah ULIRGs

oplot, s_farrah(notlimits), ne_farrah(notlimits), psym=symcat(16), color=grey
;plotsym,6,2,thick=lthick
;oplot, s_farrah[slim], ne_farrah[slim], psym=8,color=defcolor

; Labels

if keyword_set(label) then begin
	xyouts,s_ohm * 1.1, ne_ohm, neII_tag_ohm, color=red
	xyouts,s_arch * 1.1, ne_arch, neII_tag_arch, color=red
	xyouts,s_control * 1.1, ne_control, neII_tag_control, color=blue
	xyouts,s_agn * 1.1, ne_agn, [tommasin_names, sturm_names], color=grey
	xyouts,s_starburst * 1.1, ne_starburst, verma_names, color=grey
	xyouts,s_farrah[notlimits] * 1.1, ne_farrah[notlimits], gal[notlimits], color=grey
endif

legend, /top, /left, ['OHM','non-masing','AGN','starburst','ULIRG'], psym=[14,15,18,17,16], thick=lthick, charthick=cthick, $
	color=[red, blue, grey, grey, grey], charsize = cs

; Fits from Dale et al. 

;oplot,xarr,10^(sfrslope*alog10(xarr)+sfrb),linestyle=1,thick=lthick
;oplot,xarr,10^(seyslope*alog10(xarr)+seyb),linestyle=2,thick=lthick
;;oplot,xarr,10^(fit(1)*alog10(xarr)+fit(0)),linestyle=0,thick=lthick	; Likely too few points to plot a meaningful fit
;legend, /top,/left, ['Star-forming','Seyfert'], linestyle=[1,2], thick=lthick, charthick=cthick, charsize = cs

fit_ohm = mpfitexpr(expr,alog10([s_ohm,s_arch]),alog10([ne_ohm,ne_arch]),alog10(0.1 * [ne_ohm,ne_arch]),start,/quiet, perr=perr_ohm)
fit_con = mpfitexpr(expr,alog10(s_control),alog10(ne_control),alog10(0.1 * ne_control),start,/quiet, perr=perr_con)
fit_ulirg = mpfitexpr(expr,alog10(s_farrah),alog10(ne_farrah),alog10(0.1 * ne_farrah),start,/quiet, perr=perr_ulirg)
fit_starburst = mpfitexpr(expr,alog10(s_starburst),alog10(ne_starburst),alog10(0.1 * ne_starburst),start,/quiet, perr=perr_starburst)
fit_agn = mpfitexpr(expr,alog10(s_agn),alog10(ne_agn),alog10(0.1 * ne_agn),start,/quiet, perr=perr_agn)

;oplot, xarr, 10^(fit_ohm[1]*alog10(xarr)+fit_ohm[0]), color=red, linestyle=2, thick=lthick
;oplot, xarr, 10^(fit_con[1]*alog10(xarr)+fit_con[0]), color=blue, linestyle=2, thick=lthick
;oplot, xarr, 10^(fit_ulirg[1]*alog10(xarr)+fit_ulirg[0]), color=defcolor, linestyle=2, thick=lthick
;oplot, xarr, 10^(fit_starburst[1]*alog10(xarr)+fit_starburst[0]), color=green, linestyle=2, thick=lthick
;oplot, xarr, 10^(fit_agn[1]*alog10(xarr)+fit_agn[0]), color=orange, linestyle=2, thick=lthick

print,''
print, 'OHMs: ',fit_ohm, perr_ohm
print, 'Non-masing: ',fit_con, perr_con
print, 'AGN: ',fit_ulirg, perr_ulirg
print, 'Starbursts: ',fit_starburst, perr_starburst
print, 'ULIRG: ',fit_ulirg, perr_ulirg
print,''


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(ohmstop) then stop

; Plot data for CSO paper

; Load in CSO data

	restore,'~/Astronomy/Research/Spitzer/cso/linedata/neII.sav'
	neII_cso = line
	restore,'~/Astronomy/Research/Spitzer/cso/linedata/neIII.sav'
	neIII_cso = line
	restore,'~/Astronomy/Research/Spitzer/cso/linedata/sIII.sav'
	sIII_cso = line
	restore,'~/Astronomy/Research/Spitzer/cso/linedata/sIV.sav'
	sIV_cso = line
	
	; Isolate targets where both neon lines are detected
	
	match, neii_cso.tag, neiii_cso.tag, a, b
	match, siii_cso.tag, siv_cso.tag, c, d
	
	neii_a = neii_cso.tag[a] & neiii_b = neiii_cso.tag[b]
	siii_c = siii_cso.tag[c] & siv_d = siv_cso.tag[d]
	
	match, neii_a, siii_c, aa, cc
	match, neiii_b, siv_d, bb, dd
	
	match, neii_cso.tag, neii_a[aa], aaa, junk
	match, neiii_cso.tag, neiii_b[bb], bbb, junk
	match, siii_cso.tag, siii_c[cc], ccc, junk
	match, siv_cso.tag, siv_d[dd], ddd, junk
	
	neII_hr = neii_cso.flux[aaa]
	neIII_hr = neiii_cso.flux[bbb]
	sIII_hr = siii_cso.flux[ccc]
	sIV_hr = siv_cso.flux[ddd]
	
	neII_tag = neii_cso.tag[aaa]
	neIII_tag = neiii_cso.tag[bbb]
	sIII_tag = siii_cso.tag[ccc]
	sIV_tag = siv_cso.tag[ddd]

	neII_hr_cso = neII_hr[where(strmid(neII_tag,0,3) eq 'cso')]
	neIII_hr_cso = neIII_hr[where(strmid(neIII_tag,0,3) eq 'cso')]
	sIII_hr_cso = sIII_hr[where(strmid(sIII_tag,0,3) eq 'cso')]
	sIV_hr_cso = sIV_hr[where(strmid(sIV_tag,0,3) eq 'cso')]
	
	neII_tag_cso = neII_tag[where(strmid(neII_tag,0,3) eq 'cso')]
	neIII_tag_cso = neIII_tag[where(strmid(neIII_tag,0,3) eq 'cso')]
	sIII_tag_cso = sIII_tag[where(strmid(sIII_tag,0,3) eq 'cso')]
	sIV_tag_cso = sIV_tag[where(strmid(sIV_tag,0,3) eq 'cso')]
	
; For the ancillary data, now include ALL objects with line ratios (CSOs span a large range in L_IR)

; Sturm et al 2002 (AGN, ISO)

sturm_names = ['Cen A','Circinus','NGC 1068','NGC 1365','NGC 4151','NGC 5506','NGC 7469','NGC 7582']
sturm_neII = [22.1,90.0,70.0,40.9,11.8,5.9,22.6,14.8] * 1d-1	; Sturm measures fluxes in 10^-20 W/cm^2
sturm_neIII = [14.1,33.5,160.0,7.7,20.7,5.8,2.2,6.7] * 1d-1	
sturm_sIII = [6.4,35.2,40.0,13.5,5.4,3.0,9.2,5.2] * 1d-1
sturm_sIV = [1.4,12.7,58.0,2.6,11.3,5.4,0.9,1.8] * 1d-1

; Verma et al 2003 (starbursts, ISO)

verma_names = ['II Zw 40','M 82','NGC 3256','NGC 3690A','NGC 3690B/C','NGC 4038','NGC 5236','NGC 5253','NGC 7552']
verma_neII = [1.7,714.0,89.2,31.9,27.7,7.7,133.9,7.9,68.0] * 1d-1		; Verma measures fluxes in 10^-20 W/cm^2
verma_neIII = [17.0,126.0,14.2,9.3,20.0,6.5,6.8,28.9,5.1] * 1d-1
verma_sIII = [5.1,252.0,32.5,8.2,18.4,7.3,54.4,15.9,24.6] * 1d-1
verma_sIV = [19.8,14.9,0.9,1.0,3.5,2.7,0.9,34.0,0.3] * 1d-1

; Tommasin et al 2008 (AGN, IRS)

tommasin_names = ['IRAS 00198-7926','ESO 012-G021','ESO 541-IG012','NGC 424','NGC 513','Mrk 1034 NED02',$
	'ESO 545-G013','NGC 931','NGC 1125','ESO 033-G002','Mrk 6','Mrk 9','NGC 3516','NGC 3660',$
	'TOLOLO 1238-364','NGC 4748','NGC 4968','Mrk 817','IRAS F15091-2107','ESO 141-G055','NGC 6890',$
	'IRAS F22017+0319','NGC 7496']
tommasin_neII = [6.19,11.95,1.87,8.70,12.76,34.99,10.05,5.47,16.37,2.13,28.00,3.23,8.07,6.51,$
	45.15,7.37,24.90,3.83,11.52,2.24,11.32,5.95,48.08]
tommasin_neIII = [14.03,6.42,2.02,18.45,4.43,3.61,10.50,15.41,15.55,9.22,49.34,1.90,17.72,1.49,27.00,$
	15.93,33.80,4.58,16.29,5.62,6.57,14.07,6.67]
tommasin_sIII = [5.41,5.63,1.65,6.96,6.76,9.07,3.59,4.86,10.99,3.99,14.10,2.38,5.86,3.66,16.32,$
	7.61,15.10,2.76,9.95,1.75,4.34,6.10,23.48]			; Used SH measurement where available
tommasin_sIV = [8.10,2.5,2.03,8.98,2.77,1.89,5.32,10.70,6.07,5.54,16.69,2.37,13.33,1.48,5.70,$
	9.87,9.63,1.53,7.29,3.45,2.92,10.31,1.3]

s_cso = sIV_hr_cso / sIII_hr_cso
ne_cso = neIII_hr_cso / neII_hr_cso

s_agn = [tommasin_sIV, sturm_sIV] / [tommasin_sIII, sturm_sIII]
ne_agn = [tommasin_neIII, sturm_neIII] / [tommasin_neII, sturm_neII]

s_starburst = verma_sIV / verma_sIII
ne_starburst = verma_neIII / verma_neII

; Add limits for the three CSOs without SIV detections

csolim = ['cso001','cso002','cso004']
cso_neII_lim = dblarr(n_elements(csolim))
cso_neIII_lim = dblarr(n_elements(csolim))
cso_sIII_lim = dblarr(n_elements(csolim))
cso_sIV_lim = dblarr(n_elements(csolim))

for i = 0, n_elements(csolim)-1 do begin
	junk1 = getlineflux(csolim[i],'sIII')
	cso_sIII_lim[i] = junk1[0]
	junk2 = getlineflux(csolim[i],'neIII')
	cso_neIII_lim[i] = junk2[0]
	junk3 = getlineflux(csolim[i],'neII')
	cso_neII_lim[i] = junk3[0]
	if csolim[i] eq 'cso002' then width = 0.05 else width = 0.2
	cso_sIV_lim[i] = linelim(csolim[i],'sIV',/noplot,width = width) * 1d-21
endfor

s_csolim = cso_sIV_lim / cso_sIII_lim
ne_csolim = cso_neIII_lim / cso_neII_lim

plotname='~/Astronomy/Research/Spitzer/cso/papers/excitation_cso.ps'

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color, /portrait
	cthick = 4
	lthick = 4
	arrowcolor=fsc_color("Grey")
	defcolor=fsc_color("Black")
endif else begin
	cthick = 1
	lthick = 1
	arrowcolor=fsc_color("purple")
	defcolor=fsc_color("White")
endelse

plot, sIV/sIII, neIII/neII, $
	/nodata, $
	xtitle = '[SIV]/[SIII]', $
	ytitle = '[NeIII]/[NeII]', $
	xr=[0.015,5], /xlog, /xstyle, $
	yr=[0.075,21], /ylog, /ystyle, $
	;xr=[0.015,2.5], /xlog, /xstyle, $
	;yr=[0.075,5], /ylog, /ystyle, $
	charsize = cs, $
	charthick = cthick, $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick

; Farrah

oplot, s_farrah(notlimits), ne_farrah(notlimits), psym=symcat(16), color=defcolor

; Darling and archived OHMs

oplot, s_ohm, ne_ohm, psym=symcat(16), color=defcolor
oplot, s_arch, ne_arch, psym=symcat(16), color=defcolor

; Control sample

oplot, s_control, ne_control, psym=symcat(16), color=defcolor

; AGN

oplot, s_agn, ne_agn, psym=symcat(15), color=orange

; Starbursts

oplot, s_starburst, ne_starburst, psym=symcat(14), color=green

; CSOs

oplot, s_cso, ne_cso, psym = symcat(18), color = red, symsize = 1.5
oplot, s_csolim, ne_csolim, psym = symcat(18), color=red, symsize = 1.5
arrow, s_csolim, ne_csolim, s_csolim * 0.8, ne_csolim, color=red,/data, thick = lthick

; Labels

if keyword_set(label) then begin
	xyouts,s_ohm * 1.1, ne_ohm, neII_tag_ohm, color=defcolor
	xyouts,s_arch * 1.1, ne_arch, neII_tag_arch, color=defcolor
	xyouts,s_control * 1.1, ne_control, neII_tag_control, color=defcolor
	xyouts,s_agn * 1.1, ne_agn, [tommasin_names, sturm_names], color=orange
	xyouts,s_starburst * 1.1, ne_starburst, verma_names, color=green
	xyouts,s_cso * 1.1, ne_cso, neII_tag_cso, color=purple
	xyouts,s_csolim * 1.1, ne_csolim, csolim, color=purple
endif

legend, /bottom,/right, $
	['AGN','Starburst','ULIRG','CSO'], $
	psym=[15,14,16,18], $
	color=[orange, green, defcolor, red], $
	thick=lthick, charthick=cthick, charsize = cs

; Fits from Dale et al. 

oplot,xarr,10^(sfrslope*alog10(xarr)+sfrb),linestyle=1,thick=lthick
oplot,xarr,10^(seyslope*alog10(xarr)+seyb),linestyle=2,thick=lthick
;oplot,xarr,10^(fit(1)*alog10(xarr)+fit(0)),linestyle=0,thick=lthick	; Likely too few points to plot a meaningful fit
legend, /top,/left, ['Star-forming','Seyfert'], linestyle=[1,2], thick=lthick, charthick=cthick, charsize = cs

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop
end
