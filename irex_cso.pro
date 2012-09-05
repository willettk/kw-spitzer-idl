pro irex_cso, ps = ps, label = label, stop = stop
;+
; NAME:
;       
;	IREX_CSO
;
; PURPOSE:
;
;	Plot the excitation level from IRS data using neon and sulfur data
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
;	IDL> irex_cso
;
; NOTES:
;
; 	Adapted from IREX.pro to use only the neon data
;
; REVISION HISTORY
;       Written by K. Willett                Jan 09
;	Added data from Cycle 5 CSOs
;	Corrected math on mean value vertical lines in plots - Aug 09
;	Removed NGC 5793, 1245+676 - Nov 09
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

	neII_hr_cso = neII_hr[where(strmid(neII_tag,0,3) eq 'cso' and strmid(neII_tag,0,6) ne 'cso003' and strmid(neII_tag,0,6) ne 'cso010')]
	neIII_hr_cso = neIII_hr[where(strmid(neIII_tag,0,3) eq 'cso' and strmid(neIII_tag,0,6) ne 'cso003' and strmid(neII_tag,0,6) ne 'cso010')]
	
	neII_tag_cso = neII_tag[where(strmid(neII_tag,0,3) eq 'cso' and strmid(neII_tag,0,6) ne 'cso003' and strmid(neII_tag,0,6) ne 'cso010')]
	neIII_tag_cso = neIII_tag[where(strmid(neIII_tag,0,3) eq 'cso' and strmid(neIII_tag,0,6) ne 'cso003' and strmid(neII_tag,0,6) ne 'cso010')]
	
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

plotname='~/Astronomy/Research/Spitzer/cso/papers/irex_cso.ps'
!p.multi=[0,1,4]

;erase

x0 = 0.15
x1 = 0.8
x2 = 0.95

y0 = 0.15
y1 = 0.8
y2 = 0.95

!x.style = 1
!y.style = 1

;plot, pah62ew_ohm, sil_ohm, $
;	/nodata, $
;	xtitle = '6.2 !7l!3m PAH EW [!7l!3m]', $
;	ytitle = 'Silicate strength', $
;	charsize = cs, $
;	thick = ls, $
;	xthick = ls, $
;	ythick = ls, $
;	charthick = ls, $
;	/xlog, /xstyle, $
;	ystyle = 1, $
;	xr = [1d-3, 3], $
;	yr = [1, -5], $
;	yticks = 6, yminor = 0, ytickv = reverse(fillarr(1,-5,1)), ytickname = ['1','0','-1','-2','-3','-4',' '], $
;	position = [x0,y0,x1,y1], $
;	color = defcolor


if keyword_set(ps) then begin
;	!p.font=0
	set_plot,'ps'
	device, filename = plotname, /color, /portrait, xs = 18, ys = 20, xoff = 1, yoff = 1
	cthick = 4
	lthick = 6
	labelsize = 1.5
	cs = 4
	arrowcolor=fsc_color("Grey")
	defcolor=fsc_color("Black")
endif else begin
	cthick = 1
	lthick = 1
	labelsize = 2
	arrowcolor=fsc_color("purple")
	defcolor=fsc_color("White")
endelse

bs = 0.4

; AGN

plothist, alog10(ne_agn), /halfbin, $
	position = [0.1, 0.73, 0.95, 0.94], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,20], $
	xticks = 5, xtickv = fillarr(0.5,-2.0,1.5), xtickname = replicate(' ',6), $
	yticks = 5, yminor = 0, ytickv = fillarr(5,0,15), ytickname = ['0','5','10','15'], $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

xyouts, 0.15, 0.87, 'AGN', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10(ne_agn)), linestyle = 2, thick = lthick

; ULIRGs - Farrah, Darling OHMs and control

plothist, alog10([ne_farrah[notlimits], ne_ohm, ne_arch, ne_control]), /halfbin, $
	position = [0.1, 0.52, 0.95, 0.73], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,40], $
	xticks = 5, xtickv = fillarr(0.5,-2.0,1.5), xtickname = replicate(' ',6), $
	yticks = 5, yminor = 0, ytickv = fillarr(10,0,40), ytickname = ['0','10','20','30','40'], $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

xyouts, 0.15, 0.66, 'ULIRGs', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10([ne_farrah[notlimits], ne_ohm, ne_arch, ne_control])), linestyle = 2, thick = lthick

; Starbursts

plothist, alog10(ne_starburst), /halfbin, $
	position = [0.1, 0.31, 0.95, 0.52], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,5], $
	xticks = 5, xtickv = fillarr(0.5,-2.0,1.5), xtickname = replicate(' ',6), $
	yticks = 6, yminor = 0, ytickv = fillarr(2,0,6), ytickname = ['0','2','4',' '], $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

xyouts, 0.15, 0.45, 'Starbursts', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10(ne_starburst)), linestyle = 2, thick = lthick

; CSOs

plothist, alog10(ne_cso), /halfbin, $
	position = [0.1, 0.1, 0.95, 0.31], $
	xr = [-2,1.5], /xstyle, $
	yr = [0,6], $
	yticks = 6, yminor = 0, ytickv = fillarr(2,0,6), ytickname = ['0','2','4',' '], $
	xtitle = 'log ([Ne III]/[Ne II])', $
;	ytitle = 'Count', $
	color=defcolor, $
	bin = bs, $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charsize = cs, $
	charthick = cthick

xyouts, 0.15, 0.24, 'CSOs', /normal, charsize = labelsize, charthick = cthick
ver, mean(alog10(ne_cso)), linestyle = 2, thick = lthick

xyouts, 0.0, 0.46, 'Count', /normal, orientation=90, charthick=cthick, charsize = 2

if keyword_set(ps) then begin
;	!p.font=-1
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
