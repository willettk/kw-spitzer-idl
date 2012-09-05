;pro ncp2, ps=ps

;+
; NAME:
;       
;	
;
; PURPOSE:
;
;	Select target sample of OHMs and non-masing galaxies from OH surveys; make plots
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
;	IDL> .r ncp2
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                
; 	Plot limits for max L_OH for new objects w/Spitzer spectra - KW, 26 Feb 08
;	Added header - Mar 09
;	Black and white figures for paper - Feb 10
;	Eliminate the five OHMs with L_OH < 10^2.3 L_sun - Jun 10
;-

ps = 1

; Physical constants

c = 299792.458d ; km/s
h0 = 75 ; km/s/Mpc
mpc2cm = 3.09d24 ; cm/Mpc
mjy2cgs = 1d-26 ; erg s^-1 cm^-2 Hz^-1 / mJy
lsun = 3.826d33 ; erg/s

; Reading in distance, OH rms data

; Kent et al 2003

kfile = '~/Astronomy/Research/Spitzer/control/kent2.lst'
readcol,kfile,kname,kz,krms,format='a,f,f',/silent
kvel = c * kz
k_dl = (kvel / h0) * (1d + 0.5d * kz)
klohmax = (4 * !dpi * (k_dl * mpc2cm)^2 * 1.5 * (krms * mjy2cgs) * (150/c) * (1665d6 / (1+kz))) / lsun

; Staveley-Smith et al 1992

ssfile = '~/Astronomy/Research/Spitzer/control/ss.lst'
readcol,ssfile,ssname,ssvel,ssrms,format='a,f,f',/silent
ssz = ssvel / c 
ss_dl = (ssvel / h0) * (1d + 0.5d * ssz)
sslohmax = (4 * !dpi * (ss_dl * mpc2cm)^2 * 1.5 * (ssrms * mjy2cgs) * (150/c) * (1665d6 / (1+ssz))) / lsun

; Darling & Giovanelli I,II,III
dfile = '~/Astronomy/Research/Spitzer/control/darling.lst'
readcol,dfile,dname,dz,drms,format='a,f,f',/silent
dname = 'IRAS'+dname
dvel = c * dz
d_dl = (dvel / h0) * (1d + 0.5d * dz)
dlohmax = (4 * !dpi * (d_dl * mpc2cm)^2 * 1.5 * (drms * mjy2cgs) * (150/c) * (1665d6 / (1+dz))) / lsun

; Baan et al 1992

bfile = '~/Astronomy/Research/Spitzer/control/baan.lst'
readcol,bfile,bname,bvel,brms,format='a,f,f',/silent
bz = bvel / c 
b_dl = (bvel / h0) * (1d + 0.5d * bz)
blohmax = (4 * !dpi * (b_dl * mpc2cm)^2 * 1.5 * (brms * mjy2cgs) * (150/c) * (1665d6 / (1+bz))) / lsun

; OHM data

oh_reg = float(transpose(ohmdat('logoh')))
oh_arch = float(transpose(archdat('logoh')))
ohz = float(transpose(ohmdat('redshift')))
archz = float(transpose(archdat('redshift')))
aname = transpose(archdat('obj'))

oh_iras = ohmdat('iras')
oh_vel = c * ohz
oh_dl = (oh_vel / h0) * (1d + 0.5d * ohz)
oh_lfir = lfir(oh_iras(2,*),oh_iras(3,*),oh_dl)

a_vel = c * archz
a_dl = (a_vel / h0) * (1d + 0.5d * archz)

mysym = symcat(16)
green=fsc_color("Green")
orange=fsc_color("Orange")
red=fsc_color("Red")
blue=fsc_color("Blue")

; Histogram the data

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/control/plots/ncp2.ps',/color,/landscape
;		/inches,xsize=7,ysize=10,xoffset=0.5,yoffset=0.5
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")


!p.multi=[0,1,1]
!p.charthick=2
 
; OHMs
plothist,[bz,dz,ssz,kz],bin=0.02,/fill,forient=-45,fcolor=blue, /fline, $
	xrange=[0,0.3],yrange=[0,55],/ystyle,$
	xtitle='z',ytitle='Frequency'
plothist,[ohz,archz],/fill,forient=45,fcolor=red,bin=0.02,/overplot, /fline
legend, /top,/right,['Control sample','OHMs'],psym=[replicate(mysym,2)], $
	color=[blue,red]
xyouts,0.2,30,'N!IOHM!N = '+string(n_elements([ohz,archz]),format='(i3)'),charsize=2.5,/data
xyouts,0.2,20,'N!Icon!N = '+string(n_elements([bz,dz,ssz,kz]),format='(i3)'),charsize=2.5,/data

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Read in the IRAS data for each object as obtained from NED

; Baan

bfile_ir='~/Astronomy/Research/Spitzer/control/ned_baan.lst'
openr,lun,bfile_ir,/get_lun
nbaan = numlines(bfile_ir)  & nbsource = nbaan/5
if nbaan mod 5 ne 0 then print,'Wrong number of files in Baan list'
baandat = strarr(nbaan)
readf,lun,baandat
close,lun

bname_ir = baandat(indgen(nbsource)*5)
bname_ir = strmid(bname_ir,10,100)

for i = 0, nbsource - 1 do begin
	temp2 = strsplit(bname_ir(i),/extract)
	if temp2(0) eq 'IRAS' then temp2(0) = 'IR'
	bname_ir(i) = strupcase(strjoin(temp2))
endfor

b12 = baandat(indgen(nbsource)*5 + 1)
b12 = strmid(b12,27,11)
b25 = baandat(indgen(nbsource)*5 + 2)
b25 = strmid(b25,27,11)
b60 = baandat(indgen(nbsource)*5 + 3)
b60 = strmid(b60,27,11)
b100 = baandat(indgen(nbsource)*5 + 4)
b100 = strmid(b100,27,11)

j=0
bmatch = intarr(n_elements(bname))
for i = 0, n_elements(bname)-1 do begin
	temp1 = where(bname(i) eq bname_ir)
	if temp1 ne -1 then begin
		bmatch(j) = temp1
		j=j+1
	endif
endfor

bmatch = bmatch(0:j-1)
bnames_both = bname_ir(bmatch)
bdl_both = fltarr(n_elements(bnames_both))
bz_both = fltarr(n_elements(bnames_both))
blohmax_both = fltarr(n_elements(bnames_both))

b12_both = fltarr(n_elements(bnames_both))
b25_both = fltarr(n_elements(bnames_both))
b60_both = fltarr(n_elements(bnames_both))
b100_both = fltarr(n_elements(bnames_both))

for i = 0, n_elements(bnames_both) - 1 do begin
	temp3 = where(bnames_both(i) eq bname)
	bdl_both(i) = b_dl(temp3)
	bz_both(i) = bz(temp3)
	blohmax_both(i) = blohmax(temp3)

	temp_irname = where(bnames_both(i) eq bname_ir)
	b12_both(i) = b12(temp_irname)
	b25_both(i) = b25(temp_irname)
	b60_both(i) = b60(temp_irname)
	b100_both(i) = b100(temp_irname)
endfor

b_lir = lir(b12_both,b25_both,b60_both,b100_both,bdl_both)
b_lfir = lfir(b60_both,b100_both,bdl_both)

; Darling

dfile_ir='~/Astronomy/Research/Spitzer/control/ned_darling.lst'
openr,lun,dfile_ir,/get_lun
ndarling = numlines(dfile_ir)  & ndsource = ndarling/5
if ndarling mod 5 ne 0 then print,'Wrong number of files in Darling list'
darlingdat = strarr(ndarling)
readf,lun,darlingdat
close,lun

dname_ir = darlingdat(indgen(ndsource)*5)
dname_ir = strmid(dname_ir,10,100)

for i = 0, ndsource - 1 do begin
	temp2 = strsplit(dname_ir(i),/extract)
;	if temp2(0) eq 'IRAS' then temp2(0) = ''
	dname_ir(i) = strupcase(strjoin(temp2))
endfor

d12 = darlingdat(indgen(ndsource)*5 + 1)
d12 = strmid(d12,27,11)
d25 = darlingdat(indgen(ndsource)*5 + 2)
d25 = strmid(d25,27,11)
d60 = darlingdat(indgen(ndsource)*5 + 3)
d60 = strmid(d60,27,11)
d100 = darlingdat(indgen(ndsource)*5 + 4)
d100 = strmid(d100,27,11)

j=0
dmatch = intarr(n_elements(dname))
for i = 0, n_elements(dname)-1 do begin
	temp1 = where(dname(i) eq dname_ir)
	if temp1 ne -1 then begin
		dmatch(j) = temp1
		j=j+1
	endif
endfor

dmatch = dmatch(0:j-1)
dnames_both = dname_ir(dmatch)
ddl_both = fltarr(n_elements(dnames_both))
dz_both = fltarr(n_elements(dnames_both))
dlohmax_both = fltarr(n_elements(dnames_both))

for i = 0, n_elements(dnames_both) - 1 do begin
	temp3 = where(dnames_both(i) eq dname)
	ddl_both(i) = d_dl(temp3)
	dz_both(i) = dz(temp3)
	dlohmax_both(i) = dlohmax(temp3)
endfor

d_lir = lir(d12,d25,d60,d100,ddl_both)
d_lfir = lfir(d60,d100,ddl_both)

; Archived

afile_ir='~/Astronomy/Research/Spitzer/archived/ned_arch.lst'
openr,lun,afile_ir,/get_lun
narch = numlines(afile_ir)  & nasource = narch/5
if narch mod 5 ne 0 then print,'Wrong number of files in Darling list'
archir = strarr(narch)
readf,lun,archir
close,lun

aname_ir = archir(indgen(nasource)*5)
aname_ir = strmid(aname_ir,10,100)

for i = 0, nasource - 1 do begin
	temp2 = strsplit(aname_ir(i),/extract)
	if temp2(0) eq 'IRAS' then temp2(0) = 'IRAS '
	aname_ir(i) = strupcase(strjoin(temp2))
endfor

a12 = archir(indgen(nasource)*5 + 1)
a12 = strmid(a12,27,11)
a25 = archir(indgen(nasource)*5 + 2)
a25 = strmid(a25,27,11)
a60 = archir(indgen(nasource)*5 + 3)
a60 = strmid(a60,27,11)
a100 = archir(indgen(nasource)*5 + 4)
a100 = strmid(a100,27,11)

a_lir = lir(a12,a25,a60,a100,a_dl)
a_lfir = lfir(a60,a100,a_dl)

goodz = archz
good_lfir = a_lfir
goodnames = aname_ir

; Plot the data

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/control/plots/ncp2_ir.ps',/landscape,/color
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

mysym = symcat(16)
ohmsym = symcat(15)
red=fsc_color("Red")
green=fsc_color("Green")
blue=fsc_color("Blue")
orange=fsc_color("Orange")
yellow=fsc_color("Yellow")
black=fsc_color("Black")

!p.multi=0

plot, bz_both, b_lfir, $
	xrange=[0,0.35], $
	yrange=[8,14], $
	xtitle='z', $
	ytitle='L!IFIR!N!E(max)!N', $
	charsize=2, $
	psym=mysym, $
	title='FIR luminosity of OHMs and non-detections', $
	/nodata


oplot,bz_both,b_lfir,psym=mysym,color=blue
oplot,dz_both,d_lfir,psym=mysym,color=defcolor
oplot,ohz,oh_lfir,psym=ohmsym,color=red
oplot,archz,a_lfir,psym=ohmsym,color=red
lowbaan=where(bz_both lt 0.02 and b_lfir gt 11)
goodbaan=where(bz_both gt 0.02 and b_lfir gt 11)
gbaan=where(bz_both gt 0.02 and bz_both lt 0.04 and b_lfir gt 11)
gbaan2=where(bz_both gt 0.04 and b_lfir gt 11)

legend, /bottom,/right,['Baan control','Darling control','OHMs'],psym=[replicate(mysym,3)], $
	color=[blue,defcolor,red]

ver,0.1,linestyle=2
hor,11,linestyle=2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

allgoodbaanind=[lowbaan(5),gbaan([0,1,2,4]),gbaan2]

; Histogram again

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/control/plots/ncp2_hist2.ps',/landscape,/color
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

plothist, ohz, $
	xtitle='z', $
	ytitle='Frequency', $
	yrange=[0,12], /ystyle, $
	xrange=[0,0.35], /xstyle, $
	title='Redshift distribution of OHMs and control objects', $
	/nodata

plothist,[ohz,goodz], bin=0.02, color=red, linethick=2, /overplot
plothist,[bz_both(allgoodbaanind),dz_both], bin=0.02, color=blue, linethick=2, /overplot

finalnames = [bnames_both(allgoodbaanind),dnames_both]
legend,/top,/right,['OHMs','Control'],color=[red,blue],psym=[replicate(mysym,2)]

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Final z-L_FIR plot

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/control/plots/ncp2_lfir_final.ps',/landscape,/color
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

plot, bz_both, b_lfir, $
	xrange=[0,0.35], $
	yrange=[8,14], $
	xtitle='z', $
	ytitle='L!IFIR!N!E(max)!N', $
	charsize=2, $
	psym=mysym, $
	title='Final z-L!IFIR!N dist. of Spitzer OHMs and control', $
	/nodata


oplot,bz_both(allgoodbaanind),b_lfir(allgoodbaanind),psym=symcat(14),color=blue
oplot,dz_both,d_lfir,psym=symcat(15),color=blue
oplot,ohz,oh_lfir,psym=symcat(16),color=red
oplot,goodz,good_lfir,psym=symcat(17),color=red

legend, /bottom,/right,['Baan','Darling','OHM (Darling)', 'OHM (archive)'],psym=[14,15,16,17], $
	color=[blue,blue,red,red]

ver,0.1,linestyle=2
hor,11,linestyle=2


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Final z-L_OH^max plot

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/control/plots/ncp2_lohmax_final.ps',/landscape,/color
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

plot, bz_both, blohmax_both, $
	xrange=[0,0.3], /xstyle, $
	yrange=[-1,4], /ystyle, $
	xtitle='z', $
	ytitle='log L!IOH!N!E(max)!N', $
	charsize=2, $
	psym=mysym, $
	title='Final z-L!IOH!N!E(max)!N distribution', $
	/nodata


oplot,bz_both(allgoodbaanind),alog10(blohmax_both(allgoodbaanind)),psym=symcat(14),color=blue
oplot,dz_both,alog10(dlohmax_both),psym=symcat(15),color=blue
oplot,goodz,oh_arch,psym=symcat(16),color=red
arrow,bz_both(allgoodbaanind),alog10(blohmax_both(allgoodbaanind)),bz_both(allgoodbaanind),alog10(blohmax_both(allgoodbaanind))-0.2, $
	color=blue, /data
arrow,dz_both,alog10(dlohmax_both),dz_both,alog10(dlohmax_both)-0.2,color=blue,/data

legend, /bottom,/right,['Baan','Darling','Archived OHMs'],psym=[14,15,16], $
	color=[blue,blue,red]

ver,0.1,linestyle=2
hor,2,linestyle=2


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; z-L_OH plot for paper, including OHMs from our original sample

if keyword_set(ps) then begin
	set_plot,'ps'
;	device,filename='~/Astronomy/thesis/figures/ncutplot_paper.ps',/color,/portrait
	device,filename='~/Astronomy/Research/Spitzer/papers/ncutplot_paper.ps',/portrait
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

ohmcolor = defcolor
concolor = defcolor
;ohmcolor = red
;concolor = blue

plot, bz_both, blohmax_both, $
	xrange=[0,0.25], /xstyle, $
	yrange=[0,4], /ystyle, $
	xtitle='redshift', $
	ytitle='log L!IOH!N!E!N [L'+sunsymbol()+']', $
	charsize=1.5, $
	psym=mysym, $
	charthick=4, $
	thick=4, $
	xthick=4,$
	ythick=4,$
	/nodata

; Control sample

oplot,bz_both(allgoodbaanind),alog10(blohmax_both(allgoodbaanind)),psym=symcat(7),color=concolor, thick=4
oplot,dz_both,alog10(dlohmax_both),psym=symcat(7),color=concolor, thick=4
arrow,bz_both(allgoodbaanind),alog10(blohmax_both(allgoodbaanind)),$
	bz_both(allgoodbaanind),alog10(blohmax_both(allgoodbaanind))-0.3, $
	color=concolor, /data, thick=3, hsize=300
arrow,dz_both,alog10(dlohmax_both),dz_both,alog10(dlohmax_both)-0.3,color=concolor,/data,thick=3, hsize=300

; OHMs

; Eliminate OHMs with L_OH < 10^2.3 L_sun. This means the non-masing and OHM samples will be completely separate - KW, Jun 2010

arch_newlim = where(oh_arch gt 2.3)
oh_newlim   = where(oh_reg  gt 2.3)

oplot, archz[arch_newlim], oh_arch[arch_newlim], psym=symcat(14), color=ohmcolor
oplot, ohz[oh_newlim],     oh_reg[oh_newlim],    psym=symcat(14), color=ohmcolor

legend, /bottom,/right,['OHMs','Non-masing galaxies'],psym=[14,7], $
	color=[ohmcolor,concolor], charsize=1.2, charthick=4, thick=4

hor,2.3,linestyle=2, thick=3, color=defcolor

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; End

close,/all

end
