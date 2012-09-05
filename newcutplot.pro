pro newcutplot, ps=ps

; NEWCUTPLOT

; Plot limits for max L_OH for new objects w/Spitzer spectra
; KW, 26 Feb 08

; Data is from Kent, Braatz, and Darling high-z survey of OHMs at GBT

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
oh_arch = float(archdat('logoh'))
ohz = float(transpose(ohmdat('redshift')))
archz = float(archdat('redshift'))

; Plot the data

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/control/newcutplot.ps',/landscape,/color
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

mysym = symcat(16)
red=fsc_color("Red")
green=fsc_color("Green")
blue=fsc_color("Blue")
orange=fsc_color("Orange")

!p.multi=0

plot, kz, klohmax, $
	xrange=[0,0.25], $
	yrange=[1d-4,1d5],/ylog, $
	xtitle='z', $
	ytitle='L!IOH!N!Emax!N', $
	charsize=2, $
	psym=mysym, $
	title='Non-detections of OH with full IRS spectra', $
	/nodata


oplot, dz, dlohmax, color=defcolor, psym=mysym
oplot, bz, blohmax, color=blue, psym=mysym
oplot, ssz, sslohmax, color=green, psym=mysym
oplot, kz, klohmax, color=orange, psym=mysym

oplot, ohz, 10^oh_reg, color=red, psym=symcat(15)
oplot, archz, 10^oh_arch, color=red, psym=symcat(15)

legend, /bottom,/right,['Baan','Staveley-Smith','Darling','Kent','OHMs'],psym=[replicate(mysym,4),symcat(15)], $
	color=[blue,green,defcolor,orange,red]

ver, 0.1, linestyle=2
hor, 1d2, linestyle=2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

;stop
end
