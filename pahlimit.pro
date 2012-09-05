;pro pahlimit
;+
; NAME:
;       
;	PAHLIMIT
;
; PURPOSE:
;
;	Sets upper limit on the equivalent width of the 6.2 um PAH feature in Spitzer LR spectra
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
;	IDL> 
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                
;-

; Color definitions

red = fsc_color("Red")
blue = fsc_color("Dodger Blue")
defcolor = fsc_color("White")
green = fsc_color("Green")
orange = fsc_color("Orange")

; Read in the spectra

fname = 'cso005'

datapath = '~/Astronomy/Research/Spitzer/cso/data/structures/'
restore, datapath+fname+'.sav'

wave = sed.wave_lr
flux = sed.flux_lr
err = sed.err_lr

; Pivot points for 6.2 um spline fit

pahpiv = [closetomed(wave,flux,5.55),closetomed(wave,flux,5.95),closetomed(wave,flux,6.55),closetomed(wave,flux,7.1)]

bpah = closetomed(wave,flux,5.95)
epah = closetomed(wave,flux,6.55)

; Spline fit to PAH continuum

wp = wave(pahpiv)
fp = flux(pahpiv)
pspl1 = spl_init(wp,fp)
pahspl = spl_interp(wp,fp,pspl1,wave)


newpah = flux - pahspl
addpah = fltarr(2,epah - bpah + 1)
for i = bpah, epah do begin
	addpah(0,i-bpah) = newpah(i)
	addpah(1,i-bpah) = (wave(i+1) - wave(i)) * 3d14 / (wave(i))^2
endfor

; Also try a polynomial fit between 5.95 and 6.55 um

expr = 'p[0] + x*p[1] + x*x*p[2]'
start = [1d-2,1d-2,1d-2]
range = where(wave gt 5.5 and wave lt 6.55)
result = mpfitexpr(expr,wave(range),flux(range),err(range),start,/quiet)

fit = result(0) + wave*result(1) + wave^2 * result(2)
polyfit = flux - fit

; Drude profile

baseline = pahspl(closetomed(wave,flux,6.2))
lambdar = 6.22
gammar = 0.030
br = 0.018 - baseline
drude = -1 * br * gammar^2 / ((wave/lambdar - lambdar/wave)^2+gammar^2)

; Crude approximation

newpah2 = drude 
addpah = fltarr(2,epah - bpah + 1)
for i = bpah, epah do begin
	addpah(0,i-bpah) = newpah2(i)
	addpah(1,i-bpah) = (wave(i+1) - wave(i)) * 3d14 / (wave(i))^2
endfor



pah62fluxarr = addpah(0,*) * addpah(1,*) 								; Flux density in Jy, d_nu in Hz
pah62flux = total(pah62fluxarr) * 1d-30									; Flux in W/cm^2
pahew62 = pah62flux / pahspl(closetomed(wave,flux,6.2)) * (6.2)^2 / 3d14 * 1d30				; EW in um
baseerr62 = stddev(flux(closetomed(wave,flux,6.1):closetomed(wave,flux,6.3)))				; Estimate error in baseline flux density [Jy]
range = where(wave gt 6.1 and wave lt 6.3)
baseerr62 = max(flux(range)) - min(flux(range))
pahew62err = baseerr62 * pah62flux / (pahspl(closetomed(wave,flux,6.2)))^2 * (6.2)^2/ 3d14*1d30		; Error in equivalent width

print,'PAH 6.2 EW [um] = ',pahew62,' +- ',pahew62err

; Plot the spectrum and theoretical profile with best fits

!p.multi=[0,2,1]

plot, wave, flux, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
;	/xlog, $
	/ylog, $
	xrange = [5,8], $
	yrange = [5d-3,5d-2], /ystyle, $
	/xstyle, $
	charsize = cs, $
	thick = ls, $
	charthick = ls, $
	title = sed.obj

oplot, wp, fp, psym = 4, symsize = 2, color=orange, thick = ls
oplot, wave, pahspl, color=green, thick = ls
oplot, wave, fit, color=red, thick = ls
oplot, wave, drude + baseline, color=blue, thick = ls
ver, wave(bpah), linestyle = 1
ver, wave(epah), linestyle = 1
legend,/top,/right,$
	['Spectrum','Spline fit','Polynomial fit','Drude'],$
	color=[defcolor,green,red,blue], $
	linestyle = intarr(4)


;xyouts, 5.2, 0.07, 'PAH 6.2 um EW  : '+string(pahew62,format='(f6.3)')+' um', /data, charsize = cs
;xyouts, 5.2, 0.05, 'PAH 6.2 um flux : '+string(pah62flux,format='(e9.3)')+' W/cm^2', /data, charsize = cs

; Plot the baseline-subtracted spectrum and theoretical profile

plot, wave, newpah, $
	/nodata, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
;	/xlog, $
;	/ylog, $
	xrange = [5,8], $
	yrange = [-1d-2,1d-2], $
	/xstyle, $
	charsize = cs, $
	thick = ls, $
	color = defcolor, $
	charthick = ls, $
	title = sed.obj

oplot, wave, newpah, color=green, thick=ls
oplot, wave, drude, color=blue, thick=ls
oplot, wave, polyfit, color=red, thick = ls

legend,/top,/right,$
	['Spline subtraction','Polynomial subtraction','Drude'],$
	color=[green,red,blue], $
	linestyle = intarr(3)

ver,wave(bpah),linestyle=1
ver,wave(epah),linestyle=1

stop
end
