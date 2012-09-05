pro pahspl, fname, pahew62, pah62err, pahew11, pah11err, pah62flux, pah11flux, spl = spl, ps = ps, noplot = noplot
;+
; NAME:
;       PAHSPL
;
; PURPOSE:
;	Measure the 6.2 EW from calibrated lo-res IRS spectra
;
; INPUTS:
;
;	FNAME - 	tag of object to reduce (eg, 'mega005')
;
;	SPL - 		integer giving type of spline fit for used for the 9.7 um continuum
;
;			1: continuum-dominated
;			2: PAH-dominated
;			3: absorption-dominated (DEFAULT)
;			4: absorption-dominated w/no PAH flux at 8 um
;			5: custom spline at user-defined pivot points 
;
;	- calibrated lo-res spectra from all four modules via SPEC2PUFLUX
;
; OUTPUTS:
;
; KEYWORDS:
;
;	PS - hard copies of spectra with spline fits and continuum removed
;
; EXAMPLE:
;	IDL> 
;
; REQUIRES:
;
;	TAG.pro
;	READCOL,pro
;
; NOTES:
;
; REVISION HISTORY
;       Written by K. Willett                Sep 2007
;-


device, window_state = state

; Read in the full LR spectrum

tag,fname,dirtag

readcol, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_sl1_cal.tbl', $
	sl1_order, sl1_wave, sl1_flux, sl1_err, sl1_bit, format = 'i,f,f,f,i', /silent

readcol, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_sl2_cal.tbl', $
	sl2_order, sl2_wave, sl2_flux, sl2_err, sl2_bit, format = 'i,f,f,f,i', /silent

readcol, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_sl3_cal.tbl', $
	sl3_order, sl3_wave, sl3_flux, sl3_err, sl3_bit, format = 'i,f,f,f,i', /silent

readcol, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_ll1_cal.tbl', $
	ll1_order, ll1_wave, ll1_flux, ll1_err, ll1_bit, format = 'i,f,f,f,i', /silent

readcol, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_ll2_cal.tbl', $
	ll2_order, ll2_wave, ll2_flux, ll2_err, ll2_bit, format = 'i,f,f,f,i', /silent

readcol, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'+fname+'_ll3_cal.tbl', $
	ll3_order, ll3_wave, ll3_flux, ll3_err, ll3_bit, format = 'i,f,f,f,i', /silent


; Stitch together lo-res spectra into one 

wave1 = [sl2_wave, sl3_wave, sl1_wave, ll2_wave, ll3_wave, ll1_wave]
flux1 = [sl2_flux, sl3_flux, sl1_flux, ll2_flux, ll3_flux, ll1_flux]

waves = sort(wave1)

wave = wave1(waves)
flux = flux1(waves)

; Shift to rest-frame wavelengths

targets, fname, redshift, obj
wave = wave / (1. + redshift)


; Spline fits for the 9.7 um continuum

if not keyword_set(spl) then spl = 2
case spl of
	1: begin
		m1 = alog10(flux(closetomed(wave,flux,7.5)) / flux(closetomed(wave,flux,5.0))) $
			/ alog10(wave(closetomed(wave,flux,7.5)) / wave(closetomed(wave,flux,5.0)))
	   	b1 = alog10(flux(closetomed(wave,flux,7.5))) - m1 * alog10(wave(closetomed(wave,flux,7.5)))
	   	powerfit1 = 10^(m1 * alog10(wave(0:closetomed(wave,flux,7.5))) + b1)
		m2 = alog10(flux(closetomed(wave,flux,31.5)) / flux(closetomed(wave,flux,26.1))) $
			/ alog10(wave(closetomed(wave,flux,31.5)) / wave(closetomed(wave,flux,26.1)))
	   	b2 = alog10(flux(closetomed(wave,flux,31.5))) - m2 * alog10(wave(closetomed(wave,flux,31.5)))
	   	powerfit2 = 10^(m2 * alog10(wave(closetomed(wave,flux,26.1):n_elements(wave)-1)) + b2)
		wcont = [wave(closetomed(wave,flux,7.5)),wave(closetomed(wave,flux,14.0)),wave(closetomed(wave,flux,26.1))]
		fcont = [flux(closetomed(wave,flux,7.5)),flux(closetomed(wave,flux,14.0)),flux(closetomed(wave,flux,26.1))]
		yspl = spl_init(wcont,fcont)
	   	spltemp = spl_interp(wcont, fcont, yspl, wave(closetomed(wave,flux,7.5):closetomed(wave,flux,26.1)))
		splresult = [powerfit1(0:n_elements(powerfit1)-2), spltemp, powerfit2(1:n_elements(powerfit2)-1)]
		pivotpoints = [closetomed(wave,flux,7.5),closetomed(wave,flux,14.0),closetomed(wave,flux,26.1)]
	   end
	2: begin 
		if wave(n_elements(wave)-1) lt 30.0 then endind = 28.0 else endind = 30.0
		wcont = wave([closetomed(wave,flux,14.5),closetomed(wave,flux,26.0),closetomed(wave,flux,endind)]) 
	   	fcont = flux([closetomed(wave,flux,14.5),closetomed(wave,flux,26.0),closetomed(wave,flux,endind)]) 
	   	yspl = spl_init(wcont, fcont)
	   	spltemp = spl_interp(wcont, fcont, yspl, wave(closetomed(wave,flux,14.5):n_elements(wave)-1))
	   	m = alog10(flux(closetomed(wave,flux,14.5)) / flux(closetomed(wave,flux,5.5))) $
			/ alog10(wave(closetomed(wave,flux,14.5)) / wave(closetomed(wave,flux,5.5)))
	   	b = alog10(flux(closetomed(wave,flux,14.5))) - m * alog10(wave(closetomed(wave,flux,14.5)))
	   	powerfit = 10^(m * alog10(wave(0:closetomed(wave,flux,14.5)-1)) + b)
	   	splresult = [powerfit,spltemp]
		pivotpoints = [closetomed(wave,flux,5.5),closetomed(wave,flux,14.5),closetomed(wave,flux,26.0),closetomed(wave,flux,endind)]
	   end
	3: begin
		wcont = wave([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]) 
		fcont = flux([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]) 
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]
	   end
	4: begin
		wcont = wave([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,7.8),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]) 
		fcont = flux([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,7.8),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]) 
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),closetomed(wave,flux,7.8),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]
	   end
	5: begin
		pp = [5,6,11.7,14,27]
		wcont = fltarr(n_elements(pp))
		fcont = fltarr(n_elements(pp))
		pivotpoints = fltarr(n_elements(pp))
		for j = 0, n_elements(pp)-1 do begin
			wcont(j) = wave(closetomed(wave,flux,pp(j)))
			fcont(j) = flux(closetomed(wave,flux,pp(j)))
			pivotpoints(j) = closetomed(wave,flux,pp(j))
		endfor
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
	   end
endcase


; Add feature to measure the PAH 6.2 um EW for placement in Spoon's classification scheme

; Designate continuum regions on either side

if closeto(wave,5.15) ne 0 then pahpiv = [closetomed(wave,flux,5.15),closetomed(wave,flux,5.55),$
	closetomed(wave,flux,5.95),closetomed(wave,flux,6.55),closetomed(wave,flux,7.1)] else $
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
	addpah(1,i-bpah) = wave(i+1) - wave(i)
endfor

pah62fluxarr = addpah(0,*) * addpah(1,*)
pah62flux = total(pah62fluxarr)
pahew62 = pah62flux / pahspl(closetomed(wave,flux,6.2))
baseerr62 = stddev(pahspl(closetomed(wave,flux,6.1):closetomed(wave,flux,6.3)))
pah62err = baseerr62 * pah62flux / (pahspl(closetomed(wave,flux,6.2)))^2

; Also define a PAH EW using the spline-based silicate continuum, instead of the measured continuum at 6.2 um

pahew62_alt = pah62flux / splresult(closetomed(wave,flux,6.2))

; Add feature to measure the PAH 11 um EW

pahpiv11 = [closetomed(wave,flux,10.1), closetomed(wave,flux,10.9),closetomed(wave,flux,11.8),closetomed(wave,flux,12.4)]

bpah11 = closetomed(wave,flux,10.9)
epah11 = closetomed(wave,flux,11.8)

; Spline fit to PAH continuum

wp11 = wave(pahpiv11)
fp11 = flux(pahpiv11)
pspl11 = spl_init(wp11,fp11)
pahspl11 = spl_interp(wp11,fp11,pspl11,wave)


newpah11 = flux - pahspl11
addpah11 = fltarr(2,epah11 - bpah11 + 1)
for i = bpah11, epah11 do begin
	addpah11(0,i-bpah11) = newpah11(i)
	addpah11(1,i-bpah11) = wave(i+1) - wave(i)
endfor

pah11fluxarr = addpah11(0,*) * addpah11(1,*)
pah11flux = total(pah11fluxarr)
pahew11 = pah11flux / pahspl11(closetomed(wave,flux,11.3))
baseerr11 = stddev(pahspl11(closetomed(wave,flux,11.2):closetomed(wave,flux,11.4)))
pah11err = baseerr11 * pah11flux / (pahspl11(closetomed(wave,flux,11.3)))^2

; Plot the results
if not keyword_set(noplot) then begin

!p.multi = [0,2,1]
if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/pah/'+fname+'_pah.ps', /color
	cs = 1
	ls = 2
endif else begin
	cs = 2
	ls = 1
endelse

plot, wave, flux, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
;	/xlog, $
	/ylog, $
	xrange = [5,8], $
	yrange = [1d-4,1d-1], $
	/xstyle, $
	charsize = cs, $
	thick = ls, $
	charthick = ls, $
	title = obj

oplot, wp, fp, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
oplot, wave, pahspl, color=fsc_color("Green"), thick = ls
ver, wave(bpah), linestyle = 1
ver, wave(epah), linestyle = 1

xyouts, 5.2, 0.07, 'PAH 6.2 um EW  : '+string(pahew62,format='(f6.3)')+' um', /data, charsize = cs
xyouts, 5.2, 0.05, 'PAH 6.2 um flux : '+string(pah62flux*1d3,format='(f6.3)')+' mJy', /data, charsize = cs

plot, wave, flux, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
;	/xlog, $
	/ylog, $
	xrange = [5,15], $
	yrange = [1d-4,1d0], $
	/xstyle, $
	charsize = cs, $
	thick = ls, $
	charthick = ls, $
	title = obj

oplot, wp11, fp11, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
oplot, wave, pahspl11, color=fsc_color("Green"), thick = ls
ver, wave(bpah11), linestyle = 1
ver, wave(epah11), linestyle = 1

xyouts, 5.2, 0.7, 'PAH 11um EW  : '+string(pahew11,format='(f6.3)')+' um', /data, charsize = cs
xyouts, 5.2, 0.5, 'PAH 11um flux : '+string(pah11flux*1d3,format='(f7.3)')+' mJy', /data, charsize = cs

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif
!p.multi=[0,1,1]

; Print results to screen

print,''
print,'PAH 6.2 EW [10^3 um] :  ',pahew62*1d3,' +- ',pah62err*1d3
;print,'PAH 6.2 EW (ice):  ',pahew62_alt
;print,''
print,'PAH 11.3 EW [10^3 um]:  ',pahew11*1d3,' +- ',pah11err*1d3
;print,''

endif
;stop
end

