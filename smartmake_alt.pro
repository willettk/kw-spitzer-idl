pro smartmake_alt, fname, ps = ps, write = write, nobonus = nobonus, $
	noclean = noclean, nocal = nocal, lores = lores, nospicehdr = nospicehdr, $
	xr=xr, yr=yr,stop=stop,plottitle=plottitle
;+
; NAME:
;       SMARTMAKE_ALT
;
;-

nobonus = 1

print,''
print, 'Running ',fname

; Set device to read in colors

device, decomposed = 0

	yellow = fsc_color('Yellow')
	green = fsc_color('Green')
	red = fsc_color('Red')
	blue = fsc_color('Blue')
	orange = fsc_color('Orange')

; Find directory and object name

tag, fname, dirtag 
targets, fname, redshift, obj

; Read in the spectrum

	specpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/coadd/'

	; Lores

	readcol, specpath+fname+'_sl1_coadd.tbl', $
		det_sl1, wave_sl1, flux_sl1, err_sl1, bit_sl1, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_sl2_coadd.tbl', $
		det_sl2, wave_sl2, flux_sl2, err_sl2, bit_sl2, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_sl3_coadd.tbl', $
		det_sl3, wave_sl3, flux_sl3, err_sl3, bit_sl3, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll1_coadd.tbl', $
		det_ll1, wave_ll1, flux_ll1, err_ll1, bit_ll1, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll2_coadd.tbl', $
		det_ll2, wave_ll2, flux_ll2, err_ll2, bit_ll2, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll3_coadd.tbl', $
		det_ll3, wave_ll3, flux_ll3, err_ll3, bit_ll3, format = 'i,f,f,f,i', skipline = 1, /silent

	; Hires

	readcol, specpath+fname+'_sh_coadd.tbl', $
		det_sh, wave_sh, flux_sh, err_sh, bit_sh, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_lh_coadd.tbl', $
		det_lh, wave_lh, flux_lh, err_lh, bit_lh, format = 'i,f,f,f,i', skipline = 1, /silent

; Read in the filters

pupath = '~/Astronomy/Research/Spitzer/spitzer/'

readcol, pupath+'bluePUtrans.txt', bwave, btrans, /silent
readcol, pupath+'redPUtrans.txt', rwave, rtrans, /silent


; LL1 module (18-35 um) is fixed

calflux_ll1 = flux_ll1

; Anchor LL2 to LL1

	first_pixel = wave_ll1(0)
	last_pixel= wave_ll2(n_elements(wave_ll2)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_ll2(n_elements(flux_ll2)-npix:n_elements(flux_ll2)-1)
		first_sec = calflux_ll1(0:npix-1)

		ll2_ll1_frac = mean(first_sec) / mean(last_sec)	
		
		calflux_ll2 = flux_ll2 * ll2_ll1_frac
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_ll2 ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin
			junk1 = where(abs(wave_ll2(overlap_pixels(0)+i) - wave_ll1) eq $
				min(abs(wave_ll2(overlap_pixels(0)+i) - wave_ll1)))
			overlap_bin(i) = fix(junk1(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_ll2(overlap_pixels) * (2. * i / nsteps)
			caltot(i) = abs(total(tempcal - calflux_ll1(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		ll2_ll1_frac = (2. * junkindex(0) / nsteps)

		calflux_ll2 = flux_ll2 * ll2_ll1_frac
	endelse

; Anchor SL1 to LL2

	first_pixel = wave_ll2(0)
	last_pixel= wave_sl1(n_elements(wave_sl1)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_sl1(n_elements(flux_sl1)-npix:n_elements(flux_sl1)-1)
		first_sec = calflux_ll2(0:npix-1)

		sl1_ll2_frac = mean(first_sec) / mean(last_sec)	
		
		calflux_sl1 = flux_sl1 * sl1_ll2_frac
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_sl1 ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin
			junk2 = where(abs(wave_sl1(overlap_pixels(0)+i) - wave_ll2) eq $
				min(abs(wave_sl1(overlap_pixels(0)+i) - wave_ll2)))
			overlap_bin(i) = fix(junk2(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_sl1(overlap_pixels) * (4. * i / nsteps)
			caltot(i) = abs(total(tempcal - calflux_ll2(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		sl1_ll2_frac = (4. * junkindex(0) / nsteps)

		calflux_sl1 = flux_sl1 * sl1_ll2_frac
	endelse

; Anchor SL2 to SL1

	first_pixel = wave_sl1(0)
	last_pixel= wave_sl2(n_elements(wave_sl2)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_sl2(n_elements(flux_sl2)-npix:n_elements(flux_sl2)-1)
		first_sec = calflux_sl1(0:npix-1)

		sl2_sl1_frac = mean(first_sec) / mean(last_sec)	
		
		calflux_sl2 = flux_sl2 * sl2_sl1_frac
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_sl2 ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin 
			junk3 = where(abs(wave_sl2(overlap_pixels(0)+i) - wave_sl1) eq $
				min(abs(wave_sl2(overlap_pixels(0)+i) - wave_sl1)))
			overlap_bin(i) = fix(junk3(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_sl2(overlap_pixels) * (2. * i / nsteps)
			caltot(i) = abs(total(tempcal - calflux_sl1(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		sl2_sl1_frac = (2. * junkindex(0) / nsteps)

		calflux_sl2 = flux_sl2 * sl2_sl1_frac
	endelse

; Calibrate bonus orders to first-order modules

calflux_sl3 = flux_sl3 * sl1_ll2_frac
calflux_ll3 = flux_ll3


; Align the HR modules

	first_pixel = wave_lh(0)
	last_pixel= wave_sh(n_elements(wave_sh)-1)

	if first_pixel gt last_pixel then begin		; If there is no overlap between modules, then align the last few pixels
	
		npix = 2				; Number of pixels used to align modules

		last_sec = flux_sh(n_elements(flux_sh)-npix:n_elements(flux_sh)-1)
		first_sec = flux_lh(0:npix-1)

		hr_frac = mean(first_sec) / mean(last_sec)	

		; Prevent unusually low values from setting the SH module to zeroes 

		if hr_frac eq 0 then begin
			hr_frac = 1.
			print,''
			print,'Unable to stitch hi-res modules due to low flux-scaling'
		endif
		
		hr_frac = 1d

			calflux_sh = flux_sh * hr_frac
			calflux_lh = flux_lh
	
	endif else begin				; If there is overlap, then move the modules up and down for the best fit
	
		overlap_pixels = where(wave_sh ge first_pixel)
		n_overlap = n_elements(overlap_pixels)
		overlap_bin = intarr(n_overlap)
		for i = 0, n_overlap - 1 do begin
			junk4 = where(abs(wave_sh(overlap_pixels(0)+i) - wave_lh) eq $
				min(abs(wave_sh(overlap_pixels(0)+i) - wave_lh)))
			overlap_bin(i) = fix(junk4(0))
		endfor
		nsteps = 3000
		caltot = fltarr(nsteps)
		for i = 0, nsteps - 1 do begin
			tempcal = flux_sh(overlap_pixels) * (2. * i / nsteps)
			caltot(i) = abs(total(tempcal - flux_lh(overlap_bin)))
		endfor
	
		junkindex = where(caltot eq min(caltot))
		hr_frac = (2. * junkindex(0) / nsteps)

		if hr_frac eq 0 then begin
			hr_frac = 1.
			print,''
			print,'Unable to stitch hi-res modules due to low flux-scaling'
		endif

		hr_frac = 1d

			calflux_sh = flux_sh * hr_frac
			calflux_lh = flux_lh

	endelse

; Align SH module to SL modules

junk1 = [wave_sl2,wave_sl1,wave_ll2,wave_ll1]
junk2 = [calflux_sl2,calflux_sl1,calflux_ll2,calflux_ll1]
junk3 = [err_sl1,err_sl1,err_ll2,err_ll1]
slind = sort(junk1)
slwave = junk1(slind)
slflux = junk2(slind)
slerr  = junk3(slind)

shflux = flux_sh
shwave = wave_sh
lhflux = flux_lh
lhwave = wave_lh

minlh = min(wave_lh)
maxlh = max(wave_lh)
minsh = min(wave_sh)
maxsh = max(wave_sh)
minsl = min(slwave)
maxsl = max(slwave)

minarr = max([minsh,minsl])
maxarr = min([maxsh,maxsl])

; Bin LH and LL

	binsize = 0.03
	wavearr_long = fillarr(binsize,min(lhwave),min([max(lhwave),max(slwave)]))
	llfluxbin = fltarr(n_elements(wavearr_long))
	lhfluxbin = fltarr(n_elements(wavearr_long))
	llerrbin  = fltarr(n_elements(wavearr_long))
	
	for i = 0, n_elements(wavearr_long)-1 do begin
		llfluxbin(i) = slflux(closeto(slwave,wavearr_long(i)))
		lhfluxbin(i) = lhflux(closeto(lhwave,wavearr_long(i)))
		llerrbin(i)  = slerr (closeto(slwave,wavearr_long(i)))
	endfor
	
	; Fit the high to the low res
	
	; Scaling
	
	midind = n_elements(wavearr_long) / 2
	firstguess = llfluxbin(midind) / lhfluxbin(midind)
	
	nsteps = 3000
	caltot_ll = fltarr(nsteps)
	for i = 0, nsteps - 1 do begin
		tempcal = lhfluxbin * (2. * i / nsteps * firstguess)
		caltot_ll(i) = total(((tempcal - llfluxbin)/llerrbin)^2)
	endfor
	
	junkindex_ll = where(caltot_ll eq min(caltot_ll))
	lh_frac = (2. * junkindex_ll(0) / nsteps * firstguess)
	
	; Shifting
	
	firstguess_shift = llfluxbin(midind) - lhfluxbin(midind)
	
	nsteps = 3000
	caltot_llshift = fltarr(nsteps)
	for i = 0, nsteps - 1 do begin
		tempcal = lhfluxbin + (2d * i/float(nsteps) - 1d) * firstguess_shift
		caltot_llshift(i) = total(((tempcal - llfluxbin)/llerrbin)^2)
	endfor
	
	ji_llshift = where(caltot_llshift eq min(caltot_llshift))
	lh_shift = (2d * ji_llshift(0) / nsteps - 1d) * firstguess_shift


; Bin each spectrum

binsize = 0.03
wavearr = fillarr(binsize,minarr,maxarr)
slfluxbin = fltarr(n_elements(wavearr))
shfluxbin = fltarr(n_elements(wavearr))
slerrbin  = fltarr(n_elements(wavearr))

for i = 0, n_elements(wavearr)-1 do begin
	slfluxbin(i) = slflux(closeto(slwave,wavearr(i)))
	shfluxbin(i) = shflux(closeto(shwave,wavearr(i)))
	slerrbin(i)  = slerr (closeto(slwave,wavearr(i)))
endfor
	
; Fit the high to the low res

; Scaling

midind = n_elements(wavearr) / 2
firstguess = slfluxbin(midind) / shfluxbin(midind)

nsteps = 3000
caltot = fltarr(nsteps)
for i = 0, nsteps - 1 do begin
	tempcal = shfluxbin * (2. * i / nsteps * firstguess)
	caltot(i) = total(((tempcal - slfluxbin)/slerrbin)^2)
endfor

junkindex = where(caltot eq min(caltot))
hr_frac = (2. * junkindex(0) / nsteps * firstguess)

; Shifting

firstguess_shift = slfluxbin(midind) - shfluxbin(midind)

nsteps = 3000
caltot_shift = fltarr(nsteps)
for i = 0, nsteps - 1 do begin
	tempcal = shfluxbin + (2d * i/float(nsteps) - 1d) * firstguess_shift
	caltot_shift(i) = total(((tempcal - slfluxbin)/slerrbin)^2)
endfor

ji_shift = where(caltot_shift eq min(caltot_shift))
hr_shift = (2d * ji_shift(0) / nsteps - 1d) * firstguess_shift

print,''
print,'sh_frac = ',hr_frac
print,'sh_shift = ',hr_shift
print,''
print,'lh_frac = ',lh_frac
print,'lh_shift = ',lh_shift
print,''

pspath = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/noskytest/'
if not keyword_set(plottitle) then plottitle = 'arch005_sh'

if keyword_set(ps) then begin
	set_plot,'ps'
	device, /landscape, /color, bits_per_pixel = 8, $
		filename = pspath+plottitle+'.ps'
	axiscolor = fsc_color('Black')
	bgcolor = fsc_color('White')
	cs = 1
	ct = 2
endif else begin
	axiscolor = fsc_color('White')
	bgcolor = fsc_color('Black')
	cs = 2
	ct = 1
endelse

if not keyword_set(xr) then xr=[8,20]
if not keyword_set(yr) then yr=[0,1]

plot,wavearr,slfluxbin,title=obj,xr=xr,/xstyle,yr=yr,/ystyle
oplot,wavearr,shfluxbin,color=fsc_color("Blue")
oplot,wavearr,shfluxbin * hr_frac,color=fsc_color("Red")
oplot,wavearr,shfluxbin + hr_shift,color=fsc_color("Green")
oplot,wavearr_long,llfluxbin
oplot,wavearr_long,lhfluxbin,color=fsc_color("Blue")
oplot,wavearr_long,lhfluxbin * lh_frac,color=fsc_color("Red")
oplot,wavearr_long,lhfluxbin + lh_shift,color=fsc_color("Green")
legend,/top,/left,['Low-res','Hi-res original','Hi-res scaled','Hi-res shifted'],linestyle=replicate(0,4),$
	color=[fsc_color("White"),fsc_color("Blue"),fsc_color("Red"),fsc_color("Green")]


xyouts,0.1,0.77,/normal,'Residuals',charsize=1.5
xyouts,0.1,0.70,/normal,'Scale: '+strtrim(caltot(junkindex)/n_elements(slfluxbin))
xyouts,0.1,0.65,/normal,'Shift: '+strtrim(caltot_shift(ji_shift)/n_elements(slfluxbin))
xyouts,0.1,0.60,/normal,'Scale: '+strtrim(caltot_ll(junkindex_ll)/n_elements(llfluxbin))
xyouts,0.1,0.55,/normal,'Shift: '+strtrim(caltot_llshift(ji_llshift)/n_elements(llfluxbin))
xyouts,0.27,0.67,'SH',/normal
xyouts,0.27,0.57,'LH',/normal

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

ggg = 0
if ggg eq 1 then begin
!p.multi = [0,1,1]
;!p.multi = [0,1,2]

; LORES spectra

pspath = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibration/'

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = pspath+fname+'_calflux.ps', /landscape, /color, bits_per_pixel = 8
	axiscolor = fsc_color('Black')
	bgcolor = fsc_color('White')
	cs = 1
	ct = 2
endif else begin
	axiscolor = fsc_color('White')
	bgcolor = fsc_color('Black')
	cs = 2
	ct = 1
endelse

if not keyword_set(xr) then xr=[4,40]
if not keyword_set(yr) then yr=[0,40]

plot, wave_sl1, flux_sl1, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
	xr = xr, /xstyle, $
;	yr = [1d-4,1d1], /ystyle, $
	yr = yr, /ystyle, $
;	/xlog,/ylog, $
	title = obj, $
	charsize = cs, $
	color = axiscolor, $
	background = bgcolor, $
	linestyle = 2, $
	thick = ct, $
	charthick = ct, $
	/nodata


oplot, wave_sl1, flux_sl1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_sl2, flux_sl2, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll1, flux_ll1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll2, flux_ll2, color = axiscolor, linestyle = 2, thick = ct

oplot, wave_ll2, calflux_ll2, color = blue, thick = ct
oplot, wave_ll1, calflux_ll1, color = red, thick = ct
oplot, wave_sl2, calflux_sl2, color = yellow, thick = ct
oplot, wave_sl1, calflux_sl1, color = green, thick = ct

if not keyword_set(nobonus) then begin
	oplot, wave_sl3, flux_sl3, color = axiscolor, linestyle = 2, thick = ct
	oplot, wave_ll3, flux_ll3, color = axiscolor, linestyle = 2, thick = ct
	oplot, wave_ll3, calflux_ll3, color = orange, thick = ct
	oplot, wave_sl3, calflux_sl3, color = orange, thick = ct
endif

oplot, rwave, rtrans / max(rtrans) * max(flux_ll1), linestyle = 1, color=red, thick = ct
oplot, bwave, btrans / max(btrans) * max(flux_ll1), linestyle = 1, color=blue, thick = ct


xyouts, 0.2, 0.85, 'Lo-res', /normal, charsize = cs, charthick = ct, color = axiscolor

; HIRES spectra

;plot, wave_sh, flux_sh, $
;	xtitle = 'Wavelength [!7l!3m]', $
;	ytitle = 'Flux [Jy]', $
;	xr = [5,40], $
;	yr = [0,max(flux_lh)], $
;	title = obj, $
;	charsize = cs, $
;	linestyle = 2, $
;	color = axiscolor, $
;	background = bgcolor, $
;	thick = ct, $
;	charthick = ct, $
;	/nodata

;oplot, wave_sh, flux_sh, color = axiscolor, linestyle = 1, thick = ct, psym = 10
;oplot, wave_lh, flux_lh, color = axiscolor, linestyle = 1, thick = ct, psym = 10
oplot, wave_sh, calflux_sh+0.5, color = blue, linestyle = 0, thick = ct, psym = 10
oplot, wave_sh, calflux_sh*1.5, color = fsc_color("Cyan"), linestyle = 0, thick = ct, psym = 10
oplot, wave_lh, calflux_lh, color = red, linestyle = 0, thick = ct, psym = 10

oplot, rwave, rtrans / max(rtrans) * max(flux_lh), linestyle = 1, color=red, thick = ct
oplot, bwave, btrans / max(btrans) * max(flux_lh), linestyle = 1, color=blue, thick = ct

xyouts, 0.2, 0.35, 'Hi-res', /normal, charsize = cs, charthick = ct, color = axiscolor
xyouts, 0.8, 0.95, fname, /normal, charsize = cs, charthick = ct, color=axiscolor

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif


; Open the individual nods and calibrate them using the same scale factors as the coadded spectra

nodpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/nods/'

	; Lores

	readcol, nodpath+fname+'_sl1_1p.tbl', $
		det_sl1_1p, wave_sl1_1p, flux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl1_2p.tbl', $
		det_sl1_2p, wave_sl1_2p, flux_sl1_2p, err_sl1_2p, bit_sl1_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl2_1p.tbl', $
		det_sl2_1p, wave_sl2_1p, flux_sl2_1p, err_sl2_1p, bit_sl2_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl2_2p.tbl', $
		det_sl2_2p, wave_sl2_2p, flux_sl2_2p, err_sl2_2p, bit_sl2_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl3_1p.tbl', $
		det_sl3_1p, wave_sl3_1p, flux_sl3_1p, err_sl3_1p, bit_sl3_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl3_2p.tbl', $
		det_sl3_2p, wave_sl3_2p, flux_sl3_2p, err_sl3_2p, bit_sl3_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll1_1p.tbl', $
		det_ll1_1p, wave_ll1_1p, flux_ll1_1p, err_ll1_1p, bit_ll1_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll1_2p.tbl', $
		det_ll1_2p, wave_ll1_2p, flux_ll1_2p, err_ll1_2p, bit_ll1_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll2_1p.tbl', $
		det_ll2_1p, wave_ll2_1p, flux_ll2_1p, err_ll2_1p, bit_ll2_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll2_2p.tbl', $
		det_ll2_2p, wave_ll2_2p, flux_ll2_2p, err_ll2_2p, bit_ll2_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll3_1p.tbl', $
		det_ll3_1p, wave_ll3_1p, flux_ll3_1p, err_ll3_1p, bit_ll3_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll3_2p.tbl', $
		det_ll3_2p, wave_ll3_2p, flux_ll3_2p, err_ll3_2p, bit_ll3_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	; Hires

	readcol, nodpath+fname+'_sh_1p.tbl', $
		det_sh_1p, wave_sh_1p, flux_sh_1p, err_sh_1p, bit_sh_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sh_2p.tbl', $
		det_sh_2p, wave_sh_2p, flux_sh_2p, err_sh_2p, bit_sh_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_lh_1p.tbl', $
		det_lh_1p, wave_lh_1p, flux_lh_1p, err_lh_1p, bit_lh_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_lh_2p.tbl', $
		det_lh_2p, wave_lh_2p, flux_lh_2p, err_lh_2p, bit_lh_2p, format = 'i,f,f,f,i', skipline = 1, /silent

if not keyword_set(nocal) then begin

; Calibrate the nods to peakups

; Lo-res

	calflux_sl1_1p = flux_sl1_1p * sl1_ll2_frac
	calflux_sl2_1p = flux_sl2_1p * sl2_sl1_frac
	calflux_sl3_1p = flux_sl3_1p * sl1_ll2_frac
	calflux_ll1_1p = flux_ll1_1p
	calflux_ll3_1p = flux_ll3_1p 
	calflux_ll2_1p = flux_ll2_1p * ll2_ll1_frac
	
	calflux_sl1_2p = flux_sl1_2p * sl1_ll2_frac
	calflux_sl2_2p = flux_sl2_2p * sl2_sl1_frac
	calflux_sl3_2p = flux_ll1_2p * sl1_ll2_frac
	calflux_ll1_2p = flux_ll1_2p 
	calflux_ll3_2p = flux_ll3_2p 
	calflux_ll2_2p = flux_ll2_2p * ll2_ll1_frac
	
; Hi-res

calflux_sh_1p = flux_sh_1p
calflux_lh_1p = flux_lh_1p

calflux_sh_2p = flux_sh_2p
calflux_lh_2p = flux_lh_2p

; Option for non-calibrated spectra (or when PUs are suspect)

endif

if keyword_set(nocal) then begin

	calflux_sl1 = flux_sl1
	calflux_sl2 = flux_sl2
	calflux_sl3 = flux_sl3
	calflux_ll1 = flux_ll1
	calflux_ll2 = flux_ll2
	calflux_ll3 = flux_ll3
	calflux_sh  = flux_sh
	calflux_lh  = flux_lh
	
	calflux_sl1_1p = flux_sl1_1p
	calflux_sl2_1p = flux_sl2_1p
	calflux_sl3_1p = flux_sl3_1p
	calflux_ll1_1p = flux_ll1_1p
	calflux_ll2_1p = flux_ll2_1p
	calflux_ll3_1p = flux_ll3_1p
	calflux_sh_1p  = flux_sh_1p
	calflux_lh_1p  = flux_lh_1p
	
	calflux_sl1_2p = flux_sl1_2p
	calflux_sl2_2p = flux_sl2_2p
	calflux_sl3_2p = flux_sl3_2p
	calflux_ll1_2p = flux_ll1_2p
	calflux_ll2_2p = flux_ll2_2p
	calflux_ll3_2p = flux_ll3_2p
	calflux_sh_2p  = flux_sh_2p
	calflux_lh_2p  = flux_lh_2p

endif

if keyword_set(write) then begin

	; Writing stitched spectra to disk
	
	writepath_coadd = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/stitched/coadd/'
	writepath_nod   = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/stitched/nods/'

	; Coadded

	forprint, det_sl1, wave_sl1, calflux_sl1, err_sl1, bit_sl1, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl1_cal.tbl', /silent

	forprint, det_sl2, wave_sl2, calflux_sl2, err_sl2, bit_sl2, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl2_cal.tbl', /silent

	forprint, det_sl3, wave_sl3, calflux_sl3, err_sl3, bit_sl3, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl3_cal.tbl', /silent

	forprint, det_ll1, wave_ll1, calflux_ll1, err_ll1, bit_ll1, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll1_cal.tbl', /silent

	forprint, det_ll2, wave_ll2, calflux_ll2, err_ll2, bit_ll2, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll2_cal.tbl', /silent

	forprint, det_ll3, wave_ll3, calflux_ll3, err_ll3, bit_ll3, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll3_cal.tbl', /silent

	forprint, det_sh, wave_sh, calflux_sh, err_sh, bit_sh, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sh_cal.tbl', /silent

	forprint, det_lh, wave_lh, calflux_lh, err_lh, bit_lh, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_lh_cal.tbl', /silent

	; Nods

	forprint, det_sl1_1p, wave_sl1_1p, calflux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl1_1p_cal.tbl', /silent

	forprint, det_sl1_1p, wave_sl1_1p, calflux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl1_2p_cal.tbl', /silent

	forprint, det_sl2_1p, wave_sl2_1p, calflux_sl2_1p, err_sl2_1p, bit_sl2_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl2_1p_cal.tbl', /silent

	forprint, det_sl2_2p, wave_sl2_2p, calflux_sl2_2p, err_sl2_2p, bit_sl2_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl2_2p_cal.tbl', /silent

	forprint, det_sl3_1p, wave_sl3_1p, calflux_sl3_1p, err_sl3_1p, bit_sl3_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl3_1p_cal.tbl', /silent

	forprint, det_sl3_2p, wave_sl3_2p, calflux_sl3_2p, err_sl3_2p, bit_sl3_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl3_2p_cal.tbl', /silent

	forprint, det_ll1_1p, wave_ll1_1p, calflux_ll1_1p, err_ll1_1p, bit_ll1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll1_1p_cal.tbl', /silent

	forprint, det_ll1_2p, wave_ll1_2p, calflux_ll1_2p, err_ll1_2p, bit_ll1_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll1_2p_cal.tbl', /silent

	forprint, det_ll2_1p, wave_ll2_1p, calflux_ll2_1p, err_ll2_1p, bit_ll2_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll2_1p_cal.tbl', /silent

	forprint, det_ll2_2p, wave_ll2_2p, calflux_ll2_2p, err_ll2_2p, bit_ll2_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll2_2p_cal.tbl', /silent

	forprint, det_ll3_1p, wave_ll3_1p, calflux_ll3_1p, err_ll3_1p, bit_ll3_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll3_1p_cal.tbl', /silent

	forprint, det_ll3_2p, wave_ll3_2p, calflux_ll3_2p, err_ll3_2p, bit_ll3_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll3_2p_cal.tbl', /silent

	forprint, det_sh_1p, wave_sh_1p, calflux_sh_1p, err_sh_1p, bit_sh_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sh_1p_cal.tbl', /silent

	forprint, det_sh_2p, wave_sh_2p, calflux_sh_2p, err_sh_2p, bit_sh_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sh_2p_cal.tbl', /silent

	forprint, det_lh_1p, wave_lh_1p, calflux_lh_1p, err_lh_1p, bit_lh_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_lh_1p_cal.tbl', /silent

	forprint, det_lh_2p, wave_lh_2p, calflux_lh_2p, err_lh_2p, bit_lh_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_lh_2p_cal.tbl', /silent

endif

print, 'LL2 to LL1 scale factor: ', ll2_ll1_frac
print, 'SL1 to LL2 scale factor: ', sl1_ll2_frac
print, 'SL2 to SL1 scale factor: ', sl2_sl1_frac
print, 'SH to LH scale factor: ', hr_frac
print, 'Calibrated spectra to peakup fluxes for ',fname

;if not keyword_set(nospicehdr) then spicehdr, fname, readdir = 'stitched'

; Save the scaling fractions for later use

save, hr_frac, ll2_ll1_frac, sl1_ll2_frac, sl2_sl1_frac, $
	filename='~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/scaling_frac/'+fname+'.sav'

endif

; End program

end
