
;+
; NAME:
;       
;	SURV_BIREG_PLOT
;
; PURPOSE:
;
;	Plot results for bivariate survival regression
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
;       Written by K. Willett                Mar 09
;-

; Read in regression fit

silfile = '~/Astronomy/Research/Spitzer/surv/bisurv_regression/bireg_silicate.txt'

nlines = file_lines(silfile)
emptyarr = strarr(nlines)

openr, lun, silfile, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

intercept_line = where(strmid(emptyarr, 8, 9) eq 'Intercept')
slope_line     = where(strmid(emptyarr, 8, 5) eq 'Slope')

intercept = float(strmid(emptyarr(intercept_line),34,6))
intercept = intercept[0]

slope = float(strmid(emptyarr(slope_line),34,6))
slope_err = float(strmid(emptyarr(slope_line),63,6))
slope = slope[0]
slope_err = slope_err[0]

; Read in data

sildatafile = '~/Astronomy/Research/Spitzer/surv/bisurv_regression/bireg_data_silicate.txt'

readcol, sildatafile, censor, sil, logoh, skipline=1, format='i,f,f', /silent

censind = where(censor eq -1)
nocens = where(censor eq 0)

!p.multi=[0,1,1]

erase

x0 = 0.08
x2 = 0.95
x1 = (x2 + x0)/2

y0 = 0.10
y1 = 0.95

!x.style = 1
!y.style = 1

set_plot,'ps'
device, filename='~/Astronomy/Research/Spitzer/papers/surv_bireg_plot.ps', /color

thick = 5.0
cthick = 5.0
csize = 1.6
hsize = 250

blue = fsc_color("Blue")
red = fsc_color("Red")

x = fillarr(0.01,-10,10)

plot, sil, logoh, $
	/nodata, $
	xthick = thick, $
	ythick = thick, $
	thick = thick, $
	charthick = cthick, $
	charsize = csize, $
	xtitle='!7s!3!I9.7!N', $
	ytitle='log (L!IOH!N [L'+sunsymbol()+'])', $
	xr=[-0.5,3.9], /xstyle, $
	yr=[0,4.5], $
	ystyle = 9, $
	position = [x0,y0,x1,y1]

oplot, sil[nocens], logoh[nocens], psym=symcat(14)
oplot, sil[censind], logoh[censind], psym=7, thick = thick, symsize=1.2
arrow, sil[censind], logoh[censind], sil[censind], logoh[censind] - 0.3, /data, thick = thick, hsize=hsize

oplot, x, x*slope + intercept, thick = thick
oplot, x, x*(slope + slope_err) + intercept, linestyle=1, thick = thick
oplot, x, x*(slope - slope_err) + intercept, linestyle=1, thick = thick


; Read in regression fit

slope30file = '~/Astronomy/Research/Spitzer/surv/bisurv_regression/bireg_slope30.txt'

nlines = file_lines(slope30file)
strarr_slope = strarr(nlines)

openr, lun, slope30file, /get_lun
readf, lun, strarr_slope
close, lun & free_lun, lun

intercept_line = where(strmid(strarr_slope, 8, 9) eq 'Intercept')
slope_line     = where(strmid(strarr_slope, 8, 5) eq 'Slope')

s30_intercept = float(strmid(strarr_slope(intercept_line),34,6))
s30_intercept = s30_intercept[0]

s30_slope = float(strmid(strarr_slope(slope_line),34,6))
s30_slope_err = float(strmid(strarr_slope(slope_line),63,6))
s30_slope = s30_slope[0]
s30_slope_err = s30_slope_err[0]

; Read in data

slope30datafile = '~/Astronomy/Research/Spitzer/surv/bisurv_regression/bireg_data_slope30.txt'

readcol, slope30datafile, censor, slope30, logoh, skipline=1, format='i,f,f', /silent

censind = where(censor eq -1)
nocens = where(censor eq 0)

x = fillarr(0.01,-10,10)
plot, slope30, logoh, $
	/nodata, $
	/noerase, $
	xthick = thick, $
	ythick = thick, $
	thick = thick, $
	charthick = cthick, $
	charsize = csize, $
	xr = [-7.5,1], xstyle=1, $
	yr=[0,4.5], ystyle=1, $
	xtitle='!7a!3!I30-20!N', $
	yticks=1, $
	yminor=1, ytickname=replicate(' ',2), $
	position = [x1,y0,x2,y1]

oplot, slope30[nocens], logoh[nocens], psym=symcat(14)
oplot, slope30[censind], logoh[censind], psym=7, thick = thick, symsize=1.2
arrow, slope30[censind], logoh[censind], slope30[censind], logoh[censind] - 0.3, /data, thick = thick, hsize=hsize

oplot, x, x*s30_slope + s30_intercept, thick = thick
oplot, x, x*(s30_slope + s30_slope_err) + s30_intercept, linestyle=1, thick = thick
oplot, x, x*(s30_slope - s30_slope_err) + s30_intercept, linestyle=1, thick = thick

device,/close
set_plot,'x'

!p.multi=[0,1,1]

;plot,slope30,sil,psym=4


end
