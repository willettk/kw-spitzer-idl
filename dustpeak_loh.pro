pro dustpeak_loh, ps=ps
;+
; NAME:
;       
;	DUSTPEAK_LOH
;
; PURPOSE:
;
;	Plot the dust temperature from DUSTY/PAHFIT fits against the OH luminosity
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
;	IDL> dustpeak_loh
;
; NOTES:
;
;	CLUMPY_OHM plots the same thing - this program is likely obsolete. 
;
; REVISION HISTORY
;       Written by K. Willett                Jan 10
;-

; Load L_OH

loh = [transpose(float(archdat('logoh'))),transpose(float(ohmdat('logoh')))]
loh_con = [transpose(float(condat('logoh')))]

; Load dust temps from BB fit

;dtemp_both = [transpose(float(archdat('dtemp'))),transpose(float(ohmdat('dtemp')))]
;dtemp_both_con = [transpose(float(condat('dtemp')))]
;
;dtemp     = dtemp_both[*,0]
;dtemp_err = dtemp_both[*,1]
;
;dtemp_con     = dtemp_both_con[*,0]
;dtemp_con_err = dtemp_both_con[*,1]

; Load dust temps from the DUSTY fits

restore, '~/Astronomy/Research/Spitzer/dusty_pahfit_results.sav'	; [arch, mega]

dtemp = ohm_tdustarr
dtemp_con = con_tdustarr

; Plot

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/dustpeak_loh.ps', /portrait, /color
	defcolor=fsc_color("Black")
	cs = 1.5
	th = 5
	hsize = 250
endif else begin
	defcolor=fsc_color("White")
	cs = 2
	th = 1
	hsize = 30
endelse

erase

x0 = 0.15
x1 = 0.8
x2 = 0.95

y0 = 0.15
y1 = 0.95
y2 = 0.95

!x.style = 1
!y.style = 1

red = fsc_color("Red")
blue = fsc_color("Blue")

plot, [loh, loh_con], [dtemp, dtemp_con], $
	position = [x0,y0,x1,y1], $
	color=defcolor, $
	xr=[0,4.5], /xstyle, $
	yr=[0,130], $
	/nodata, $
	thick = th, $
	xthick = th, $
	ythick = th, $
	charthick = th, $
	charsize = cs, $
	ytitle='Dust temperature [K]', $
	xtitle='log (L!IOH!N [L'+sunsymbol()+'])'

;oploterror, loh, dtemp, dtemp_err, $
oplot, loh, dtemp, $
	psym=symcat(16), $
;	errcolor=red, $
	color=red

;oploterror, loh_con, dtemp_con, dtemp_con_err, $
oplot, loh_con, dtemp_con, $
	psym=symcat(16), $
;	errcolor=blue, $
	color=blue

;xyouts, 2.6, 110, /data, '!7q!3!IOHM!N ='+string(correlate(dtemp,loh),format='(f5.2)'), charsize=2
;xyouts, 2.6, 100,  /data, '!7q!3!Icon!N ='+string(correlate(dtemp_con,loh_con),format='(f5.2)'), charsize=2

	arrow, loh_con, dtemp_con, loh_con - 0.2, dtemp_con, /data, color=blue, thick=th, hsize=hsize

; Linear fit

xarr = fillarr(1d-1,0,200)

expr='p[0]+x*p[1]'
start = [0,1]

;fit = mpfitexpr(expr, loh, dtemp, dtemp_err, start,/quiet, perror = perror)

;oplot, xarr, fit[0] + xarr * fit[1], linestyle=2

; Silicate on side

plot, indgen(10), /nodata, $
	/xstyle, /ystyle, $
	xr = [0,30], $
	yr = [0,130], $
	yticks = 1, $
	ytickv = [0,130], $
	ytickname = replicate(' ',3), $
	xticks = 3, $
	xtickv = fillarr(10,0,30), $
	xtickname = ['0','10','20','30'], $
	thick = th, $
	charthick = th, $
	xthick = th, $
	ythick = th, $
	charsize = cs, $
	position = [x1,y0,x2,y1],$
	/noerase

bs = 10
plothist, dtemp,     obin, ocount, /noplot, bin = bs
plothist, dtemp_con, cbin, ccount, /noplot, bin = bs

verbin = bs / 2. 

; Vertical histogram 

for i = 0, n_elements(obin) - 1 do begin
	plots, [ocount[i], ocount[i]], [obin[i] - verbin, obin[i] + verbin],  color = red, thick = th
	if i ne n_elements(obin)-1 then plots, [ocount[i],ocount[i+1]], [obin[i] + verbin, obin[i] + verbin], color = red, thick = th else $
		plots, [ocount[i],0], [obin[i] + verbin, obin[i] + verbin], color = red, thick = th 
	if i ne 0 then plots, [ocount[i],ocount[i-1]], [obin[i] - verbin, obin[i] - verbin], color = red, thick = th else $
		plots, [ocount[i],0], [obin[i] - verbin, obin[i] - verbin], color = red, thick = th 
endfor

for i = 0, n_elements(cbin) - 1 do begin
	plots, [ccount[i], ccount[i]], [cbin[i] - verbin, cbin[i] + verbin],  color = blue, thick = th
	if i ne n_elements(cbin)-1 then plots, [ccount[i],ccount[i+1]], [cbin[i] + verbin, cbin[i] + verbin], color = blue, thick = th else $
		plots, [ccount[i],0], [cbin[i] + verbin, cbin[i] + verbin], color = blue, thick = th 
	if i ne 0 then plots, [ccount[i],ccount[i-1]], [cbin[i] - verbin, cbin[i] - verbin], color = blue, thick = th else $
		plots, [ccount[i],0], [cbin[i] - verbin, cbin[i] - verbin], color = blue, thick = th 
endfor

plots, [x1,x1], [y0,y1], /normal, thick=th

plots, [x1, x2], [mean(dtemp), mean(dtemp)] * (y1-y0)/130. + y0, color=red, linestyle=1, /normal, thick=th
plots, [x1, x2], [mean(dtemp_con), mean(dtemp_con)] * (y1-y0)/130. + y0, color=blue, linestyle=1, /normal, thick=th

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

end
