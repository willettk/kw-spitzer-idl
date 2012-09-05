pro specavg, ps=ps, nolabels = nolabels, stop = stop
;+
; NAME:
;       
;	SPECAVG
;
; PURPOSE:
;
;	Create overlaid IRS spectra of all LR, SH, and LH modules for both the control sample and OH megamasers.
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
;	PS - create postscript hard copy of plots
;
;	NOLABELS - do not overlay lines indicating major emission/absorption features
;
; EXAMPLE:
;
;	IDL> specavg
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Mar 08
;	Added spectral index measurements for CSO plot - Oct 08
; 	Made paper plots B&W where necessary - Feb 10
;-

!p.font = -1

;######################
;#######  LR  #########
;######################

; Load names of files from data structures

onames = ohmdat('tag')
badohm = where(onames eq 'mega034')
goodohm = setdifference(indgen(n_elements(onames)),badohm)
onames = onames(goodohm)
no = n_elements(onames)

anames = transpose(archdat('tag'))
na = n_elements(anames)

cnames = condat('tag')
badcon = where(cnames eq 'control033')
goodcon = setdifference(indgen(n_elements(cnames)),badcon)
cnames = cnames(goodcon)
nc = n_elements(cnames)

; Wavelength grid over which to plot data

grid = fillarr(0.087,5.0,30.0)		; Avg. separation between wave bins is 0.087 um

allspec = fltarr(n_elements(grid),no+na)
conspec = fltarr(n_elements(grid),nc)

!p.multi = [0,1,2]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/lravg_ind_both.ps', /color, /landscape 
;	if keyword_set(image) then device,filename='~/Astronomy/Comps2/images/lravg.ps', /color, /landscape else $
;		device,filename='~/Astronomy/Comps2/figures/lravg.ps', /color, xsize=24, ysize=18, yoff=2
	cs = 1.5
	lthick = 5
	cthick = 5
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Dark Grey")
white = fsc_color("White")
black = fsc_color("Black")

; Blank plot for OHMs

plot,indgen(50),indgen(50), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [1d-5,1], /ystyle, $
	/xlog, /ylog, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'OHMs (Darling + archived)', $
	thick = lthick, charthick = cthick, charsize = cs

; Loop over OHMs
	 
allohms = [onames,anames]

for i = 0, n_elements(allohms) - 1 do begin

	; Restore data

	tag, allohms(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	; Overplot data

	oplot,newx,newy,color=grey, thick = lthick

	allspec[*,i] = newy
endfor

; Plot the median template from all spectra

meanohm = median(allspec,dim=2)
oplot,newx,meanohm,color=red,thick=2

	; Silicate and water ice absorption

	plots,[8,12],[4d-5,4d-5],linestyle=1, thick = lthick
	xyouts, 9, 2d-5,'Silicate', charsize = labelsize, charthick = lthick
	plots,[17,20],[4d-3,4d-3], linestyle=1, thick = lthick
	xyouts, 17.5, 2d-3,'Silicate', charsize = labelsize, charthick = lthick
	plots,[5.8,6.2],[4d-5,4d-5], linestyle=1, thick = lthick
	xyouts, 5.6, 2d-5,'H!I2!NO ice', charsize = labelsize, charthick = lthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 7.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 12.6, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 13.1, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 18.7 & line_id = '[SIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 10.5 & line_id = '[SIV]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	ver, 15.0, linestyle = 2


; Control sample


plot,indgen(35), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [1d-4,1], /ystyle, $
	/xlog, /ylog, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'Non-masing', $
	thick = lthick, charthick = cthick, charsize = cs

for i = 0, nc - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/control/data/structures/'+cnames(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15
	oplot,newx,newy,color=grey, thick = lthick

	conspec(*,i) = newy
endfor

meancon = median(conspec,dim=2)
oplot,newx,meancon,color=blue,thick=2

if not keyword_set(nolabels) then begin

	; Silicate and water ice absorption

	plots,[8,12],[4d-4,4d-4],linestyle=1, thick = lthick
	xyouts, 9, 2d-4,'Silicate', charsize = labelsize, charthick = lthick
	plots,[17,20],[4d-3,4d-3], linestyle=1, thick = lthick
	xyouts, 17.5, 2d-3,'Silicate', charsize = labelsize, charthick = lthick
	plots,[5.8,6.2],[4d-4,4d-4], linestyle=1, thick = lthick
	xyouts, 5.6, 2d-4,'H!I2!NO ice', charsize = labelsize, charthick = lthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 7.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 12.6, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 13.1, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 18.7 & line_id = '[SIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 10.5 & line_id = '[SIV]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

endif

	ver, 15.0, linestyle=2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Plot mean spectra without individual objects

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/lravg_both.ps', /color, /portrait 
	cs = 1.4
	lthick = 5
	cthick = 5
	medthick = 4
	defcolor=black
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
	medthick = 1
	defcolor=white
endelse

plot,indgen(35), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [5d-4,2d-1], /ystyle, $
	/xlog, /ylog, $
	xticks = 5, $
	xtickv = [5,10,15,20,25,30], $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, charthick = cthick, charsize = cs

oplot,newx,meanohm,color=defcolor,thick=lthick, psym=10
oplot,newx,meancon,color=grey,thick=medthick, psym=10

legend, /top, /left, ['OHMs ('+strtrim(n_elements(allohms),2)+')', 'non-masing ('+strtrim(nc,2)+')'], $
	color=[defcolor,grey], $
	linestyle=[0,0], $
	thick = [lthick,medthick], $
	charsize = 1, charthick=cthick

ver, 15.0, linestyle=2, thick=lthick

zz = 0
if zz eq 1 then begin

	; Silicate and water ice absorption

	plots,[8,12],[6d-4,6d-4],linestyle=1, thick = lthick
	xyouts, 9.5, 4d-4,'Silicate', charsize = labelsize, charthick = lthick
	plots,[17,20],[7d-3,7d-3], linestyle=1, thick = lthick
	xyouts, 17.5, 4d-3,'Silicate', charsize = labelsize, charthick = lthick
	plots,[5.8,6.2],[8d-4,8d-4], linestyle=1, thick = lthick
	xyouts, 5.6, 5d-4,'H!I2!NO ice', charsize = labelsize, charthick = lthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 1d-2
	plots, [lw,lw], [st, 5d-2], color=defcolor, linestyle = 0
	xyouts, lw, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 7.7 & line_id = 'PAH' & st = 1.5d-2
	plots, [lw,lw], [st, 5d-2], color=defcolor, linestyle = 0
	xyouts, lw, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 8.6 & line_id = 'PAH' & st = 1d-2
	plots, [lw,lw], [st, 5d-2], color=defcolor, linestyle = 0
	xyouts, lw, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 5d-2], color=defcolor, linestyle = 0
	xyouts, lw, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 5d-2], color=defcolor, linestyle = 0
	xyouts, 12.6, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.555 & line_id = '[NeIII]' & st = 2d-2
	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0
	xyouts, lw, 6d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 1.8d-2
	plots, [lw,lw], [st, 6d-2], color=defcolor, linestyle = 0
	xyouts, 13.3, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 18.7 & line_id = '[SIII]' & st = 3d-2
	plots, [lw,lw], [st, 5d-2], color=defcolor, linestyle = 0
	xyouts, lw, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 10.5 & line_id = '[SIV]' & st = 6d-3
	plots, [lw,lw], [st, 9d-3], color=defcolor, linestyle = 0
	xyouts, lw, 1.5d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N S(1)' & st = 2d-2
	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0
	xyouts, lw, 7d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N S(3)' & st = 4d-3
	plots, [lw,lw], [st, 8d-3], color=defcolor, linestyle = 0
	xyouts, lw, 1.3d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

endif

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

;stop

;######################
;#######  SH  #########
;######################

onames = ohmdat('tag')
no = n_elements(onames)

cnames = condat('tag')
nc = n_elements(cnames)

anames = archdat('tag')
na = n_elements(anames)

grid = fillarr(0.01,8,18)		; Avg. separation in wavelength bins is 0.0105 um

allspec = fltarr(n_elements(grid),no+na)
conspec = fltarr(n_elements(grid),nc)

!p.multi = [0,1,2]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/shavg_ind.ps', /color, /landscape
	cs = 1.5
	lthick = 2
	cthick = 2
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse


; OHMs

allohms_sh = [transpose(onames),transpose(anames)]

plot,indgen(50),indgen(50), /nodata, $
	xr = [7,19], /xstyle, $
	yr = [-0.01,0.1], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'OHMs (Darling + archived)', $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs

	 
for i = 0, no + na - 1 do begin

	tag, allohms_sh(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms_sh(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err =  sed.err_sh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15
	oplot,newx,newy,color=grey, thick = lthick

	allspec(*,i) = newy
endfor

meanohm_sh = median(allspec,dim=2)
oplot,newx,meanohm_sh,color=red,thick=2

if not keyword_set(nolabels) then begin

	; PAH emission

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, 12.7, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 6.5d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, 12.9, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 10.5 & line_id = '[SIV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N S(1)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N S(3)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.28 & line_id = 'H!I2!N S(2)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 14.322 & line_id = '[NeV]' & st = 4d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

endif

	ver, 15.0, linestyle = 2

; SH control sample

plot,indgen(50),indgen(50), /nodata, $
	xr = [7,19], /xstyle, $
	yr = [-0.01,0.1], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'Non-masing', $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs

for i = 0, nc - 1 do begin

	tag, cnames(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err =  sed.err_sh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15
	oplot,newx,newy,color=grey, thick = lthick

	allspec(*,i) = newy
endfor

meancon_sh = median(allspec,dim=2)
oplot,newx,meancon_sh,color=blue,thick=2

if not keyword_set(nolabels) then begin

	; PAH emission

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, 12.7, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 6.5d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, 12.9, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 10.5 & line_id = '[SIV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N S(1)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N S(3)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.28 & line_id = 'H!I2!N S(2)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 14.322 & line_id = '[NeV]' & st = 4d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

endif

	ver, 15.0, linestyle = 2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Plot SH mean spectra without individual objects

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/shavg.ps', /color, /landscape 
	cs = 1.5
	lthick = 2
	cthick = 2
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse

plot,indgen(50),indgen(50), /nodata, $
	xr = [7,19], /xstyle, $
	yr = [-0.005,0.03], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = 'Mean SH IRS spectra of OHMs vs. non-masing ULIRGs', $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs

oplot,newx,meanohm_sh,color=red,thick=lthick, psym=10
oplot,newx,meancon_sh,color=blue,thick=lthick, psym=10

legend, /top, /left, ['OHMs ('+strtrim(n_elements(allohms_sh),2)+')', 'non-masing ('+strtrim(nc,2)+')'], $
	color=[red, blue], $
	linestyle=[0,0], $
	thick = [2,2], $
	charsize = cs

if not keyword_set(nolabels) then begin

	; PAH emission

	lw = 8.6 & line_id = 'PAH' & st = 4d-3
	plots, [lw,lw], [st, 7d-3], color=defcolor, linestyle = 0
	xyouts, lw, 8d-3, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 9d-3
	plots, [lw,lw], [st, 0.012], color=defcolor, linestyle = 0
	xyouts, lw, 0.013, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 0.017
	plots, [lw,lw], [st, 0.022], color=defcolor, linestyle = 0
	xyouts, 12.7, 0.023, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.555 & line_id = '[NeIII]' & st = 0.014
	plots, [lw,lw], [st, 0.017], color=defcolor, linestyle = 0
	xyouts, lw, 0.018, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.814 & line_id = '[NeII]' & st = 0.020
	plots, [lw,lw], [st, 0.022], color=defcolor, linestyle = 0
	xyouts, 13.0, 0.023, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 10.511 & line_id = '[SIV]' & st = 0.005
	plots, [lw,lw], [st, 7d-3], color=defcolor, linestyle = 0
	xyouts, lw, 8d-3, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.04 & line_id = 'H!I2!N S(1)' & st = 0.015
	plots, [lw,lw], [st, 0.017], color=defcolor, linestyle = 0
	xyouts, lw, 0.018, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.28 & line_id = 'H!I2!N S(2)' & st = 9d-3
	plots, [lw,lw], [st, 0.011], color=defcolor, linestyle = 0
	xyouts, lw, 0.012, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.67 & line_id = 'H!I2!N S(3)' & st = 5d-3
	plots, [lw,lw], [st, 7d-3], color=defcolor, linestyle = 0
	xyouts, lw, 8d-3, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 14.322 & line_id = '[NeV]' & st = 0.011
	plots, [lw,lw], [st, 0.013], color=defcolor, linestyle = 0
	xyouts, lw, 0.014, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

endif

ver, 15.0, linestyle=2

;######################
;#######  LH  #########
;######################

onames = ohmdat('tag')
no = n_elements(onames)

cnames = condat('tag')
nc = n_elements(cnames)

anames = archdat('tag')
na = n_elements(anames)

grid = fillarr(0.02,16,33)		; Avg. separation in wavelength bins is 0.022 um

allspec = fltarr(n_elements(grid),no+na)
conspec = fltarr(n_elements(grid),nc)

!p.multi = [0,1,2]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/lhavg_ind.ps', /color, /landscape
	cs = 1.5
	lthick = 2
	cthick = 2
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse


; OHMs

allohms_lh = [transpose(onames),transpose(anames)]

plot,indgen(50),indgen(50), /nodata, $
	xr = [15,34], /xstyle, $
	yr = [-0.01,0.1], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'OHMs (Darling + archived)', $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs

	 
for i = 0, no + na - 1 do begin

	tag, allohms_lh(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms_lh(i)+'.sav'

	wave = sed.wave_lh
	flux = sed.flux_lh
	err =  sed.err_lh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy25 = newy(closeto(newx,25.0))

	newy = newy * 0.01 / newy25
	oplot,newx,newy,color=grey, thick = lthick

	allspec(*,i) = newy
endfor

meanohm_lh = median(allspec,dim=2)
oplot,newx,meanohm_lh,color=red,thick=2

if not keyword_set(nolabels) then begin

	; Fine-structure emission


	lw = 18.713 & line_id = '[SIII]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 24.318 & line_id = '[NeV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 28.22 & line_id = 'H!I2!N S(0)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.04 & line_id = 'H!I2!N S(1)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 25.890 & line_id = '[OIV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	ver, 25.0, linestyle = 2


endif

	ver, 25.0, linestyle = 2

; LH control sample

plot,indgen(50),indgen(50), /nodata, $
	xr = [15,34], /xstyle, $
	yr = [-0.01,0.1], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'Non-masing', $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs

	 
for i = 0, nc - 1 do begin

	tag, cnames(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames(i)+'.sav'

	wave = sed.wave_lh
	flux = sed.flux_lh
	err =  sed.err_lh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy25 = newy(closeto(newx,25.0))

	newy = newy * 0.01 / newy25
	oplot,newx,newy,color=grey, thick = lthick

	allspec(*,i) = newy
endfor

meancon_lh = median(allspec,dim=2)
oplot,newx,meancon_lh,color=blue,thick=2

if not keyword_set(nolabels) then begin

	; Fine-structure emission


	lw = 18.713 & line_id = '[SIII]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 24.318 & line_id = '[NeV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 28.22 & line_id = 'H!I2!N S(0)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.04 & line_id = 'H!I2!N S(1)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 25.890 & line_id = '[OIV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	ver, 25.0, linestyle = 2


	endif

	ver, 25.0, linestyle = 2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Plot LH mean spectra without individual objects

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/lhavg.ps', /color, /landscape 
	cs = 1.5
	lthick = 2
	cthick = 2
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse

plot,indgen(50),indgen(50), /nodata, $
	xr = [15,34], /xstyle, $
	yr = [-0.005,0.03], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = 'Mean LH IRS spectra of OHMs vs. non-masing ULIRGs', $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs

oplot,newx,meanohm_lh,color=red,thick=lthick, psym=10
oplot,newx,meancon_lh,color=blue,thick=lthick, psym=10

legend, /top, /left, ['OHMs ('+strtrim(n_elements(allohms_lh),2)+')', 'non-masing ('+strtrim(nc,2)+')'], $
	color=[red, blue], $
	linestyle=[0,0], $
	thick = [2,2], $
	charsize = cs

if not keyword_set(nolabels) then begin

	; Fine-structure emission


	lw = 18.713 & line_id = '[SIII]' & st = 0.006
	plots, [lw,lw], [st, 0.008], color=defcolor, linestyle = 0
	xyouts, lw, 0.009, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 24.318 & line_id = '[NeV]' & st = 0.010
	plots, [lw,lw], [st, 0.012], color=defcolor, linestyle = 0
	xyouts, lw, 0.013, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 28.22 & line_id = 'H!I2!N S(0)' & st = 0.019
	plots, [lw,lw], [st, 0.021], color=defcolor, linestyle = 0
	xyouts, lw, 0.022, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.04 & line_id = 'H!I2!N S(1)' & st = 0.005
	plots, [lw,lw], [st, 0.007], color=defcolor, linestyle = 0
	xyouts, lw, 0.008, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 25.890 & line_id = '[OIV]' & st = 0.014
	plots, [lw,lw], [st, 0.016], color=defcolor, linestyle = 0
	xyouts, 25.7, 0.017, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 25.988 & line_id = '[FeII]' & st = 0.015
	plots, [lw,lw], [st, 0.017], color=defcolor, linestyle = 0
	xyouts, 26.2, 0.0175, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 22.925 & line_id = '[FeII]' & st = 0.008
	plots, [lw,lw], [st, 0.009], color=defcolor, linestyle = 0
	xyouts, lw, 0.01, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	ver, 25.0, linestyle = 2


endif

ver, 25.0, linestyle=2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Individual LR plot just for the Darling OHMs

grid = fillarr(0.087,5.0,30.0)		; Avg. separation between wave bins is 0.087 um
ohmspec = fltarr(n_elements(grid),no)

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/lravg_ind.ps', /color, /portrait 
	cs = 1.5
	lthick = 2
	cthick = 2
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Grey")

labelsize=1.2

; Blank plot for OHMs

plot,indgen(50),indgen(50), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [1d-5,1], /ystyle, $
	xticks = 5, $
	xtickv = [5,10,15,20,25,30], $
	/xlog, /ylog, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Normalized flux density [Jy]', $
	thick = 3, charthick = 3, charsize = 1.5, $
	xthick = 3, $
	ythick = 3

; Loop over OHMs
	 
allohms = [onames]

for i = 0, n_elements(allohms) - 1 do begin

	; Restore data

	tag, allohms(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	; Overplot data

	oplot,newx,newy,color=grey, thick = lthick

	ohmspec(*,i) = newy
endfor

; Plot the median template from all spectra

meanohm = median(ohmspec,dim=2)
oplot,newx,meanohm,color=red,thick=3

	; Silicate and water ice absorption

	plots,[8,12],[4d-5,4d-5],linestyle=1, thick = 3
	xyouts, 9, 2d-5,'Silicate', charsize = labelsize, charthick = lthick
	plots,[17,20],[4d-3,4d-3], linestyle=1, thick = 3
	xyouts, 17.1, 2d-3,'Silicate', charsize = labelsize, charthick = lthick
	plots,[5.8,6.2],[4d-5,4d-5], linestyle=1, thick = 3
	xyouts, 5.6, 2d-5,'H!I2!NO ice', charsize = labelsize, charthick = lthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 7.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, 12.3, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, 13.1, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
	
	lw = 18.7 & line_id = '[SIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 10.5 & line_id = '[SIV]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0, thick=2
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick

	ver, 15.0, linestyle = 2, thick=2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; CSOs

cnames = csodat('tag')

	 ; Remove NGC 5793, 1245+676

	 csoind = [0,1,3,4,5,6,7,8]
	 cnames = cnames[csoind]

nc = n_elements(cnames)

grid = fillarr(0.081,5.0,34.0)		; Avg. separation between wave bins is 0.081 um
csospec = fltarr(n_elements(grid),nc)
csospecerr = fltarr(n_elements(grid),nc)

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/specavg_cso.ps', /color, /portrait 
	cs = 1.5
	lthick = 5
	cthick = 5
	defcolor=fsc_color("Black")
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
	defcolor=fsc_color("White")
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Grey")
black = fsc_color("Black")

labelsize=1.0

;; Blank plot for CSOs
;
;plot,indgen(50),indgen(50), /nodata, $
;	xr = [4.5,35], /xstyle, $
;	yr = [5d-4,3d-1], /ystyle, $
;	xticks = 6, $
;	xtickv = [5,10,15,20,25,30,35], $
;	/xlog, /ylog, $
;	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
;	ytitle = 'Normalized flux density [Jy]', $
;	thick = lthick, $
;	charthick = cthick, $
;	xthick = lthick, $
;	ythick = lthick, $
;	charsize = cs
;
;; Loop over CSOs
;	 
;for i = 0, nc - 1 do begin			; Remove cso010
;
;	; Restore data
;
;	tag, cnames(i), dirtag
;	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames(i)+'.sav'
;
;	wave = sed.wave_lr
;	flux = sed.flux_lr
;	err =  sed.err_lr
;
;	; Sample the lo-res spectra at 0.1 um intervals
;
;	newx = grid
;	newy    = fltarr(n_elements(newx))
;	newyerr = fltarr(n_elements(newx))
;	for j = 0,n_elements(newx)-1 do begin
;		newy(j) = flux(closeto(wave,newx[j]))
;		newyerr(j) = err(closeto(wave,newx[j]))
;	endfor
;
;	; Scale spectra to 10 mJy at 15.0 um
;
;	newy15 = newy(closeto(newx,15.0))
;	newy = newy * 0.01 / newy15
;	newyerr = newyerr * 0.01 / newy15
;
;	; Overplot data
;
;	badind = where(newy eq newy[n_elements(newy)-1], count)
;	if count gt 2 then begin
;		newy[badind[1]:badind[n_elements(badind)-1]] = !values.f_nan
;		newyerr[badind[1]:badind[n_elements(badind)-1]] = !values.f_nan
;	endif
;	oplot,newx,newy,color=grey, thick = lthick
;
;;	if cnames[i] eq 'cso009' then oplot,newx,newy,color=red, thick = lthick
;
;	csospec(*,i) = newy
;	csospecerr(*,i) = newyerr
;endfor
;
;; Plot the median template from all spectra
;
;meancso = median(csospec,dim=2)
;meancsoerr = median(csospecerr,dim=2)
;oplot,newx,meancso,thick=4, psym = 10, color=defcolor
;
;; Measure the spectral indices of the template
;
;	a1_min = 5.3
;	a1_max = 14.8
;	a2_min = 20.0
;	a2_max = 30.0
;
;	sp6  = cmw(newx, meancso, meancsoerr, a1_min)
;	sp15 = cmw(newx, meancso, meancsoerr, a1_max)
;	sp20 = cmw(newx, meancso, meancsoerr, a2_min)
;	sp30 = cmw(newx, meancso, meancsoerr, a2_max)
;	
;	index1 = alog10(sp15/sp6) / alog10(a1_max / a1_min)
;	index2 = alog10(sp30/sp20) / alog10(a2_max / a2_min)
;		
;	sig_f6  = stddev(meancso(closeto(newx,a1_min - 0.4):closeto(newx,a1_min + 0.4)))
;	sig_f15 = stddev(meancso(closeto(newx,a1_max - 0.8):closeto(newx,a1_max + 0.8)))
;	sig_f20 = stddev(meancso(closeto(newx,a2_min - 1.0):closeto(newx,a2_min + 1.0)))
;	sig_f30 = stddev(meancso(closeto(newx,a2_max - 1.0):closeto(newx,a2_max + 1.0)))
;
;	sig_wave6  = abs(newx[closeto(newx,a1_min)] - newx[closeto(newx,a1_min) + 1]) 
;	sig_wave15 = abs(newx[closeto(newx,a1_max)] - newx[closeto(newx,a1_max) + 1]) 
;	sig_wave20 = abs(newx[closeto(newx,a2_min)] - newx[closeto(newx,a2_min) + 1]) 
;	sig_wave30 = abs(newx[closeto(newx,a2_max)] - newx[closeto(newx,a2_max) + 1]) 
;
;	dalpha_df15    = 1d / (alog10(a1_max / a1_min) * sp15 * alog(10))
;	dalpha_df6     = -1d / (alog10(a1_max / a1_min) * sp6 * alog(10))
;	dalpha_dwave15 = -1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a1_max)
;	dalpha_dwave6  = 1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a2_max)
;
;	dalpha_df30    = 1d / (alog10(a2_max / a2_min) * sp30 * alog(10))
;	dalpha_df20     = -1d / (alog10(a2_max / a2_min) * sp20 * alog(10))
;	dalpha_dwave30 = -1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)
;	dalpha_dwave20  = 1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)
;
;	sig_alpha1 = sqrt(sig_f15^2 * dalpha_df15^2 + sig_f6^2  * dalpha_df6^2  + sig_wave15^2 * dalpha_dwave15^2 + sig_wave6^2  * dalpha_dwave6^2)
;	sig_alpha2 = sqrt(sig_f30^2 * dalpha_df30^2 + sig_f20^2 * dalpha_df20^2 + sig_wave30^2 * dalpha_dwave30^2 + sig_wave20^2 * dalpha_dwave20^2)
;
;	;print, 'NIR index: ',index1, ' +- ',sig_alpha1
;	;print, 'MIR index: ',index2, ' +- ',sig_alpha2
;
;	if not keyword_set(ps) then begin
;
;		oplot, [a1_min], [sp6], psym = symcat(14), symsize = 1.5, color=fsc_color("Red")
;		oplot, [a1_max],[sp15], psym = symcat(14), symsize = 1.5, color=fsc_color("Red")
;		oplot, [a2_min],[sp20], psym = symcat(14), symsize = 1.5, color=fsc_color("Yellow")
;		oplot, [a2_max],[sp30], psym = symcat(14), symsize = 1.5, color=fsc_color("Yellow")
;		
;		plots,[a1_min,a1_max],[sp6,sp15], color=fsc_color("Red")
;		plots,[a2_min,a2_max],[sp20,sp30], color=fsc_color("Yellow")
;		
;		xyouts,0.2,0.8,'!7a!3!I1!N = '+string(index1,format='(f5.2)'),color=fsc_color("Red"), /normal, charsize = 2
;		xyouts,0.2,0.7,'!7a!3!I2!N = '+string(index2,format='(f5.2)'),color=fsc_color("Yellow"), /normal, charsize = 2
;
;	endif
;
;	; Silicate and water ice absorption
;
;	plots,[8,12],[1.5d-3,1.5d-3],linestyle=1, thick = 3
;	xyouts, 9, 1d-3,'Silicate', charsize = labelsize, charthick = lthick
;	plots,[17,20],[4d-3,4d-3], linestyle=1, thick = 3
;	xyouts, 17.1, 2d-3,'Silicate', charsize = labelsize, charthick = lthick
;	plots,[5.8,6.2],[9d-4,9d-4], linestyle=1, thick = 3
;	xyouts, 5.6, 7d-4,'H!I2!NO ice', charsize = labelsize, charthick = lthick
;
;	; PAH emission
;
;	lw = 6.2 & line_id = 'PAH' & st = 2d-2
;	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 5d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 7.7 & line_id = 'PAH' & st = 2.5d-2
;	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 5d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 8.6 & line_id = 'PAH' & st = 1.5d-2
;	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 5d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 11.3 & line_id = 'PAH' & st = 2d-2
;	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 5d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 12.7 & line_id = 'PAH' & st = 2.5d-2
;	plots, [lw,lw], [st, 4d-2], color=defcolor, linestyle = 0, thick=2
;	xyouts, 12.4, 5d-2, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	; Fine-structure emission
;
;
;	lw = 25.9 & line_id = '[OIV]' & st = 7d-2
;	plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	xyouts, 25.9, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 15.6 & line_id = '[NeIII]' & st = 3d-2
;	plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	xyouts, 15.8, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 12.8 & line_id = '[NeII]' & st = 6d-2
;	plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	xyouts, 13.0, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;	
;	lw = 18.7 & line_id = '[SIII]' & st = 3d-2
;	plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	;lw = 10.5 & line_id = '[SIV]' & st = 9d-3
;	;plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	;xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 17.0 & line_id = 'H!I2!N' & st = 3.5d-2
;	plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	lw = 9.7 & line_id = 'H!I2!N' & st = 2d-2
;	plots, [lw,lw], [st, 1d-1], color=defcolor, linestyle = 0, thick=2
;	xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labelsize, /data, charthick = lthick
;
;	ver, 15.0, linestyle = 2, thick=2
;
;if keyword_set(ps) then begin
;	device,/close
;	set_plot,'x'
;endif
;
;stop

; Individual LR plot for all objects in the sample

grid = fillarr(0.087,5.0,30.0)		; Avg. separation between wave bins is 0.087 um

; Load data

onames = ohmdat('tag')
badohm = where(onames eq 'mega034')
goodohm = setdifference(indgen(n_elements(onames)),badohm)
onames = onames(goodohm)
no = n_elements(onames)

anames = transpose(archdat('tag'))
na = n_elements(anames)

cnames = condat('tag')
badcon = where(cnames eq 'control033')
goodcon = setdifference(indgen(n_elements(cnames)),badcon)
cnames = cnames(goodcon)
nc = n_elements(cnames)


!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/lravg_ind_all.ps', /color, /portrait 
	cs = 1.5
	lthick = 4
	cthick = 4
	defcolor = fsc_color("Black")
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
	defcolor = fsc_color("White")
endelse

labelsize=1.2

; Blank plot for OHMs

plot,indgen(50),indgen(50), /nodata, $
	xr = [4.8,30.5], /xstyle, $
	yr = [-4, 0], /ystyle, $
	xticks = 5, $
	xtickv = [5,10,15,20,25,30], $
;	ytickname = replicate(' ',6), $
	/xlog, $
;	/ylog, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'log (normalized S!I!7m!3!N)', $
	charthick = cthick, $
	charsize = 1.5, $
	thick = cthick, $
	xthick = cthick, $
	ythick = cthick

; Loop over OHMs
	 
allobjs = [onames,anames,cnames]
allspec = fltarr(n_elements(grid),n_elements(allobjs))
medianspec = fltarr(n_elements(grid))
sigmaspec = fltarr(n_elements(grid))

for i = 0, n_elements(allobjs) - 1 do begin

	; Restore data

	tag, allobjs(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allobjs(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	; Overplot data

;	oplot,newx,newy,color=grey, thick = lthick

	allspec[*,i] = newy
endfor

; Plot the 1-sigma envelope 

for i = 0, n_elements(grid)-1 do begin
	medianspec[i] = median(allspec[i,*])
	sigmaspec[i] = stddev(allspec[i,*])
endfor

possig = alog10(medianspec + sigmaspec)
negsig = alog10(medianspec - sigmaspec)
negsig[where(finite(negsig) ne 1)] = !y.crange[0]
negsig[where(negsig lt !y.crange[0])] = !y.crange[0]
;negsig = replicate(-4,n_elements(possig))
;negsig[where(negsig lt 10^(!y.crange[0]))] = 10^(!y.crange[0])

polyfill, [min(newx),newx,reverse(newx)], $
	[possig[0],possig,reverse(negsig)],$
	color=grey

; Plot the median template from all spectra

meanall = median(allspec,dim=2)
oplot,newx,alog10(meanall),color=defcolor,thick=4

	; Silicate and water ice absorption

	plots,[8,12],alog10([4d-4,4d-4]),linestyle=1, thick = lthick
	xyouts, 9, alog10(2d-4),'Silicate', charsize = labelsize, charthick = cthick
	plots,[17,20],alog10([4d-3,4d-3]), linestyle=1, thick = lthick
	xyouts, 17.1, alog10(2d-3),'Silicate', charsize = labelsize, charthick = cthick
	plots,[5.8,6.2],alog10([4d-4,4d-4]), linestyle=1, thick = lthick
	xyouts, 5.6, alog10(2d-4),'H!I2!NO ice', charsize = labelsize, charthick = cthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 7.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 11.3 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 12.7 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, 12.5, alog10(1.5d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[Ne III]' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, 15.8, alog10(1.5d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 12.8 & line_id = '[Ne II]' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, 13.3, alog10(1.5d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick
	
	lw = 18.7 & line_id = '[S III]' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 10.5 & line_id = '[S IV]' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 17.0 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

	lw = 9.7 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], alog10([st, 8d-2]), color=defcolor, linestyle = 0, thick=lthick
	xyouts, lw, alog10(2d-1), line_id, orientation = 90, charsize = labelsize, /data, charthick = cthick

;	ver, 15.0, linestyle = 2, thick=2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
