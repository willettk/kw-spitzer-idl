;+
; NAME: 
;       LRAVG 
;
; PURPOSE:
;
;	Make a plot showing the scaled lo-res IRS spectra for both the OHMs and the control sample
;
; CATEGORY:
;	ASTRONOMY; GRAPHICS; ANALYSIS
;
; INPUTS:
;	
; OUTPUTS:
;
; KEYWORDS:
;
;	PS - 	create a PS file of the plot
;
;	IMAGE - save the ps file to my Comps 2 directory
;
; REQUIRES:
;
; EXAMPLE:
;
;	IDL> lravg
;         
; MODIFICATION HISTORY:
;
;	Written by KW - Aug 07
;	Added archived OHMs - Dec 07
;-
pro lravg_poster, ps=ps, image = image

; Load names of files from data structures

fnames = ohmdat('tag')
fnames = fnames(0:n_elements(fnames)-2)
nf = n_elements(fnames)

cnames = condat('tag')
cnames = cnames(0:n_elements(cnames)-1)
nc = n_elements(cnames)

anames = archdat('tag')
anames = anames(0:n_elements(anames)-1)
na = n_elements(anames)

; Wavelength grid over which to plot data

grid = fillarr(0.1,5.0,30.0)

allspec = fltarr(n_elements(grid),nf+na)
conspec = fltarr(n_elements(grid),nc)
archspec = fltarr(n_elements(grid),na)

!p.multi = [0,1,2]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Desktop/jila/lravg.ps', /color, /landscape
	cs = 1.5
	lthick = 3
	cthick = 3
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse


; Blank plot for OHMs


plot,indgen(50),indgen(50), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [1d-5,1], /ystyle, $
	/xlog, /ylog, $
;	xticks = 14, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30], $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'OHM sample', $
	thick = lthick, charthick = cthick, charsize = cs

; Loop over OHMs
	 
for i = 0, nf - 1 do begin

	; Restore data

	restore,'~/Astronomy/Research/Spitzer/OHM/data/structures/'+fnames(i)+'.sav'

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

	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	allspec(*,i) = newy
endfor

for i = 0, na - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/archived/data/structures/'+anames(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15
	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	allspec(*,i+nf) = newy
endfor

; Plot the median template from all spectra

meanohm = median(allspec,dim=2)
oplot,newx,meanohm,color=fsc_color("Red"),thick=4

	; Silicate and water ice absorption

	plots,[8,12],[4d-5,4d-5],linestyle=1, thick = lthick
	xyouts, 9, 2d-5,'Silicate', charsize = labsize, charthick = lthick
	plots,[17,20],[4d-3,4d-3], linestyle=1, thick = lthick
	xyouts, 17.5, 2d-3,'Silicate', charsize = labsize, charthick = lthick
	plots,[5.8,6.2],[5d-5,5d-5], linestyle=1, thick = lthick
	xyouts, 5.6, 2d-5,'H!I2!NO ice', charsize = labsize, charthick = lthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 7.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 12.6, 1.5d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 1.0d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 13.1, 1.3d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick
	
	lw = 18.7 & line_id = '[SIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 1.6d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 10.5 & line_id = '[SIV]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 1.5d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	ver, 15.0, linestyle = 2


; Control sample


plot,indgen(35), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [1d-4,1], /ystyle, $
;	xticks = 14, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30], $
	yticks = 5, ytickname=['10!E-4!N','10!E-3!N','10!E-2!N','10!E-1!N','10!E0!N','10!E1!N'], $
	/xlog, /ylog, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density', $
	title = 'Non-masing ULIRGs', $
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
	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	conspec(*,i) = newy
endfor

meancon = median(conspec,dim=2)
oplot,newx,meancon,color=fsc_color("Blue"),thick=4

	; Silicate and water ice absorption

	plots,[8,12],[4d-4,4d-4],linestyle=1, thick = lthick
	xyouts, 9, 2d-4,'Silicate', charsize = labsize, charthick = lthick
	plots,[17,20],[4d-3,4d-3], linestyle=1, thick = lthick
	xyouts, 17.5, 2d-3,'Silicate', charsize = labsize, charthick = lthick
	plots,[5.8,6.2],[4d-4,4d-4], linestyle=1, thick = lthick
	xyouts, 5.6, 2d-4,'H!I2!NO ice', charsize = labsize, charthick = lthick

	; PAH emission

	lw = 6.2 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 7.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 12.6, 1.5d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 1.8d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, 13.1, 1.5d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick
	
	lw = 18.7 & line_id = '[SIII]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 10.5 & line_id = '[SIV]' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N' & st = 5d-2
	plots, [lw,lw], [st, 8d-2], color=defcolor, linestyle = 0
	xyouts, lw, 2d-1, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick


	ver, 15.0, linestyle=2


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

stop
end
