pro lhavg, ps=ps

fnames = ohmdat('tag')
fnames = fnames(0:n_elements(fnames)-2)
nf = n_elements(fnames)

cnames = condat('tag')
cnames = cnames(0:n_elements(cnames)-1)
nc = n_elements(cnames)

anames = archdat('tag')
anames = anames(0:n_elements(anames)-1)
na = n_elements(anames)

grid = fillarr(0.1,17,32)

allspec = fltarr(n_elements(grid),nf+na)
conspec = fltarr(n_elements(grid),nc)

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Desktop/jila/lhavg.ps', /color, /landscape
	cs = 1.5
	lthick = 2
	cthick = 2
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse


; OHMs


plot,indgen(50),indgen(50), /nodata, $
	xr = [16,33], /xstyle, $
	yr = [-0.02,0.1], /ystyle, $
;	/xlog, /ylog, $
;	xticks = 14, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30], $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = 'OHM averaged LH modules', $
	thick = lthick, charthick = cthick, charsize = cs

	 
for i = 0, nf - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/OHM/data/structures/'+fnames(i)+'.sav'

	wave = sed.wave_lh
	flux = sed.flux_lh
	err =  sed.err_lh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy25 = newy(closeto(newx,25.0))

	newy = newy * 0.01 / newy25
	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	allspec(*,i) = newy
endfor

for i = 0, na - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/archived/data/structures/'+anames(i)+'.sav'

	wave = sed.wave_lh
	flux = sed.flux_lh
	err =  sed.err_lh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy25 = newy(closeto(newx,25.0))

	newy = newy * 0.01 / newy25
	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	allspec(*,i+nf) = newy
endfor

meanohm = median(allspec,dim=2)
oplot,newx,meanohm,color=fsc_color("Red"),thick=2

	; Fine-structure emission


	lw = 18.713 & line_id = '[SIII]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 24.318 & line_id = '[NeV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 28.22 & line_id = 'H!I2!N S(0)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 17.04 & line_id = 'H!I2!N S(1)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 25.890 & line_id = '[OIV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	ver, 25.0, linestyle = 2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif


end


