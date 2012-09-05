pro shavg, ps=ps

fnames = ohmdat('tag')
fnames = fnames(0:n_elements(fnames)-2)
nf = n_elements(fnames)

cnames = condat('tag')
cnames = cnames(0:n_elements(cnames)-1)
nc = n_elements(cnames)

anames = archdat('tag')
anames = anames(0:n_elements(anames)-1)
na = n_elements(anames)


grid = fillarr(0.1,8,18)

allspec = fltarr(n_elements(grid),nf+na)
conspec = fltarr(n_elements(grid),nc)

!p.multi = [0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
;	device,filename='~/Astronomy/Comps2/images/shavg.ps', /color, /landscape
	device,filename='~/Desktop/jila/shavg.ps', /color, /landscape
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
	xr = [7,19], /xstyle, $
	yr = [-0.02,0.1], /ystyle, $
;	/xlog, /ylog, $
;	xticks = 14, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30], $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = 'OHM averaged SH modules', $
	thick = lthick, charthick = cthick, charsize = cs

	 
for i = 0, nf - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/OHM/data/structures/'+fnames(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err =  sed.err_sh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15
	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	allspec(*,i) = newy
endfor

for i = 0, na - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/archived/data/structures/'+anames(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err =  sed.err_sh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15
	oplot,newx,newy,color=fsc_color("Grey"), thick = lthick

	allspec(*,i+nf) = newy
endfor

meanohm = median(allspec,dim=2)
oplot,newx,meanohm,color=fsc_color("Red"),thick=2

linelist=[6.2,7.7,8.6,11.3,12.7]
line_id = ['PAH','PAH','PAH','PAH','PAH']

	; PAH emission

	lw = 8.6 & line_id = 'PAH' & st = 2d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 11.3 & line_id = 'PAH' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.7 & line_id = 'PAH' & st = 4d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, 12.7, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	; Fine-structure emission


	lw = 15.6 & line_id = '[NeIII]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.8 & line_id = '[NeII]' & st = 6.5d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, 12.9, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick
	
	lw = 10.5 & line_id = '[SIV]' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 17.0 & line_id = 'H!I2!N S(1)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 9.7 & line_id = 'H!I2!N S(3)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 12.28 & line_id = 'H!I2!N S(2)' & st = 3d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	lw = 14.322 & line_id = '[NeV]' & st = 4d-2
	plots, [lw,lw], [st, 7d-2], color=defcolor, linestyle = 0
	xyouts, lw, 0.08, line_id, orientation = 90, charsize = labsize, /data, charthick = lthick

	ver, 15.0, linestyle = 2

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif


end

