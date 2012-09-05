; Figure 8 from Armus et al 2007 (adapted from Laurent et al 2000)


restore,'~/Astronomy/Research/Spitzer/OHM/lines/structures/pah62.sav'
restore,'~/Astronomy/Research/Spitzer/control/lines/structures/pah62_con.sav'

p62_con = pah62_con
p62 = pah62

fnames=ohmdat('tag')
cnames=condat('tag')
nf = n_elements(fnames)
nc = n_elements(cnames)
f5 = fltarr(nf)
f15 = fltarr(nf)
c5 = fltarr(nc)
c15 = fltarr(nc)

for i = 0,nf-1 do begin
	tag,fnames(i),dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fnames(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr

	bind5 = closeto(wave,5.3)
	eind5 = closeto(wave,5.8)
	addpah = fltarr(2,eind5 - bind5 + 1)
	for j = bind5, eind5 do begin
		addpah(0,j-bind5) = flux(j)
		addpah(1,j-bind5) = (wave(j+1) - wave(j)) * 3d14 / (wave(j))^2
	endfor
	f5(i) = total(addpah(0,*) * addpah(1,*)) * 1d-23

	bind15 = closeto(wave,14)
	eind15 = closeto(wave,16)
	addpah = fltarr(2,eind15 - bind15 + 1)
	for j = bind15, eind15 do begin
		addpah(0,j-bind15) = flux(j)
		addpah(1,j-bind15) = (wave(j+1) - wave(j)) * 3d14 / (wave(j))^2
	endfor
	f15(i) = total(addpah(0,*) * addpah(1,*)) * 1d-23

endfor

for i = 0,nc-1 do begin
	tag,cnames(i),dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	
	bind5 = closeto(wave,5.3)
	eind5 = closeto(wave,5.8)
	addpah = fltarr(2,eind5 - bind5 + 1)
	for j = bind5, eind5 do begin
		addpah(0,j-bind5) = flux(j)
		addpah(1,j-bind5) = (wave(j+1) - wave(j)) * 3d14 / (wave(j))^2
	endfor
	c5(i) = total(addpah(0,*) * addpah(1,*)) * 1d-23

	bind15 = closeto(wave,14)
	eind15 = closeto(wave,16)
	addpah = fltarr(2,eind15 - bind15 + 1)
	for j = bind15, eind15 do begin
		addpah(0,j-bind15) = flux(j)
		addpah(1,j-bind15) = (wave(j+1) - wave(j)) * 3d14 / (wave(j))^2
	endfor
	c15(i) = total(addpah(0,*) * addpah(1,*)) * 1d-23

endfor

p62_con = p62_con(where(p62_con gt 0))
c5 = c5(where(p62_con gt 0))
c15 = c15(where(p62_con gt 0))

psym=1

red = fsc_color("Red")
blue = fsc_color("Blue")

plotname='~/Astronomy/Comps2/figures/mir_diag.ps'
if psym eq 1 then begin
	set_plot,'ps'
	device, filename = plotname, /color,/landscape
	cs = 1
	cthick = 2
	lthick = 2
endif else begin
	cs = 2
	cthick = 1
	lthick = 1
endelse


!p.multi=[0,1,1]
plot, p62/f5, f15/f5, $
	/nodata, $
	xtitle='PAH 6.2 / f(5.5)', $
	ytitle='f(15.0)/f(5.5)', $
	title='MIR diagnostic from Laurent et al.', $
	charsize=cs,charthick=cthick,thick=lthick, $
	/xlog, /ylog, $
	xr=[0.005,11], /xstyle, $
	yr=[0.3,40], /ystyle

oplot, p62/f5,f15/f5, color=red, psym=symcat(14),thick=lthick
oplot, p62_con/c5,c15/c5, color=blue, psym=symcat(15),thick=lthick

xyouts,0.01,0.6,'AGN',/data,charsize=cs,charthick=cthick
xyouts,2,0.6,'PDR',/data,charsize=cs,charthick=cthick
xyouts,0.3,30,'HII',/data,charsize=cs,charthick=cthick

plots,[0.01,0.3],[0.6,30],linestyle=1,thick=lthick
plots,[0.01,2],[0.6,0.6],linestyle=1,thick=lthick

xyouts,0.15,20,/data,'50%',charsize=cs,charthick=cthick
xyouts,0.08,10,/data,'75%',charsize=cs,charthick=cthick
xyouts,0.03,3,/data,'90%',charsize=cs,charthick=cthick
xyouts,1.0,0.5,/data,'50%',charsize=cs,charthick=cthick
xyouts,0.7,0.5,/data,'75%',charsize=cs,charthick=cthick
xyouts,0.25,0.5,/data,'90%',charsize=cs,charthick=cthick

if psym eq 1 then begin
	device,/close
	set_plot,'x'
endif



stop
end
