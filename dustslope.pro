pro dustslope, ps = ps


; Plot the silicate strength vs. the spectral slopes; try to see where the dust lies (NIR or MIR)

fsil = ohmdat('sil')
csil = condat('sil')
asil = archdat('sil')

fslope = ohmdat('spindex')
cslope = condat('spindex')
aslope = archdat('spindex')

red = fsc_color("Red")
blue = fsc_color("Blue")

plotname='~/Astronomy/Research/Spitzer/papers/dustslope.ps'

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color, /portrait
	cs = 1.5
	cthick = 4
	lthick = 4
endif else begin
	cs = 2
	cthick = 1
	lthick = 1
endelse

plot, fsil(*,0),fslope(*,1), $
	/nodata, $
	xtitle = '9.7 !7l!3m S!Isil!N', $
	ytitle = '!7a!3!I30-20!N', $
	xrange = [0.5,-4], /xstyle, $
	yrange = [0,7.5], /ystyle, $
	charsize = cs, $
	charthick = cthick, $
	thick = lthick, $
	xthick = lthick, $
	ythick = lthick

oplot, fsil(0,*),fslope(1,*), psym = symcat(14), color = red, symsize=1.5
oplot, csil(0,*),cslope(1,*), psym = symcat(15), color = blue, symsize=1.5
oplot, asil(0,*),aslope(1,*), psym = symcat(14), color = red, symsize=1.5

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif


end
