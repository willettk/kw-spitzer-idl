; Compute the S/N in the region surrounding posited HCN absorption at 14.0 um for a comparison with Lee's sample. 

; Easiest way is to find stddev and then compare to the depth quoted by Armus et al. 2007. 

fnames=ohmdat('tag')
nf=n_elements(fnames)
sigarr = fltarr(nf)
starr = fltarr(nf)
mflux = fltarr(nf)

for i = 0, nf - 1 do begin
	restore,'~/Astronomy/Research/Spitzer/OHM/data/structures/'+fnames(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err = sed.err_sh

	bind = closeto(wave,13.0)
	eind = closeto(wave,15.0)

	plot,wave,flux, $
		xr=[13.5,15.5],/xstyle, $
		yr=[0,0.2],/ystyle, $
		xtitle = 'Wavelength (rest frame) [!7l!3m]', $
		title = 'Optical depth in 14 !7l!3m region - '+sed.obj
	ver,13.7,color=fsc_color("Yellow")
	ver,14.02,color=fsc_color("Yellow")
	ver,15.0,color=fsc_color("Yellow")

	width = 0.2

	expr = 'p[0]+x*p[1]'
	start = [flux(bind),0.]

	fit = mpfitexpr(expr,wave(bind:eind),flux(bind:eind),err(bind:eind),start,/quiet)
	xarr=fillarr(1d-3,12,16)

	oplot,xarr,fit(0)+fit(1)*xarr,color=fsc_color("Red")

	resid = flux - (fit(0) + fit(1)*wave)

	sigma = stddev(resid(bind:eind))
	
	xyouts,14.5,0.05,'!7r!3 = '+string(sigma,format='(f7.5)'),/data,charsize=2
	sigarr(i) = sigma
	mflux(i) = flux(closetomed(wave,flux,14.0))

	aa=1
	if aa eq 1 then begin
	; Normalize and plot in optical depth

	tau = -1d*alog(resid+1)
	sigtau = -1d*alog(sigma+1)
	starr(i) = sigtau
	plot,wave,resid+1, $
		yr=[0.5,1.1],/ystyle, $
		xr=[13.5,15.5],/xstyle, $
		xtitle = 'Wavelength (rest frame) [!7l!3m]', $
		ytitle = '!7s!3', $
		title = 'Optical depth in 14 !7l!3m region - '+sed.obj
	ver,13.7,color=fsc_color("Yellow")
	ver,14.02,color=fsc_color("Yellow")
	ver,15.0,color=fsc_color("Yellow")
	hor,0.7,color=fsc_color("Green")
	hor,0.95,color=fsc_color("Green")

	endif
	junk=get_kbrd()

endfor




stop
end
