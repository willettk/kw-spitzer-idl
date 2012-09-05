function cmw, wave, flux, err, var

; Find the flux at a given wavelength by taking a 1 um weighted average of the flux in that area

wstart = closeto(wave,var - 0.5)
wend = closeto(wave, var + 0.5)

top = total(flux(wstart:wend) * err(wstart:wend))
bottom = total(err(wstart:wend))

return,top/bottom

end

restore,'~/Astronomy/Research/Spitzer/control/conname.sav'
nf = n_elements(conname)

alpha1 = fltarr(nf)
alpha2 = fltarr(nf)
alphaerr=fltarr(nf,2)

for i = 0, nf - 1 do begin

	tag,conname(i),dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+conname(i)+'.sav'
	
	; Spectral index program
	
	sp6  = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, 5.3)
	sp15 = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, 14.8)
	sp20 = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, 20.0)
	sp30 = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, 30.0)
	
	plot,sed.wave_lr, sed.flux_lr, $
		/xlog, /ylog, $
		xrange = [4,32], /xstyle, $
		xtitle = 'Wavelength [um]', $
		ytitle = 'Flux [Jy]', $
		title = sed.obj
	
	oplot, [5.3], [sp6],  psym = symcat(14), symsize = 1.5, color=fsc_color("Red")
	oplot, [14.8],[sp15], psym = symcat(14), symsize = 1.5, color=fsc_color("Red")
	oplot, [20.0],[sp20], psym = symcat(14), symsize = 1.5, color=fsc_color("Yellow")
	oplot, [30.0],[sp30], psym = symcat(14), symsize = 1.5, color=fsc_color("Yellow")
	
	plots,[5.3,14.8],[sp6,sp15], color=fsc_color("Red")
	plots,[20.0,30.0],[sp20,sp30], color=fsc_color("Yellow")
	
	index1 = alog10(sp15/sp6) / alog10(14.8 / 5.3)
	index2 = alog10(sp30/sp20) / alog10(30.0 / 20.0)
	
	xyouts,0.2,0.8,'!7a!3!I1!N = '+string(index1,format='(f5.2)'),color=fsc_color("Red"), /normal, charsize = 2
	xyouts,0.2,0.7,'!7a!3!I2!N = '+string(index2,format='(f5.2)'),color=fsc_color("Yellow"), /normal, charsize = 2
	
	alpha1(i) = index1
	alpha2(i) = index2	

	d6  = stddev(sed.flux_lr(closeto(sed.wave_lr,4.9):closeto(sed.wave_lr,6)))
	d15 = stddev(sed.flux_lr(closeto(sed.wave_lr,14):closeto(sed.wave_lr,15.5)))
	d20 = stddev(sed.flux_lr(closeto(sed.wave_lr,19):closeto(sed.wave_lr,21)))
	d30 = stddev(sed.flux_lr(closeto(sed.wave_lr,29):closeto(sed.wave_lr,31)))

	err1 = sqrt(d15^2 + d6^2)
	err2 = sqrt(d30^2 + d20^2)
	alphaerr(i,*) = [err1,err2]

endfor

save,alpha1,alpha2,alphaerr,file='~/Astronomy/Research/Spitzer/'+dirtag+'/spindex_con.sav'

end

