
;+
; NAME:
;       
;	MOL_SENSITIVITY
;
; PURPOSE:
;
;	Assess the sensitivity to absorption in C2H2, HCN, and CO2 in the OHM spectra
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
;	IDL> .r mol_sensitivity
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Sep 08
;-

tags = ohmdat('tag')
ntag = n_elements(tags)

noise = dblarr(3, ntag)
cont = dblarr(3, ntag)

for i = 0, ntag - 1 do begin

	fname = tags[i]
	tag,fname,dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
	
	wave = sed.wave_sh
	flux = sed.flux_sh
	err = sed.err_sh
	
	pivarr = fillarr(0.5,13.0,16.0)
	npiv = n_elements(pivarr)
	piv = fltarr(npiv)
	for j=0,npiv - 1 do piv[j] = closetomed(wave,flux,pivarr[j])
	
	bind = closetomed(wave,flux,6.7)
	eind = closetomed(wave,flux,7.0)
	
	; Linear fit to continuum in three spots

	xarr = fillarr(1d-2,0,30)
	
	expr='p[0]+x*p[1]'
	start = [0,1]

	linwidth = 0.3

	bc2h2 = closeto(wave,13.7 - linwidth) & ec2h2 = closeto(wave,13.7 + linwidth)	
	bhcn = closeto(wave,14.02 - linwidth) & ehcn = closeto(wave,14.02 + linwidth)	
	bco2 = closeto(wave,15.0 - linwidth) & eco2 = closeto(wave,15.0 + linwidth)	

	fit_c2h2 = mpfitexpr(expr,wave[bc2h2:ec2h2], flux[bc2h2:ec2h2],0.1*flux[bc2h2:ec2h2],start,/quiet, perror = errc2h2)
	fitlin_c2h2 = fit_c2h2[0] + wave * fit_c2h2[1]
	fit_hcn = mpfitexpr(expr,wave[bhcn:ehcn], flux[bhcn:ehcn],0.1*flux[bhcn:ehcn],start,/quiet, perror = errhcn)
	fitlin_hcn = fit_hcn[0] + wave * fit_hcn[1]
	fit_co2 = mpfitexpr(expr,wave[bco2:eco2], flux[bco2:eco2],0.1*flux[bco2:eco2],start,/quiet, perror = errco2)
	fitlin_co2 = fit_co2[0] + wave * fit_co2[1]

	norm_c2h2 = flux - fitlin_c2h2
	norm_hcn =  flux - fitlin_hcn
	norm_co2 =  flux - fitlin_co2
	
	sigwidth = 0.3

	bc2h2 = closeto(wave,13.7 - sigwidth) & ec2h2 = closeto(wave,13.7 + sigwidth)	
	bhcn = closeto(wave,14.02 - sigwidth) & ehcn = closeto(wave,14.02 + sigwidth)	
	bco2 = closeto(wave,15.0 - sigwidth) & eco2 = closeto(wave,15.0 + sigwidth)	

	sig_c2h2 = stddev(norm_c2h2[bc2h2:ec2h2])
	sig_hcn = stddev(norm_hcn[bhcn:ehcn])
	sig_co2 = stddev(norm_co2[bco2:eco2])

	!p.multi = [0,1,3]

	plot, wave, norm_c2h2, $
		xr = [13.2,14.2], /xstyle, $
		psym = 10, $
		charsize = 1.5, $
		xtitle = 'Wavelength', $
		ytitle = 'C2H2'

	ver, 13.7, linestyle = 1

	plot, wave, norm_hcn, $
		xr = [13.52,14.52], /xstyle, $
		psym = 10, $
		charsize = 1.5, $
		xtitle = 'Wavelength', $
		ytitle = 'HCN'

	ver, 14.02, linestyle = 1

	plot, wave, norm_co2, $
		xr = [14.5,15.5], /xstyle, $
		psym = 10, $
		charsize = 1.5, $
		xtitle = 'Wavelength', $
		ytitle = 'CO2'

	ver, 15.0, linestyle = 1

;	print,''
;	print,'C2H2: ',sig_c2h2
;	print,'HCN: ',sig_hcn
;	print,'CO2: ',sig_co2
;	print,''

	cont[0,i] = flux(closetomed(wave,flux,13.7))
	cont[1,i] = flux(closetomed(wave,flux,14.02))
	cont[2,i] = flux(closetomed(wave,flux,15.0))

	noise[0,i] = sig_c2h2
	noise[1,i] = sig_hcn
	noise[2,i] = sig_co2

endfor

limfac = 3.		; N of detection (assume 3-sigma)

print,''
print, 'Average C2H2 continuum: ', string(mean(cont[0,*]), format = '(f5.3)'), ' +- ',string(stddev(cont[0,*]), format = '(f5.3)')
print, 'Average HCN  continuum: ', string(mean(cont[1,*]), format = '(f5.3)'), ' +- ',string(stddev(cont[1,*]), format = '(f5.3)')
print, 'Average CO2  continuum: ', string(mean(cont[2,*]), format = '(f5.3)'), ' +- ',string(stddev(cont[2,*]), format = '(f5.3)')
print,''
print, 'Average C2H2 S/N: ', string(mean(cont[0,*] / noise[0,*]), format = '(f6.3)'), ' +- ', $ 
	string(stddev(cont[0,*] / noise[0,*]), format = '(f6.3)')
print, 'Average HCN  S/N: ', string(mean(cont[1,*] / noise[1,*]), format = '(f6.3)'), ' +- ', $
	string(stddev(cont[1,*] / noise[1,*]), format = '(f6.3)')
print, 'Average CO2  S/N: ', string(mean(cont[2,*] / noise[2,*]), format = '(f6.3)'), ' +- ', $
	string(stddev(cont[2,*] / noise[2,*]), format = '(f6.3)')
print,''

end
