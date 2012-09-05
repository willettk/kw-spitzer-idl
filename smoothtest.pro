pro smoothtest

; Measure the S/N as a function of various smoothing routines performed in SMART

; KW, 21 Aug 2007


; NO SMOOTHING

fname1 = '~/Desktop/mega023_nosmooth.txt'

readcol, fname1, wave_ns, flux_ns, err_ns, flag_ns, format = 'f,f,f,i', /silent

pix1 = where(abs(wave_ns - 14.6) eq min(abs(wave_ns - 14.6)))
pix2 = where(abs(wave_ns - 15.2) eq min(abs(wave_ns - 15.2)))

plot, wave_ns, flux_ns, xr = [13,15]
oplot, wave_ns(pix1:pix2), flux_ns(pix1:pix2), color = fsc_color("Red")

sigma_ns = stddev(flux_ns(pix1:pix2))

; BOXCAR

fname2 = '~/Desktop/mega023_boxcar.txt'

readcol, fname2, wave_bc, flux_bc, err_bc, flag_bc, format = 'f,f,f,i', /silent

pix1 = where(abs(wave_bc - 14.6) eq min(abs(wave_bc - 14.6)))
pix2 = where(abs(wave_bc - 15.2) eq min(abs(wave_bc - 15.2)))

oplot, wave_bc, flux_bc, color = fsc_color("Yellow")
sigma_bc = stddev(flux_bc(pix1:pix2))

; HANNING

fname3 = '~/Desktop/mega023_hanning.txt'

readcol, fname3, wave_ha, flux_ha, err_ha, flag_ha, format = 'f,f,f,i', /silent

pix1 = where(abs(wave_ha - 14.6) eq min(abs(wave_ha - 14.6)))
pix2 = where(abs(wave_ha - 15.2) eq min(abs(wave_ha - 15.2)))

oplot, wave_ha, flux_ha, color = fsc_color("Blue")
sigma_ha = stddev(flux_ha(pix1:pix2))

; GAUSSIAN

fname4 = '~/Desktop/mega023_gaussian.txt'

readcol, fname4, wave_ga, flux_ga, err_ga, flag_ga, format = 'f,f,f,i', /silent

pix1 = where(abs(wave_ga - 14.6) eq min(abs(wave_ga - 14.6)))
pix2 = where(abs(wave_ga - 15.2) eq min(abs(wave_ga - 15.2)))

oplot, wave_ga, flux_ga, color = fsc_color("Green")
sigma_ga = stddev(flux_ga(pix1:pix2))

print,''
print, 'Non-smoothed 1-sigma [mJy]: ', sigma_ns * 1d3
print, 'Boxcar-smoothed 1-sigma: ', sigma_bc * 1d3
print, 'Hanning-smoothed 1-sigma: ', sigma_ha * 1d3
print, 'Gaussian-smoothed 1-sigma: ', sigma_ga * 1d3
print,''




stop
end
