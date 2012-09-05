
;+
; NAME:
;       
;	NEV_COMP
;
; PURPOSE:
;
;	Compare OHMs, non-masing galaxies depending on if they show Ne V or not
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
;	IDL> 
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Aug 2010
;-

!p.multi=[0,2,1]

; Get names of galaxies

ohms_nev14 = getlinelist('NeV')
ohms_nev14 = transpose(ohms_nev14[0,*])
ohms_nev24 = getlinelist('NeV24')
ohms_nev24 = transpose(ohms_nev24[0,*])

ohms_nev = cmset_op(ohms_nev14, 'or', ohms_nev24)

allmegatag = ohmdat('tag')
allarchtag = archdat('tag')
alltag = [transpose(allmegatag),transpose(allarchtag)]

ohms_no_nev = cmset_op(alltag, 'and', /not2, ohms_nev)

grid = fillarr(0.087,5.0,30.0)		; Avg. separation between wave bins is 0.087 um
spec_ohms_nev = fltarr(n_elements(grid), n_elements(ohms_nev))
spec_ohms_no_nev = fltarr(n_elements(grid), n_elements(ohms_no_nev))

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Dark Grey")
white = fsc_color("White")
black = fsc_color("Black")

for i=0, n_elements(ohms_nev) - 1 do begin

	tag, ohms_nev[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+ohms_nev[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy[j] = flux(closeto(wave,newx[j]))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	spec_ohms_nev[*,i] = newy

endfor

; Plot the median template 

meanohm = median(spec_ohms_nev,dim=2)
plot,newx,spec_ohms_nev,thick=2, title='OHMs', /xlog, /ylog, psym=10, /xstyle, /ystyle, yr=[1d-4,5d-1], xr=[4,31]
;ver, 15.0, linestyle=2

for i=0, n_elements(ohms_no_nev) - 1 do begin

	tag, ohms_no_nev[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+ohms_no_nev[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy[j] = flux(closeto(wave,newx[j]))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	spec_ohms_no_nev[*,i] = newy

endfor

; Plot the median template from all spectra

meanohm = median(spec_ohms_no_nev,dim=2)
oplot,newx,spec_ohms_no_nev,thick=2, color=red, psym=10

xyouts, /data, charsize=1.5, 5, 1d-1, 'White: OHMs with NeV '+string(n_elements(ohms_nev),format='(i2)')
xyouts, /data, charsize=1.5, 5, 6d-2, 'Red:  OHMs w/o  NeV '+string(n_elements(ohms_no_nev),format='(i2)'), color=red
ver, 14.4, linestyle=1 & ver, 24.4, linestyle=1

; Measure average 30-20 um slope for both samples

	a1_min = 5.3
	a1_max = 14.8
	a2_min = 20.0
	if max(newx) lt 30. then a2_max = 28.0 else a2_max = 30.0

	; Spectral index program
	
	 sp6_ind = closeto(newx, a1_min) & sp6 = spec_ohms_nev[sp6_ind]
	sp15_ind = closeto(newx, a1_max) & sp15 = spec_ohms_nev[sp15_ind]
	sp20_ind = closeto(newx, a2_min) & sp20 = spec_ohms_nev[sp20_ind]
	sp30_ind = closeto(newx, a2_max) & sp30 = spec_ohms_nev[sp30_ind]
	
	; Calculate spectral indices

	index1 = alog10(sp15/sp6) / alog10(a1_max / a1_min)
	index2 = alog10(sp30/sp20) / alog10(a2_max / a2_min)
		
	; Find errors in flux measurements

	sig_f6  = stddev(spec_ohms_nev(closeto(newx,a1_min - 0.4):closeto(newx,a1_min + 0.4)))
	sig_f15 = stddev(spec_ohms_nev(closeto(newx,a1_max - 0.8):closeto(newx,a1_max + 0.8)))
	sig_f20 = stddev(spec_ohms_nev(closeto(newx,a2_min - 1.0):closeto(newx,a2_min + 1.0)))
	sig_f30 = stddev(spec_ohms_nev(closeto(newx,a2_max - 1.0):closeto(newx,a2_max + 1.0)))

	; Calculate errors in slope

	sig_wave6  = abs(newx[closeto(newx,a1_min)] - newx[closeto(newx,a1_min) + 1]) 
	sig_wave15 = abs(newx[closeto(newx,a1_max)] - newx[closeto(newx,a1_max) + 1]) 
	sig_wave20 = abs(newx[closeto(newx,a2_min)] - newx[closeto(newx,a2_min) + 1]) 
	sig_wave30 = abs(newx[closeto(newx,a2_max)] - newx[closeto(newx,a2_max) + 1]) 

	dalpha_df15    = 1d / (alog10(a1_max / a1_min) * sp15 * alog(10))
	dalpha_df6     = -1d / (alog10(a1_max / a1_min) * sp6 * alog(10))
	dalpha_dwave15 = -1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a1_max)
	dalpha_dwave6  = 1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a2_max)

	dalpha_df30    = 1d / (alog10(a2_max / a2_min) * sp30 * alog(10))
	dalpha_df20     = -1d / (alog10(a2_max / a2_min) * sp20 * alog(10))
	dalpha_dwave30 = -1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)
	dalpha_dwave20  = 1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)

	; Propagate errors from calculation of indices

	sig_alpha1 = sqrt(sig_f15^2 * dalpha_df15^2 + sig_f6^2  * dalpha_df6^2  + sig_wave15^2 * dalpha_dwave15^2 + sig_wave6^2  * dalpha_dwave6^2)
	sig_alpha2 = sqrt(sig_f30^2 * dalpha_df30^2 + sig_f20^2 * dalpha_df20^2 + sig_wave30^2 * dalpha_dwave30^2 + sig_wave20^2 * dalpha_dwave20^2)

	print,''
	print, '15-6  um slope (   NeV OHMs): ', string(index1,format='(f5.2)'), ' +- ', string(sig_alpha1,format='(f5.2)')
	print, '30-20 um slope (   NeV OHMs): ', string(index2,format='(f5.2)'), ' +- ', string(sig_alpha2,format='(f5.2)')

	; Spectral index program
	
	 sp6  = spec_ohms_no_nev[sp6_ind]
	 sp15 = spec_ohms_no_nev[sp15_ind]
	 sp20 = spec_ohms_no_nev[sp20_ind]
	 sp30 = spec_ohms_no_nev[sp30_ind]
	
	; Calculate spectral indices

	index1 = alog10(sp15/sp6) / alog10(a1_max / a1_min)
	index2 = alog10(sp30/sp20) / alog10(a2_max / a2_min)
		
	; Find errors in flux measurements

	sig_f6  = stddev(spec_ohms_nev(closeto(newx,a1_min - 0.4):closeto(newx,a1_min + 0.4)))
	sig_f15 = stddev(spec_ohms_nev(closeto(newx,a1_max - 0.8):closeto(newx,a1_max + 0.8)))
	sig_f20 = stddev(spec_ohms_nev(closeto(newx,a2_min - 1.0):closeto(newx,a2_min + 1.0)))
	sig_f30 = stddev(spec_ohms_nev(closeto(newx,a2_max - 1.0):closeto(newx,a2_max + 1.0)))

	; Calculate errors in slope

	sig_wave6  = abs(newx[closeto(newx,a1_min)] - newx[closeto(newx,a1_min) + 1]) 
	sig_wave15 = abs(newx[closeto(newx,a1_max)] - newx[closeto(newx,a1_max) + 1]) 
	sig_wave20 = abs(newx[closeto(newx,a2_min)] - newx[closeto(newx,a2_min) + 1]) 
	sig_wave30 = abs(newx[closeto(newx,a2_max)] - newx[closeto(newx,a2_max) + 1]) 

	dalpha_df15    = 1d / (alog10(a1_max / a1_min) * sp15 * alog(10))
	dalpha_df6     = -1d / (alog10(a1_max / a1_min) * sp6 * alog(10))
	dalpha_dwave15 = -1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a1_max)
	dalpha_dwave6  = 1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a2_max)

	dalpha_df30    = 1d / (alog10(a2_max / a2_min) * sp30 * alog(10))
	dalpha_df20     = -1d / (alog10(a2_max / a2_min) * sp20 * alog(10))
	dalpha_dwave30 = -1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)
	dalpha_dwave20  = 1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)

	; Propagate errors from calculation of indices

	sig_alpha1 = sqrt(sig_f15^2 * dalpha_df15^2 + sig_f6^2  * dalpha_df6^2  + sig_wave15^2 * dalpha_dwave15^2 + sig_wave6^2  * dalpha_dwave6^2)
	sig_alpha2 = sqrt(sig_f30^2 * dalpha_df30^2 + sig_f20^2 * dalpha_df20^2 + sig_wave30^2 * dalpha_dwave30^2 + sig_wave20^2 * dalpha_dwave20^2)

	print, '15-6  um slope (no NeV OHMs): ', string(index1,format='(f5.2)'), ' +- ', string(sig_alpha1,format='(f5.2)')
	print, '30-20 um slope (no NeV OHMs): ', string(index2,format='(f5.2)'), ' +- ', string(sig_alpha2,format='(f5.2)')

; KS tests for OH properties of NeV vs. non-NeV OHMs

	; L_OH
	
	aloh = archdat('logoh',/ver)
	oloh = ohmdat('logoh',/ver)
	logoh_names = [transpose(aloh[0,*]),transpose(oloh[0,*])]
	logoh = [transpose(aloh[1,*]),transpose(oloh[1,*])]
	match, logoh_names, ohms_nev, anev, nev_oh
	match, logoh_names, ohms_no_nev, anonev, no_nev_oh
	
	logoh_nev = float(logoh[anev])
	logoh_nonev = float(logoh[anonev])
	
	kstwo, logoh_nev, logoh_nonev, D_nir, prob_nir
	gauss_nir = sqrt(2d) * inverf(1d - prob_nir)
	
	print,''
	print,'D_KS    for log L_OH: '+string(D_nir,format='(f7.3)')
	print,'KS-prob for log L_OH: '+string(prob_nir,format='(f7.3)')
	print,'Gaussian probability:  '+string(gauss_nir,format='(f7.3)')+' sig'
	print,'Average value (Ne V) :  '+string(mean(logoh_nev),format='(f7.2)')+' +- '+string(stddev(logoh_nev),format='(f7.2)')
	print,'Average value (no Ne V) :  '+string(mean(logoh_nonev),format='(f7.2)')+' +- '+string(stddev(logoh_nonev),format='(f7.2)')

	; F1667
	
	aloh = archdat('f1667',/ver)
	oloh = ohmdat('f1667',/ver)
	f1667_names = [transpose(aloh[0,*]),transpose(oloh[0,*])]
	f1667 = [transpose(aloh[1,*]),transpose(oloh[1,*])]
	match, f1667_names, ohms_nev, anev, nev_oh
	match, f1667_names, ohms_no_nev, anonev, no_nev_oh
	
	f1667_nev = float(f1667[anev])
	f1667_nonev = float(f1667[anonev])
	
	kstwo, f1667_nev, f1667_nonev, D_nir, prob_nir
	gauss_nir = sqrt(2d) * inverf(1d - prob_nir)
	
	print,''
	print,'D_KS    for 1667 MHz flux: '+string(D_nir,format='(f7.3)')
	print,'KS-prob for 1667 MHz flux: '+string(prob_nir,format='(f7.3)')
	print,'Gaussian probability:  '+string(gauss_nir,format='(f7.3)')+' sig'
	print,'Average value (Ne V) :  '+string(mean(f1667_nev),format='(f7.2)')+' +- '+string(stddev(f1667_nev),format='(f7.2)')
	print,'Average value (no Ne V) :  '+string(mean(f1667_nonev),format='(f7.2)')+' +- '+string(stddev(f1667_nonev),format='(f7.2)')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get names of non-masing galaxies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

con_nev14 = getlinelist('NeV',/con)
con_nev14 = transpose(con_nev14[0,*])
con_nev24 = getlinelist('NeV24',/con)
con_nev24 = transpose(con_nev24[0,*])

con_nev = cmset_op(con_nev14, 'or', con_nev24)

allcontag = transpose(condat('tag'))

con_no_nev = cmset_op(allcontag, 'and', /not2, con_nev)

grid = fillarr(0.087,5.0,30.0)		; Avg. separation between wave bins is 0.087 um
spec_con_nev = fltarr(n_elements(grid), n_elements(con_nev))
spec_con_no_nev = fltarr(n_elements(grid), n_elements(con_no_nev))

for i=0, n_elements(con_nev) - 1 do begin

	tag, con_nev[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+con_nev[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy[j] = flux(closeto(wave,newx[j]))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	spec_con_nev[*,i] = newy

endfor

; Plot the median template 

meanohm = median(spec_con_nev,dim=2)
plot,newx,spec_con_nev,thick=2, title='Non-masing', /xlog, /ylog, psym=10, /xstyle, /ystyle, yr=[1d-4,5d-1], xr=[4,31]
;ver, 15.0, linestyle=2

for i=0, n_elements(con_no_nev) - 1 do begin

	tag, con_no_nev[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+con_no_nev[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy[j] = flux(closeto(wave,newx[j]))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,15.0))
	newy = newy * 0.01 / newy15

	spec_con_no_nev[*,i] = newy

endfor

; Plot the median template from all spectra

meanohm = median(spec_con_no_nev,dim=2)
oplot,newx,spec_con_no_nev,thick=2, color=red, psym=10

xyouts, /data, charsize=1.5, 5, 1d-1, 'White: Non-masing with NeV '+string(n_elements(con_nev),format='(i2)')
xyouts, /data, charsize=1.5, 5, 6d-2, 'Red:  Non-masing w/o  NeV '+string(n_elements(con_no_nev),format='(i2)'), color=red
ver, 14.4, linestyle=1 & ver, 24.4, linestyle=1

	print,''

; Measure average 30-20 um slope for both samples

	a1_min = 5.3
	a1_max = 14.8
	a2_min = 20.0
	if max(newx) lt 30. then a2_max = 28.0 else a2_max = 30.0

	; Spectral index program
	
	 sp6_ind = closeto(newx, a1_min) & sp6 = spec_con_nev[sp6_ind]
	sp15_ind = closeto(newx, a1_max) & sp15 = spec_con_nev[sp15_ind]
	sp20_ind = closeto(newx, a2_min) & sp20 = spec_con_nev[sp20_ind]
	sp30_ind = closeto(newx, a2_max) & sp30 = spec_con_nev[sp30_ind]
	
	; Calculate spectral indices

	index1 = alog10(sp15/sp6) / alog10(a1_max / a1_min)
	index2 = alog10(sp30/sp20) / alog10(a2_max / a2_min)
		
	; Find errors in flux measurements

	sig_f6  = stddev(spec_con_nev(closeto(newx,a1_min - 0.4):closeto(newx,a1_min + 0.4)))
	sig_f15 = stddev(spec_con_nev(closeto(newx,a1_max - 0.8):closeto(newx,a1_max + 0.8)))
	sig_f20 = stddev(spec_con_nev(closeto(newx,a2_min - 1.0):closeto(newx,a2_min + 1.0)))
	sig_f30 = stddev(spec_con_nev(closeto(newx,a2_max - 1.0):closeto(newx,a2_max + 1.0)))

	; Calculate errors in slope

	sig_wave6  = abs(newx[closeto(newx,a1_min)] - newx[closeto(newx,a1_min) + 1]) 
	sig_wave15 = abs(newx[closeto(newx,a1_max)] - newx[closeto(newx,a1_max) + 1]) 
	sig_wave20 = abs(newx[closeto(newx,a2_min)] - newx[closeto(newx,a2_min) + 1]) 
	sig_wave30 = abs(newx[closeto(newx,a2_max)] - newx[closeto(newx,a2_max) + 1]) 

	dalpha_df15    = 1d / (alog10(a1_max / a1_min) * sp15 * alog(10))
	dalpha_df6     = -1d / (alog10(a1_max / a1_min) * sp6 * alog(10))
	dalpha_dwave15 = -1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a1_max)
	dalpha_dwave6  = 1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a2_max)

	dalpha_df30    = 1d / (alog10(a2_max / a2_min) * sp30 * alog(10))
	dalpha_df20     = -1d / (alog10(a2_max / a2_min) * sp20 * alog(10))
	dalpha_dwave30 = -1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)
	dalpha_dwave20  = 1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)

	; Propagate errors from calculation of indices

	sig_alpha1 = sqrt(sig_f15^2 * dalpha_df15^2 + sig_f6^2  * dalpha_df6^2  + sig_wave15^2 * dalpha_dwave15^2 + sig_wave6^2  * dalpha_dwave6^2)
	sig_alpha2 = sqrt(sig_f30^2 * dalpha_df30^2 + sig_f20^2 * dalpha_df20^2 + sig_wave30^2 * dalpha_dwave30^2 + sig_wave20^2 * dalpha_dwave20^2)

	print,''
	print, '15-6  um slope (   NeV non-masing): ', string(index1,format='(f5.2)'), ' +- ', string(sig_alpha1,format='(f5.2)')
	print, '30-20 um slope (   NeV non-masing): ', string(index2,format='(f5.2)'), ' +- ', string(sig_alpha2,format='(f5.2)')

	; Spectral index program
	
	 sp6  = spec_con_no_nev[sp6_ind]
	 sp15 = spec_con_no_nev[sp15_ind]
	 sp20 = spec_con_no_nev[sp20_ind]
	 sp30 = spec_con_no_nev[sp30_ind]
	
	; Calculate spectral indices

	index1 = alog10(sp15/sp6) / alog10(a1_max / a1_min)
	index2 = alog10(sp30/sp20) / alog10(a2_max / a2_min)
		
	; Find errors in flux measurements

	sig_f6  = stddev(spec_con_no_nev(closeto(newx,a1_min - 0.4):closeto(newx,a1_min + 0.4)))
	sig_f15 = stddev(spec_con_no_nev(closeto(newx,a1_max - 0.8):closeto(newx,a1_max + 0.8)))
	sig_f20 = stddev(spec_con_no_nev(closeto(newx,a2_min - 1.0):closeto(newx,a2_min + 1.0)))
	sig_f30 = stddev(spec_con_no_nev(closeto(newx,a2_max - 1.0):closeto(newx,a2_max + 1.0)))

	; Calculate errors in slope

	sig_wave6  = abs(newx[closeto(newx,a1_min)] - newx[closeto(newx,a1_min) + 1]) 
	sig_wave15 = abs(newx[closeto(newx,a1_max)] - newx[closeto(newx,a1_max) + 1]) 
	sig_wave20 = abs(newx[closeto(newx,a2_min)] - newx[closeto(newx,a2_min) + 1]) 
	sig_wave30 = abs(newx[closeto(newx,a2_max)] - newx[closeto(newx,a2_max) + 1]) 

	dalpha_df15    = 1d / (alog10(a1_max / a1_min) * sp15 * alog(10))
	dalpha_df6     = -1d / (alog10(a1_max / a1_min) * sp6 * alog(10))
	dalpha_dwave15 = -1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a1_max)
	dalpha_dwave6  = 1d * (alog10(sp15/sp6) * alog(10)) / ((alog10(a1_max / a1_min))^2 * a2_max)

	dalpha_df30    = 1d / (alog10(a2_max / a2_min) * sp30 * alog(10))
	dalpha_df20     = -1d / (alog10(a2_max / a2_min) * sp20 * alog(10))
	dalpha_dwave30 = -1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)
	dalpha_dwave20  = 1d * (alog10(sp30/sp20) * alog(10)) / ((alog10(a2_max / a2_min))^2 * a2_max)

	; Propagate errors from calculation of indices

	sig_alpha1 = sqrt(sig_f15^2 * dalpha_df15^2 + sig_f6^2  * dalpha_df6^2  + sig_wave15^2 * dalpha_dwave15^2 + sig_wave6^2  * dalpha_dwave6^2)
	sig_alpha2 = sqrt(sig_f30^2 * dalpha_df30^2 + sig_f20^2 * dalpha_df20^2 + sig_wave30^2 * dalpha_dwave30^2 + sig_wave20^2 * dalpha_dwave20^2)

	print, '15-6  um slope (no NeV non-masing): ', string(index1,format='(f5.2)'), ' +- ', string(sig_alpha1,format='(f5.2)')
	print, '30-20 um slope (no NeV non-masing): ', string(index2,format='(f5.2)'), ' +- ', string(sig_alpha2,format='(f5.2)')

; Make figure for Paper II (referee and Jeremy suggestion)

set_plot,'ps'
device,filename='~/Astronomy/Research/Spitzer/papers/nev.ps', /color, decomposed=1

thick = 5
cthick = 5
csize = 1.6

nonev_color = fsc_color("Red")

plot,newx, alog10(spec_ohms_nev),$
	thick=thick, xthick = thick, ythick = thick, $ 
	charthick = cthick, $
	charsize = csize, $
	psym=10, $
	xstyle=5, $
	/ystyle, $
	/xlog, $
	yr=[-3.5,-0.5], $
	xr=[4,32], $
	xticks = 8, xtickv=[5,6,7,8,9,10,15,20,30], xtickname=replicate(' ',9), $
	yticks = 2, ytickv=[-3,-2,-1], $
	position=[0.13, 0.5, 0.90, 0.95]

axis, xaxis = 0, xstyle=1, xthick=thick, xticks = 8, xtickv=[5,6,7,8,9,10,15,20,30], xtickname=replicate(' ',9)
axis, xaxis = 1, xstyle=1, xthick=thick, xticks = 8, xtickv=[5,6,7,8,9,10,15,20,30], xtickname=replicate(' ',9)

oplot,newx,alog10(spec_ohms_no_nev),thick=thick, color=nonev_color, psym=10

legend, /top, /left, ['[Ne V]','no [Ne V]'], charthick=4, thick=[5,5], linestyle=[0,0], color=[fsc_color("Black"),nonev_color], number=0.5, charsize = 1.2

xyouts, /data, 20, -3, 'OHMs', charthick = 4, charsize = 1.5
ver, 15, linestyle=1, thick=3


plot,newx, alog10(spec_con_nev),$
	thick=thick, xthick = thick, ythick = thick, $ 
	charthick = cthick, $
	charsize = csize, $
	psym=10, $
	/ystyle, $
	xstyle = 9, $
	/xlog, $
	yr=[-4,-1], $
	xr=[4,32], $
	xticks = 8, xtickv=[5,6,7,8,9,10,15,20,30], $
	yticks = 3, ytickv=[-4,-3,-2,-1], $
	xtitle='Wavelength [!7l!3m]', $
	position=[0.13, 0.12, 0.90, 0.5]

oplot,newx,alog10(spec_con_no_nev),thick=thick, color=nonev_color, psym=10

xyouts, /data, 17, -3.5, 'Non-masing', charthick = 4, charsize = 1.5
ver, 15, linestyle=1, thick=3

	; Add Y title

	xyouts, orientation=90, /normal, 0.05, 0.3, 'log (normalized S!I!7m!3!N)', charthick=4, charsize=csize

device, /close
set_plot,'x'

save, filename='~/Astronomy/Research/Spitzer/nev_comp.dat', newx, spec_ohms_nev, spec_ohms_no_nev, spec_con_nev, spec_con_no_nev

; Results from measuring silicate in non-masing galaxies

;restore, '~/Astronomy/Research/Spitzer/nev_comp.dat';, newx, spec_ohms_nev, spec_ohms_no_nev, spec_con_nev, spec_con_no_nev
;wave = newx
;flux = spec_con_nev	; S_9.7 = 0.49		spl=2
;flux = spec_con_no_nev	; S_9.7 = 0.72		spl=7

end
