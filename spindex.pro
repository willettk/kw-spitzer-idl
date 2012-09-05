pro spindex, control = control, arch = arch, cso = cso, noplot = noplot, wait = wait
;+
; NAME: 
;       SPINDEX 
;
; PURPOSE:
;
;	Measure the spectral indices from 5-15 um and 20-30 um for Spitzer lo-res IRS spectra
;
; CATEGORY:
;	ASTRONOMY; DATABASE
;
; INPUTS:
;	
; OUTPUTS:
;
; KEYWORDS:
;
;	CONTROL - 	runs on the control sample objects (default is Darling sample of Spitzer OHMs)
;
;	ARCH - 		runs on the archived Spitzer OHMs
;
;	CSO - 		runs on the CSOs
;
;	NOPLOT - 	do not plot data with overlaid indices
;
;	WAIT - 		string giving pause between spectra (default is zero)
;
; REQUIRES:
;
;	CMW.pro
;
; EXAMPLE:
;
;	IDL> spindex,/arch
;
; NOTES:
;
; MODIFICATION HISTORY:
;
;	Written - KW, Sep/Oct 07
; 	Added keywords for control and archived samples - KW, Dec 07
;	Added CSO keyword		Feb 08
;	Properly computed error bars using propagation of errors - Jul 09
;-

if keyword_set(control) then fname = condat('tag') else $
	if keyword_set(arch) then fname = archdat('tag') else $
	if keyword_set(cso) then fname = csodat('tag') else $
	fname = ohmdat('tag')

nf = n_elements(fname)

alpha1 = fltarr(nf)
alpha2 = fltarr(nf)
alphaerr=fltarr(nf,2)

spindex_fname = fname

for i = 0, nf - 1 do begin

	tag,fname(i),dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname(i)+'.sav'
	

	a1_min = 5.3
	a1_max = 14.8
	a2_min = 20.0
	if max(sed.wave_lr) lt 30. then begin
		print,sed.tag+' redshifted past 30 um - using 28 um for a2_max'
		a2_max = 28.0 
	endif else a2_max = 30.0

	; Spectral index program
	
	; Locate flux measurements near desired wavelengths

	sp6  = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, a1_min)
	sp15 = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, a1_max)
	sp20 = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, a2_min)
	sp30 = cmw(sed.wave_lr, sed.flux_lr, sed.err_lr, a2_max)
	
	; Calculate spectral indices

	index1 = alog10(sp15/sp6) / alog10(a1_max / a1_min)
	index2 = alog10(sp30/sp20) / alog10(a2_max / a2_min)
		
	alpha1[i] = index1
	alpha2[i] = index2
	
	; Find errors in flux measurements

	sig_f6  = stddev(sed.flux_lr(closeto(sed.wave_lr,a1_min - 0.4):closeto(sed.wave_lr,a1_min + 0.4)))
	sig_f15 = stddev(sed.flux_lr(closeto(sed.wave_lr,a1_max - 0.8):closeto(sed.wave_lr,a1_max + 0.8)))
	sig_f20 = stddev(sed.flux_lr(closeto(sed.wave_lr,a2_min - 1.0):closeto(sed.wave_lr,a2_min + 1.0)))
	sig_f30 = stddev(sed.flux_lr(closeto(sed.wave_lr,a2_max - 1.0):closeto(sed.wave_lr,a2_max + 1.0)))

	; Calculate errors in slope

	sig_wave6  = abs(sed.wave_lr[closeto(sed.wave_lr,a1_min)] - sed.wave_lr[closeto(sed.wave_lr,a1_min) + 1]) 
	sig_wave15 = abs(sed.wave_lr[closeto(sed.wave_lr,a1_max)] - sed.wave_lr[closeto(sed.wave_lr,a1_max) + 1]) 
	sig_wave20 = abs(sed.wave_lr[closeto(sed.wave_lr,a2_min)] - sed.wave_lr[closeto(sed.wave_lr,a2_min) + 1]) 
	sig_wave30 = abs(sed.wave_lr[closeto(sed.wave_lr,a2_max)] - sed.wave_lr[closeto(sed.wave_lr,a2_max) + 1]) 

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

	alphaerr[i,*] = [sig_alpha1, sig_alpha2]

	if not keyword_set(noplot) then begin

		plot,sed.wave_lr, sed.flux_lr, $
			/xlog, /ylog, $
			xrange = [4,40], /xstyle, $
			xtitle = 'Wavelength [um]', $
			ytitle = 'Flux [Jy]', $
			title = sed.obj
			
		oplot, [a1_min], [sp6],  psym = symcat(14), symsize = 1.5, color=fsc_color("Red")
		oplot, [a1_max],[sp15], psym = symcat(14), symsize = 1.5, color=fsc_color("Red")
		oplot, [a2_min],[sp20], psym = symcat(14), symsize = 1.5, color=fsc_color("Yellow")
		oplot, [a2_max],[sp30], psym = symcat(14), symsize = 1.5, color=fsc_color("Yellow")
		
		plots,[a1_min,a1_max],[sp6,sp15], color=fsc_color("Red")
		plots,[a2_min,a2_max],[sp20,sp30], color=fsc_color("Yellow")
		
		xyouts,0.2,0.8,'!7a!3!I1!N = '+string(index1,format='(f5.2)'),color=fsc_color("Red"), /normal, charsize = 2
		xyouts,0.2,0.7,'!7a!3!I2!N = '+string(index2,format='(f5.2)'),color=fsc_color("Yellow"), /normal, charsize = 2

		xyouts,0.85,0.2,/normal,sed.tag,charsize=1.5
	
	endif
	
	if n_elements(wait) ne 0 then wait, wait
endfor

if keyword_set(control) then save,alpha1,alpha2,alphaerr,spindex_fname,file='~/Astronomy/Research/Spitzer/control/data/idl_sav/spindex_con.sav' else $
	if keyword_set(arch) then save,alpha1,alpha2,alphaerr,spindex_fname,file='~/Astronomy/Research/Spitzer/archived/data/idl_sav/spindex_arch.sav' else $
	if keyword_set(cso) then save,alpha1,alpha2,alphaerr,spindex_fname,file='~/Astronomy/Research/Spitzer/cso/data/idl_sav/spindex_cso.sav' else $
	save,alpha1,alpha2,alphaerr,spindex_fname,file='~/Astronomy/Research/Spitzer/OHM/spindex.sav' 

end
