pro speccal, fname, lr_frac, ps = ps, nowrite = nowrite, stop = stop, nospicehdr = nospicehdr, nobonus = nobonus, quiet = quiet
;+
; NAME:
;       SPECCAL
;
; PURPOSE:
;
;	Calibrate the stitched LR spectra on an absolute photometric scale, using available data
;
; INPUTS:
;
;	artificial_pu.csv - 	Excel spreadsheet containing measured artificial peakup
;				fluxes and true values for entire sample. 
;
; OUTPUTS:
;
;	LR_FRAC - 	multiplicative scaling factor used on LR spectra
;
; KEYWORDS:
;
;	PS - hard copies of overlaid raw and calibrated spectra
;
;	NOBONUS - do not plot bonus orders
;
;	NOWRITE - do not write text files of spectra to directory
;
;	NOSPICEHDR - do not run SPICEHDR to create files readable by SMART
;
;	STOP - stop program after running
;
; EXAMPLE:
;
;	IDL> speccal, 'mega001'
;
; REQUIRES:
;
;	TAG.pro
;	TARGETS.pro
;	READCOL,pro
;
; NOTES:
;
;	Program will multiplicatively scale the LR spectra as a single unit using available photometry. 
;
;	Priority of available photometry is:
;		- IRS dedicated peakup on target (exists only for Darling OHM sample)
;		- IRS target acquisition peakup on target (exists for portions of archive samples)
;		- IRAS 25 um flux (exists for portions of archive samples)
;		- no scaling
;
; REVISION HISTORY
;       Adapted from SPEC2PU_FLUX by K. Willett                Feb 2008
;	Make sure there are no non-data lines in peakup2.csv beyond the 3rd row - KW, Mar 08
;	Added switch for OHM/CSOs with no PU flux - KW, Mar 08
; 	Default is now to perform no scaling on the HR spectra (per sugg. from Lee and Henrik) - KW, May 08
;-

; Set device to read in colors

device, decomposed = 0

; Find directory and object name

tag, fname, dirtag 
targets, fname, redshift, obj

; Read in the spectrum

	specpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/stitched/coadd/'

	; Lores

	readcol, specpath+fname+'_sl1_cal.tbl', $
		det_sl1, wave_sl1, flux_sl1, err_sl1, bit_sl1, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_sl2_cal.tbl', $
		det_sl2, wave_sl2, flux_sl2, err_sl2, bit_sl2, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_sl3_cal.tbl', $
		det_sl3, wave_sl3, flux_sl3, err_sl3, bit_sl3, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll1_cal.tbl', $
		det_ll1, wave_ll1, flux_ll1, err_ll1, bit_ll1, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll2_cal.tbl', $
		det_ll2, wave_ll2, flux_ll2, err_ll2, bit_ll2, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_ll3_cal.tbl', $
		det_ll3, wave_ll3, flux_ll3, err_ll3, bit_ll3, format = 'i,f,f,f,i', skipline = 1, /silent

	; Hires

	readcol, specpath+fname+'_sh_cal.tbl', $
		det_sh, wave_sh, flux_sh, err_sh, bit_sh, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, specpath+fname+'_lh_cal.tbl', $
		det_lh, wave_lh, flux_lh, err_lh, bit_lh, format = 'i,f,f,f,i', skipline = 1, /silent

; Read in the measured and artificial peakup fluxes

apupath = '~/Astronomy/Research/Spitzer/'

readcol, apupath+'peakups3.csv', atag, aobj, apu_blue_lo, apu_red_lo, apu_iras12_lo, apu_iras25_lo, $
;	apu_blue_hi, apu_red_hi, my_bluepu, my_redpu, $
	format = 'a,a,a,a,a,a,', skipline = 3, /silent

artind = where(fname eq atag)
if artind eq -1 then begin
	message,'File not found in peakup spreadsheet'
	stop
endif

; Tests for best available photometric scale

; Dedicated (SUR) IRS peakups

datadir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+fname+'_peakup'
dirtest = file_search(datadir)

if dirtest(0) ne '' then istherepusur = 1 else istherepusur = 0

; Target acquisition (DCS) IRS PU

putargetfile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/pu_fluxes.sav'
putargettest = file_search(putargetfile)

isthereputarget = 0

if putargettest(0) ne '' and fname ne 'arch012' then begin	; arch012 did not peak up on a target

	restore,putargettest[0]

	tempind = where(fname eq pufilelist)
	if (tempind lt 0 or n_elements(tempind) ne 1) then begin
		message, 'Did not find fluxes in IDL sav file for target acquisition peakups'
		stop
	endif

	; Check if existing peakups are on science target

	targetind = where(pufilelist eq fname)
	if targetind(0) ne -1 then begin
		ispupos = pu_type(targetind(0))
		puflux_dcs = pu_flux(targetind(0)) * 1d-3
		this_filter = pu_filter(targetind(0))
	
		if ispupos eq 0 then isthereputarget = 1
;		if ispupos eq 0 then print, fname+' is an offset peakup (on science target)' else $
;			print, fname+' is a position peakup (using nearby star)'
	
	endif

endif else isthereputarget = 0

; IRAS 25 micron flux

irasfile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/'+dirtag+'_iras.sav'
irastest = file_search(irasfile)

if irastest(0) ne '' then begin

	restore, irastest(0)

	irasind = where(aname_ir eq strtrim(obj,2))
	if irasind(0) ne -1 then iras25 = a25(irasind(0)) else iras25 = 0

	if iras25 ne 0. then isthereiras = 1 else isthereiras = 0

endif else isthereiras = 0


; Apply the calibration based on available photometry

pu22_conv = apu_red_lo(artind) * 1d-3		; Jy
pu16_conv = apu_blue_lo(artind) * 1d-3		; Jy

if istherepusur then begin

	pufile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/pu_sur_'+dirtag+'.sav'
	restore, pufile

	surind = where(pufilelist eq fname)

	if surind(0) ne -1 then begin

		pu22_sur = pu22(surind)			; Jy

		lr_frac = pu22_sur(0) / pu22_conv(0)

	endif else begin
		message,fname+' not found in SUR peakup saved file'
		stop
	endelse

	if not keyword_set(quiet) then begin
		print,'Using SUR IRS 22 um peakup to calibrate '+fname
		print,'SUR / artificial = ',lr_frac
	endif
	
	wave_used = [22.]
	pu_used = [pu22_sur]

endif else if isthereputarget then begin

	if this_filter eq 'Red' then begin
		lr_frac = puflux_dcs / pu22_conv
		wave_used = [22.]
		pu_used = [puflux_dcs]
	endif else if this_filter eq 'Blu' then begin
		lr_frac = puflux_dcs / pu16_conv
		wave_used = [16.]
		pu_used = [puflux_dcs]
	endif else begin
		message,'Did not find correct filter for DCS peakup calibration'
		stop
	endelse

	lr_frac = lr_frac(0)

	if not keyword_set(quiet) then begin
		print,'Using target acquisition IRS peakups for '+fname
		print,'DCS / artificial = ',lr_frac,' for '+this_filter+' filter'
	endif

endif else if isthereiras then begin

	pu25_conv = apu_iras25_lo(artind) * 1d-3	; Jy

	lr_frac = iras25(0) / pu25_conv(0)

	if not keyword_set(quiet) then begin
		print,'Using IRAS photometry for '+fname
		print,'IRAS 25 um / artificial = ',lr_frac
	endif

	wave_used = [25.]
	pu_used = [iras25(0)]

endif else begin

	lr_frac = 1d
	if not keyword_set(quiet) then print,'No photometric calibration found for '+fname

	wave_used = [0]
	pu_used = [0]
endelse


; Open the individual nods

nodpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/stitched/nods/'

	; Lores

	readcol, nodpath+fname+'_sl1_1p_cal.tbl', $
		det_sl1_1p, wave_sl1_1p, flux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl1_2p_cal.tbl', $
		det_sl1_2p, wave_sl1_2p, flux_sl1_2p, err_sl1_2p, bit_sl1_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl2_1p_cal.tbl', $
		det_sl2_1p, wave_sl2_1p, flux_sl2_1p, err_sl2_1p, bit_sl2_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl2_2p_cal.tbl', $
		det_sl2_2p, wave_sl2_2p, flux_sl2_2p, err_sl2_2p, bit_sl2_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl3_1p_cal.tbl', $
		det_sl3_1p, wave_sl3_1p, flux_sl3_1p, err_sl3_1p, bit_sl3_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sl3_2p_cal.tbl', $
		det_sl3_2p, wave_sl3_2p, flux_sl3_2p, err_sl3_2p, bit_sl3_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll1_1p_cal.tbl', $
		det_ll1_1p, wave_ll1_1p, flux_ll1_1p, err_ll1_1p, bit_ll1_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll1_2p_cal.tbl', $
		det_ll1_2p, wave_ll1_2p, flux_ll1_2p, err_ll1_2p, bit_ll1_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll2_1p_cal.tbl', $
		det_ll2_1p, wave_ll2_1p, flux_ll2_1p, err_ll2_1p, bit_ll2_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll2_2p_cal.tbl', $
		det_ll2_2p, wave_ll2_2p, flux_ll2_2p, err_ll2_2p, bit_ll2_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll3_1p_cal.tbl', $
		det_ll3_1p, wave_ll3_1p, flux_ll3_1p, err_ll3_1p, bit_ll3_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_ll3_2p_cal.tbl', $
		det_ll3_2p, wave_ll3_2p, flux_ll3_2p, err_ll3_2p, bit_ll3_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	; Hires

	readcol, nodpath+fname+'_sh_1p_cal.tbl', $
		det_sh_1p, wave_sh_1p, flux_sh_1p, err_sh_1p, bit_sh_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_sh_2p_cal.tbl', $
		det_sh_2p, wave_sh_2p, flux_sh_2p, err_sh_2p, bit_sh_2p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_lh_1p_cal.tbl', $
		det_lh_1p, wave_lh_1p, flux_lh_1p, err_lh_1p, bit_lh_1p, format = 'i,f,f,f,i', skipline = 1, /silent

	readcol, nodpath+fname+'_lh_2p_cal.tbl', $
		det_lh_2p, wave_lh_2p, flux_lh_2p, err_lh_2p, bit_lh_2p, format = 'i,f,f,f,i', skipline = 1, /silent

; Multiply the spectra by the scaling fraction

calflux_sl1 = flux_sl1 * lr_frac
calflux_sl2 = flux_sl2 * lr_frac
calflux_sl3 = flux_sl3 * lr_frac
calflux_ll1 = flux_ll1 * lr_frac
calflux_ll2 = flux_ll2 * lr_frac
calflux_ll3 = flux_ll3 * lr_frac

calflux_sh = flux_sh
calflux_lh = flux_lh


calflux_sl1_1p = flux_sl1_1p * lr_frac
calflux_sl2_1p = flux_sl2_1p * lr_frac
calflux_sl3_1p = flux_sl3_1p * lr_frac
calflux_ll1_1p = flux_ll1_1p * lr_frac
calflux_ll2_1p = flux_ll2_1p * lr_frac
calflux_ll3_1p = flux_ll3_1p * lr_frac

calflux_sh_1p = flux_sh_1p
calflux_lh_1p = flux_lh_1p


calflux_sl1_2p = flux_sl1_2p * lr_frac
calflux_sl2_2p = flux_sl2_2p * lr_frac
calflux_sl3_2p = flux_sl3_2p * lr_frac
calflux_ll1_2p = flux_ll1_2p * lr_frac
calflux_ll2_2p = flux_ll2_2p * lr_frac
calflux_ll3_2p = flux_ll3_2p * lr_frac

calflux_sh_2p = flux_sh_2p
calflux_lh_2p = flux_lh_2p

; Plot results

!p.multi = [0,1,2]

; LORES spectra

pspath = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/calibration/'

yellow = fsc_color('Yellow')
green = fsc_color('Green')
red = fsc_color('Red')
blue = fsc_color('Blue')
orange = fsc_color('Orange')

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = pspath+fname+'_calflux.ps', /landscape, /color, bits_per_pixel = 8
	axiscolor = fsc_color('Black')
	bgcolor = fsc_color('White')
	cs = 1
	ct = 2
endif else begin
	axiscolor = fsc_color('White')
	bgcolor = fsc_color('Black')
	cs = 2
	ct = 1
endelse

plot, wave_sl1, flux_sl1, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
	xr = [4,40], /xstyle, $
	yr = [1d-4,1d1], /ystyle, $
	/xlog,/ylog, $
	title = obj, $
	charsize = cs, $
	color = axiscolor, $
	background = bgcolor, $
	linestyle = 2, $
	thick = ct, $
	charthick = ct, $
	/nodata


oplot, wave_sl1, flux_sl1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_sl2, flux_sl2, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll1, flux_ll1, color = axiscolor, linestyle = 2, thick = ct
oplot, wave_ll2, flux_ll2, color = axiscolor, linestyle = 2, thick = ct

oplot, wave_ll2, calflux_ll2, color = blue, thick = ct
oplot, wave_ll1, calflux_ll1, color = red, thick = ct
oplot, wave_sl2, calflux_sl2, color = yellow, thick = ct
oplot, wave_sl1, calflux_sl1, color = green, thick = ct

if not keyword_set(nobonus) then begin
	oplot, wave_sl3, flux_sl3, color = axiscolor, linestyle = 2, thick = ct
	oplot, wave_ll3, flux_ll3, color = axiscolor, linestyle = 2, thick = ct
	oplot, wave_ll3, calflux_ll3, color = orange, thick = ct
	oplot, wave_sl3, calflux_sl3, color = orange, thick = ct
endif

if wave_used ne 0 then oplot, wave_used, pu_used,   psym = symcat(14), color = red,  symsize = 2
oplot, [16], [pu16_conv],  psym = 4, color = axiscolor,  symsize = 2
oplot, [22], [pu22_conv],   psym = 4, color = axiscolor,  symsize = 2

xyouts, 0.2, 0.85, 'Lo-res', /normal, charsize = cs, charthick = ct, color = axiscolor

; HIRES spectra

plot, wave_sh, flux_sh, $
	xtitle = 'Wavelength [!7l!3m]', $
	ytitle = 'Flux [Jy]', $
	xr = [5,40], $
	yr = [0,max(flux_lh)], $
	title = obj, $
	charsize = cs, $
	linestyle = 2, $
	color = axiscolor, $
	background = bgcolor, $
	thick = ct, $
	charthick = ct, $
	/nodata

oplot, wave_sh, flux_sh, color = axiscolor, linestyle = 1, thick = ct, psym = 10
oplot, wave_lh, flux_lh, color = axiscolor, linestyle = 1, thick = ct, psym = 10

oplot, wave_sh, calflux_sh, color = blue, linestyle = 0, thick = ct, psym = 10
oplot, wave_lh, calflux_lh, color = red, linestyle = 0, thick = ct, psym = 10

xyouts, 0.2, 0.35, 'Hi-res', /normal, charsize = cs, charthick = ct, color = axiscolor
oplot, [16], [pu16_conv],  psym = 4, color = axiscolor,  symsize = 2
oplot, [22], [pu22_conv],   psym = 4, color = axiscolor,  symsize = 2
legend,/top,/right,['LR conv. PU'], psym=4

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

if not keyword_set(nowrite) then begin

	; Writing calibrated spectra to disk
	
	writepath_coadd = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
	writepath_nod   = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'

	; Coadded

	forprint, det_sl1, wave_sl1, calflux_sl1, err_sl1, bit_sl1, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl1_cal.tbl', /silent

	forprint, det_sl2, wave_sl2, calflux_sl2, err_sl2, bit_sl2, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl2_cal.tbl', /silent

	forprint, det_sl3, wave_sl3, calflux_sl3, err_sl3, bit_sl3, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sl3_cal.tbl', /silent

	forprint, det_ll1, wave_ll1, calflux_ll1, err_ll1, bit_ll1, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll1_cal.tbl', /silent

	forprint, det_ll2, wave_ll2, calflux_ll2, err_ll2, bit_ll2, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll2_cal.tbl', /silent

	forprint, det_ll3, wave_ll3, calflux_ll3, err_ll3, bit_ll3, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_ll3_cal.tbl', /silent

	forprint, det_sh, wave_sh, calflux_sh, err_sh, bit_sh, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_sh_cal.tbl', /silent

	forprint, det_lh, wave_lh, calflux_lh, err_lh, bit_lh, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_coadd+fname+'_lh_cal.tbl', /silent

	; Nods

	forprint, det_sl1_1p, wave_sl1_1p, calflux_sl1_1p, err_sl1_1p, bit_sl1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl1_1p_cal.tbl', /silent

	forprint, det_sl1_2p, wave_sl1_2p, calflux_sl1_2p, err_sl1_2p, bit_sl1_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl1_2p_cal.tbl', /silent

	forprint, det_sl2_1p, wave_sl2_1p, calflux_sl2_1p, err_sl2_1p, bit_sl2_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl2_1p_cal.tbl', /silent

	forprint, det_sl2_2p, wave_sl2_2p, calflux_sl2_2p, err_sl2_2p, bit_sl2_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl2_2p_cal.tbl', /silent

	forprint, det_sl3_1p, wave_sl3_1p, calflux_sl3_1p, err_sl3_1p, bit_sl3_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl3_1p_cal.tbl', /silent

	forprint, det_sl3_2p, wave_sl3_2p, calflux_sl3_2p, err_sl3_2p, bit_sl3_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sl3_2p_cal.tbl', /silent

	forprint, det_ll1_1p, wave_ll1_1p, calflux_ll1_1p, err_ll1_1p, bit_ll1_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll1_1p_cal.tbl', /silent

	forprint, det_ll1_2p, wave_ll1_2p, calflux_ll1_2p, err_ll1_2p, bit_ll1_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll1_2p_cal.tbl', /silent

	forprint, det_ll2_1p, wave_ll2_1p, calflux_ll2_1p, err_ll2_1p, bit_ll2_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll2_1p_cal.tbl', /silent

	forprint, det_ll2_2p, wave_ll2_2p, calflux_ll2_2p, err_ll2_2p, bit_ll2_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll2_2p_cal.tbl', /silent

	forprint, det_ll3_1p, wave_ll3_1p, calflux_ll3_1p, err_ll3_1p, bit_ll3_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll3_1p_cal.tbl', /silent

	forprint, det_ll3_2p, wave_ll3_2p, calflux_ll3_2p, err_ll3_2p, bit_ll3_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_ll3_2p_cal.tbl', /silent

	forprint, det_sh_1p, wave_sh_1p, calflux_sh_1p, err_sh_1p, bit_sh_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sh_1p_cal.tbl', /silent

	forprint, det_sh_2p, wave_sh_2p, calflux_sh_2p, err_sh_2p, bit_sh_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_sh_2p_cal.tbl', /silent

	forprint, det_lh_1p, wave_lh_1p, calflux_lh_1p, err_lh_1p, bit_lh_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_lh_1p_cal.tbl', /silent

	forprint, det_lh_2p, wave_lh_2p, calflux_lh_2p, err_lh_2p, bit_lh_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [Jy]         Uncert.        Bit type', $
		textout = writepath_nod+fname+'_lh_2p_cal.tbl', /silent

endif

if not keyword_set(nospicehdr) then spicehdr, fname, readdir = 'calibrated'

if not keyword_set(quiet) then begin
	print, 'Completed SPECCAL for ',fname
	print,''
endif

if keyword_set(stop) then stop

end
