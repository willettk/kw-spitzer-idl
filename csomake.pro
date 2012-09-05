pro csomake

;+
; NAME:
;       CSOMAKE
;
; PURPOSE:
;	Create IDL data structures with pertinent info for IRS data from CSOs
;
; OUTPUTS:
;
;	IDL structures for each CSO found in target directory
;
; EXAMPLE:
;
;	IDL> csomake
;
; NOTES:
;
;	Modeled off OHMMAKE.pro, the similar routine used for OHM data.
;
; REVISION HISTORY
;       Written by K. Willett                Oct 2007
;	Added refs. to IDL sav files instead of reading from spreadsheets - Apr 08
; 	Saved raw data from each nod - Jan 09
;	Silicate data changed to both 10 and 18 um measurements; errors now separate field - Jul 09
;-

; List of objects

fname = ['cso001', $
	'cso002', $
	'cso003', $
	'cso004', $
	'cso005', $
	'cso006',$
	'cso007',$
	'cso008',$
	'cso009',$
	'cso010']

nf = n_elements(fname)

; Peakup and IRAS photometry

	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/pu_sur_cso.sav'
	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/cso_iras.sav'

; Dust temperatures

	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/dustcso.sav'

; Spectral indices

	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/spindex_cso.sav'

; Silicate strengths and PAH EWs

	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/sil_pah_cso.sav'
	
; H2 temperatures and masses

	restore,'~/Astronomy/Research/Spitzer/cso/h2_cso.sav'

; D_L error

	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/dlerr.sav'

for i = 0, nf - 1 do begin

	; PAHFIT results

	restore,'~/Astronomy/Research/Spitzer/cso/pahfit/'+fname[i]+'_fit.sav'
	restore,'~/Astronomy/Research/Spitzer/cso/data/idl_sav/silicate/'+fname[i]+'_sil.sav'

	h2ind = where(h2.tag eq fname[i])

	; Read in the full LR spectrum
	
	tag,fname(i),dirtag
	targets,fname(i), z, obj, dl
	
	specdir_or = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
	specdir_bp = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/badpixelsremoved/'
	a = file_search(specdir_bp+'*'+fname[i]+'*')
	if a[0] ne '' then specdir = specdir_bp else specdir = specdir_or

;	specdir = specdir_or

	readcol, specdir+fname(i)+'_sl1_cal.tbl', $
			sl1_order, sl1_wave, sl1_flux, sl1_err, sl1_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname(i)+'_sl2_cal.tbl', $
		sl2_order, sl2_wave, sl2_flux, sl2_err, sl2_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname(i)+'_sl3_cal.tbl', $
		sl3_order, sl3_wave, sl3_flux, sl3_err, sl3_bit, format = 'i,f,f,f,i', /silent

	readcol, specdir+fname(i)+'_ll1_cal.tbl', $
		ll1_order, ll1_wave, ll1_flux, ll1_err, ll1_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname(i)+'_ll2_cal.tbl', $
		ll2_order, ll2_wave, ll2_flux, ll2_err, ll2_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname(i)+'_ll3_cal.tbl', $
		ll3_order, ll3_wave, ll3_flux, ll3_err, ll3_bit, format = 'i,f,f,f,i', /silent


	readcol, specdir+fname(i)+'_sh_cal.tbl', $
		sh_order, sh_wave, sh_flux, sh_err, sh_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname(i)+'_lh_cal.tbl', $
		lh_order, lh_wave, lh_flux, lh_err, lh_bit, format = 'i,f,f,f,i', /silent

	wave1   = [sl2_wave, sl1_wave, ll2_wave, ll1_wave]
	flux1   = [sl2_flux, sl1_flux, ll2_flux, ll1_flux]
	err1    = [sl2_err, sl1_err, ll2_err, ll1_err]
	order1  = [sl2_order, sl1_order, ll2_order, ll1_order]
	bit1    = [sl2_bit, sl1_bit, ll2_bit, ll1_bit]
	module1 = [replicate('SL2',n_elements(sl2_order)), replicate('SL1',n_elements(sl1_order)), $
		replicate('LL2',n_elements(ll2_order)), replicate('LL1',n_elements(ll1_order))]

	waves1 = sort(wave1)
	
	wave_lr = wave1(waves1)
	flux_lr = flux1(waves1)
	err_lr  = err1(waves1)
	order_lr  = order1(waves1)
	bit_lr  = bit1(waves1)
	module_lr  = module1(waves1)

	nlr = n_elements(flux_lr)
	nsh = n_elements(sh_flux)
	nlh = n_elements(lh_flux)
	
	; Nods

;		specdir_nod = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'
;	
;		readcol, specdir_nod+fname(i)+'_sl1_1p_cal.tbl', $
;			sl1_order_1p, sl1_wave_1p, sl1_flux_1p, sl1_err_1p, sl1_bit_1p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_sl2_1p_cal.tbl', $
;			sl2_order_1p, sl2_wave_1p, sl2_flux_1p, sl2_err_1p, sl2_bit_1p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_sl3_1p_cal.tbl', $
;			sl3_order_1p, sl3_wave_1p, sl3_flux_1p, sl3_err_1p, sl3_bit_1p, format = 'i,f,f,f,i', /silent
;	
;		readcol, specdir_nod+fname(i)+'_ll1_1p_cal.tbl', $
;			ll1_order_1p, ll1_wave_1p, ll1_flux_1p, ll1_err_1p, ll1_bit_1p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_ll2_1p_cal.tbl', $
;			ll2_order_1p, ll2_wave_1p, ll2_flux_1p, ll2_err_1p, ll2_bit_1p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_ll3_1p_cal.tbl', $
;			ll3_order_1p, ll3_wave_1p, ll3_flux_1p, ll3_err_1p, ll3_bit_1p, format = 'i,f,f,f,i', /silent
;	
;	
;		readcol, specdir_nod+fname(i)+'_sh_1p_cal.tbl', $
;			sh_order_1p, sh_wave_1p, sh_flux_1p, sh_err_1p, sh_bit_1p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_lh_1p_cal.tbl', $
;			lh_order_1p, lh_wave_1p, lh_flux_1p, lh_err_1p, lh_bit_1p, format = 'i,f,f,f,i', /silent
;	
;		wave1   = [sl2_wave_1p, sl1_wave_1p, ll2_wave_1p, ll1_wave_1p]
;		flux1   = [sl2_flux_1p, sl1_flux_1p, ll2_flux_1p, ll1_flux_1p]
;		err1    = [sl2_err_1p, sl1_err_1p, ll2_err_1p, ll1_err_1p]
;		order1  = [sl2_order, sl1_order, ll2_order, ll1_order]
;		bit1    = [sl2_bit_1p, sl1_bit_1p, ll2_bit_1p, ll1_bit_1p]
;		module1 = [replicate('SL2',n_elements(sl2_order_1p)), replicate('SL1',n_elements(sl1_order_1p)), $
;			replicate('LL2',n_elements(ll2_order_1p)), replicate('LL1',n_elements(ll1_order_1p))]
;	
;		waves1 = sort(wave1)
;		
;		wave_lr_1p = wave1(waves1)
;		flux_lr_1p = flux1(waves1)
;		err_lr_1p  = err1(waves1)
;		order_lr_1p  = order1(waves1)
;		bit_lr_1p  = bit1(waves1)
;		module_lr_1p  = module1(waves1)
;
;		readcol, specdir_nod+fname(i)+'_sl1_2p_cal.tbl', $
;			sl1_order_2p, sl1_wave_2p, sl1_flux_2p, sl1_err_2p, sl1_bit_2p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_sl2_2p_cal.tbl', $
;			sl2_order_2p, sl2_wave_2p, sl2_flux_2p, sl2_err_2p, sl2_bit_2p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_sl3_2p_cal.tbl', $
;			sl3_order_2p, sl3_wave_2p, sl3_flux_2p, sl3_err_2p, sl3_bit_2p, format = 'i,f,f,f,i', /silent
;	
;		readcol, specdir_nod+fname(i)+'_ll1_2p_cal.tbl', $
;			ll1_order_2p, ll1_wave_2p, ll1_flux_2p, ll1_err_2p, ll1_bit_2p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_ll2_2p_cal.tbl', $
;			ll2_order_2p, ll2_wave_2p, ll2_flux_2p, ll2_err_2p, ll2_bit_2p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_ll3_2p_cal.tbl', $
;			ll3_order_2p, ll3_wave_2p, ll3_flux_2p, ll3_err_2p, ll3_bit_2p, format = 'i,f,f,f,i', /silent
;	
;	
;		readcol, specdir_nod+fname(i)+'_sh_2p_cal.tbl', $
;			sh_order_2p, sh_wave_2p, sh_flux_2p, sh_err_2p, sh_bit_2p, format = 'i,f,f,f,i', /silent
;		
;		readcol, specdir_nod+fname(i)+'_lh_2p_cal.tbl', $
;			lh_order_2p, lh_wave_2p, lh_flux_2p, lh_err_2p, lh_bit_2p, format = 'i,f,f,f,i', /silent
;	
;		wave1   = [sl2_wave_2p, sl1_wave_2p, ll2_wave_2p, ll1_wave_2p]
;		flux1   = [sl2_flux_2p, sl1_flux_2p, ll2_flux_2p, ll1_flux_2p]
;		err1    = [sl2_err_2p, sl1_err_2p, ll2_err_2p, ll1_err_2p]
;		order1  = [sl2_order, sl1_order, ll2_order, ll1_order]
;		bit1    = [sl2_bit_2p, sl1_bit_2p, ll2_bit_2p, ll1_bit_2p]
;		module1 = [replicate('SL2',n_elements(sl2_order_2p)), replicate('SL1',n_elements(sl1_order_2p)), $
;			replicate('LL2',n_elements(ll2_order_2p)), replicate('LL1',n_elements(ll1_order_2p))]
;	
;		waves1 = sort(wave1)
;		
;		wave_lr_2p = wave1(waves1)
;		flux_lr_2p = flux1(waves1)
;		err_lr_2p  = err1(waves1)
;		order_lr_2p  = order1(waves1)
;		bit_lr_2p  = bit1(waves1)
;		module_lr_2p  = module1(waves1)

	sed = {obj:'', tag:'', date:'', inst:'', frame:'', redshift:0d, dl:0d, dlerr:0d, $
		dtemp:fltarr(2), angsize:fltarr(2), physize:fltarr(2), $
		h2temp:fltarr(2), h2mass:fltarr(2), h2temp_err: fltarr(2), $
		spindex:fltarr(2), spindexerr:fltarr(2), $
		sil:fltarr(2), $
		silerr:fltarr(2), $
		ice:fltarr(2), $
		pah62ew:fltarr(2), pah62ew_ice:fltarr(2), pah11ew:fltarr(2), $
		pah62flux:0d, pah11flux:0d, pah62lum:0d, pah11lum:0d, $
		pah62flux_lim:0d, pah62ew_lim:0d, pah11flux_lim:0d, pah11ew_lim:0d, pah62lum_lim:0d, pah11lum_lim:0d, $
		pahfit62ew:fltarr(2), pahfit62ew_ice:fltarr(2), pahfit11ew:fltarr(2), $
		pahfit62flux:0d, pahfit11flux:0d, pahfit62lum:0d, pahfit11lum:0d, $
		peakup:fltarr(2), peakuperr:fltarr(2), iras:fltarr(4), iraserr:fltarr(4), $
		lir:0d, lfir:0d, $
		wave_lr:fltarr(nlr), flux_lr:fltarr(nlr), err_lr:fltarr(nlr), $
		order_lr:intarr(nlr), bit_lr:intarr(nlr), module_lr:strarr(nlr), $
		wave_sh:fltarr(nsh), flux_sh:fltarr(nsh), err_sh:fltarr(nsh), $
		order_sh:intarr(nsh), bit_sh:intarr(nsh), module_sh:strarr(nsh), $
		wave_lh:fltarr(nlh), flux_lh:fltarr(nlh), err_lh:fltarr(nlh), $
		order_lh:intarr(nlh), bit_lh:intarr(nlh), module_lh:strarr(nlh)}
;		wave_lr_1p:fltarr(nlr), flux_lr_1p:fltarr(nlr), err_lr_1p:fltarr(nlr), $
;		order_lr_1p:intarr(nlr), bit_lr_1p:intarr(nlr), module_lr_1p:strarr(nlr), $
;		wave_sh_1p:fltarr(nsh), flux_sh_1p:fltarr(nsh), err_sh_1p:fltarr(nsh), $
;		order_sh_1p:intarr(nsh), bit_sh_1p:intarr(nsh), module_sh_1p:strarr(nsh), $
;		wave_lh_1p:fltarr(nlh), flux_lh_1p:fltarr(nlh), err_lh_1p:fltarr(nlh), $
;		order_lh_1p:intarr(nlh), bit_lh_1p:intarr(nlh), module_lh_1p:strarr(nlh), $
;		wave_lr_2p:fltarr(nlr), flux_lr_2p:fltarr(nlr), err_lr_2p:fltarr(nlr), $
;		order_lr_2p:intarr(nlr), bit_lr_2p:intarr(nlr), module_lr_2p:strarr(nlr), $
;		wave_sh_2p:fltarr(nsh), flux_sh_2p:fltarr(nsh), err_sh_2p:fltarr(nsh), $
;		order_sh_2p:intarr(nsh), bit_sh_2p:intarr(nsh), module_sh_2p:strarr(nsh), $
;		wave_lh_2p:fltarr(nlh), flux_lh_2p:fltarr(nlh), err_lh_2p:fltarr(nlh), $
;		order_lh_2p:intarr(nlh), bit_lh_2p:intarr(nlh), module_lh_2p:strarr(nlh)}
	
	irasind = where(fname[i] eq atag_ir, icount) & if icount ne 1 then $
		message,'Did not find correct IRAS fluxes in CSOMAKE for '+fname[i],/info
	puind = where(fname[i] eq pufilelist, pcount)& if pcount ne 1 then $
		message,'Did not find correct peakup fluxes in CSOMAKE for '+fname[i],/info
	silpahind = where(fname[i] eq csofiles_silpah, scount)& if scount ne 1 then $
		message,'Did not find correct silicate/PAH data in CSOMAKE for '+fname[i],/info


	; Shift to rest-frame wavelengths
	
	sed.obj = obj
	sed.tag = fname(i)
	sed.date = systime()
	sed.inst = 'IRS'
	sed.frame = 'Rest frame'
	sed.redshift = z
	sed.dl = dl
	sed.dlerr = dlerrarr[1,i]

	sed.dtemp = [dtemp(i,0),dtemp(i,1)]
	sed.angsize = [angsize(i,0),angsize(i,1)]
	sed.physize = [physize(i,0),physize(i,1)]

  	sed.h2temp = [h2.temp_warm[h2ind],h2.temp_hot[h2ind]]
  	sed.h2mass = [h2.mass_warm[h2ind],h2.mass_hot[h2ind]]
  	sed.h2temp_err = [h2.temp_warmerr[h2ind],h2.temp_hoterr[h2ind]]

	sed.spindex  = [alpha1(i), alpha2(i)]
	sed.spindexerr  = alphaerr(i,*)

	sed.sil    = [sil10[silpahind],sil18[silpahind]]
	sed.silerr = [sil10err[silpahind],sil18err[silpahind]]
	sed.ice = [ice[silpahind],iceerr[silpahind]]

	sed.pah62ew_ice = [pah62ew_ice[silpahind],pah62ew_iceerr[silpahind]]
	sed.pah62ew = [pah62ew[silpahind],pah62ew_err[silpahind]]
	sed.pah11ew = [pah11ew[silpahind],pah11ew_err[silpahind]]
	sed.pah62flux = pah62flux[silpahind]
	sed.pah11flux = pah11flux[silpahind]
	sed.pah62lum = alog10(pah62flux[silpahind] * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)
	sed.pah11lum = alog10(pah11flux[silpahind] * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)
	sed.pah62flux_lim = pah62flux_lim[silpahind]
	sed.pah11flux_lim = pah11flux_lim[silpahind]
	sed.pah62lum_lim = alog10(sed.pah62flux_lim * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)		
	sed.pah11lum_lim = alog10(sed.pah11flux_lim * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)
	sed.pah62ew_lim = pah62ew_lim[silpahind]
	sed.pah11ew_lim = pah11ew_lim[silpahind]

	sed.pahfit62ew_ice = [fit.dust_features[2].int_strength / sil_spline[closetomed(sil_wave,sil_spline,6.22)] * (6.22)^2 / 3d14, 0]
	sed.pahfit62ew = [fit.dust_features[2].eqw , 0]					; um
	sed.pahfit11ew = [(fit.dust_features[10].eqw + fit.dust_features[11].eqw) , 0]
	sed.pahfit62flux = fit.dust_features[2].int_strength * 1d-30
	sed.pahfit11flux = (fit.dust_features[10].int_strength + fit.dust_features[11].int_strength) * 1d-30	; W/cm^2
	sed.pahfit62lum = alog10(sed.pahfit62flux * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)		; log L_sun
	sed.pahfit11lum = alog10(sed.pahfit11flux * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)

	sed.peakup = [pu16[puind],pu22[puind]]
	sed.peakuperr = [pu16_err[puind],pu22_err[puind]]
	sed.iras = [a12[irasind],a25[irasind],a60[irasind],a100[irasind]]
	sed.iraserr = [a12_err[irasind],a25_err[irasind],a60_err[irasind],a100_err[irasind]]
	sed.lir = lir(a12[irasind], a25[irasind],a60[irasind],a100[irasind],dl)
	sed.lfir = lfir(a60[irasind],a100[irasind],dl)

	sed.wave_lr = wave_lr / (1d + z)
	sed.flux_lr = flux_lr
	sed.err_lr = err_lr
	sed.bit_lr = bit_lr
	sed.order_lr = order_lr
	sed.module_lr = module_lr
	sed.wave_sh = sh_wave / (1d + z)
	sed.flux_sh = sh_flux
	sed.err_sh = sh_err
	sed.bit_sh = sh_bit
	sed.order_sh = sh_order
	sed.wave_lh = lh_wave / (1d + z)
	sed.flux_lh = lh_flux
	sed.err_lh = lh_err
	sed.bit_lh = lh_bit
	sed.order_lh = lh_order

;		sed.wave_lr_1p = wave_lr_1p / (1d + z)
;		sed.flux_lr_1p = flux_lr_1p
;		sed.err_lr_1p = err_lr_1p
;		sed.bit_lr_1p = bit_lr_1p
;		sed.order_lr_1p = order_lr_1p
;		sed.module_lr_1p = module_lr_1p
;		sed.wave_sh_1p = sh_wave_1p / (1d + z)
;		sed.flux_sh_1p = sh_flux_1p
;		sed.err_sh_1p = sh_err_1p
;		sed.bit_sh_1p = sh_bit_1p
;		sed.order_sh_1p = sh_order_1p
;		sed.wave_lh_1p = lh_wave_1p / (1d + z)
;		sed.flux_lh_1p = lh_flux_1p
;		sed.err_lh_1p = lh_err_1p
;		sed.bit_lh_1p = lh_bit_1p
;		sed.order_lh_1p = lh_order_1p

;		sed.wave_lr_2p = wave_lr_2p / (1d + z)
;		sed.flux_lr_2p = flux_lr_2p
;		sed.err_lr_2p = err_lr_2p
;		sed.bit_lr_2p = bit_lr_2p
;		sed.order_lr_2p = order_lr_2p
;		sed.module_lr_2p = module_lr_2p
;		sed.wave_sh_2p = sh_wave_2p / (1d + z)
;		sed.flux_sh_2p = sh_flux_2p
;		sed.err_sh_2p = sh_err_2p
;		sed.bit_sh_2p = sh_bit_2p
;		sed.order_sh_2p = sh_order_2p
;		sed.wave_lh_2p = lh_wave_2p / (1d + z)
;		sed.flux_lh_2p = lh_flux_2p
;		sed.err_lh_2p = lh_err_2p
;		sed.bit_lh_2p = lh_bit_2p
;		sed.order_lh_2p = lh_order_2p

	save,sed,filename='~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname[i]+'.sav'

endfor

; Save a file with the sizes of each array so I don't have to enter it manually in OHMDAT

sedsize = {obj:1, tag:1, date:1, inst:1, frame:1, redshift:1, dl:1, dlerr:1, $
		dtemp:2, angsize:2, physize:2, $
		h2temp:2, h2mass:2, $
		h2temp_err:2, $
		spindex:2, spindexerr:2, $
		sil:2, silerr:2, $
		ice:2, $
		pah62ew:2, pah62ew_ice:2, pah11ew:2, pah62flux:1, pah11flux:1, pah62lum:1, pah11lum:1, $
		pah62flux_lim:1, pah62ew_lim:1, pah11flux_lim:1, pah11ew_lim:1, pah62lum_lim:1, pah11lum_lim:1, $
		pahfit62ew:2, pahfit62ew_ice:2, pahfit11ew:2, pahfit62flux:1, pahfit11flux:1, pahfit62lum:1, pahfit11lum:1, $
		peakup:2, peakuperr:2, iras:4, iraserr:4, $
		lir:1, lfir:1, $
		wave_lr:nlr, flux_lr:nlr, err_lr:nlr, order_lr:nlr, bit_lr:nlr, module_lr:nlr, $
		wave_sh:nsh, flux_sh:nsh, err_sh:nsh, order_sh:nsh, bit_sh:nsh, module_sh:nsh, $
		wave_lh:nlh, flux_lh:nlh, err_lh:nlh, order_lh:nlh, bit_lh:nlh, module_lh:nlh}
;		wave_lr_1p:nlr, flux_lr_1p:nlr, err_lr_1p:nlr, order_lr_1p:nlr, bit_lr_1p:nlr, module_lr_1p:nlr, $
;		wave_sh_1p:nsh, flux_sh_1p:nsh, err_sh_1p:nsh, order_sh_1p:nsh, bit_sh_1p:nsh, module_sh_1p:nsh, $
;		wave_lh_1p:nlh, flux_lh_1p:nlh, err_lh_1p:nlh, order_lh_1p:nlh, bit_lh_1p:nlh, module_lh_1p:nlh, $
;		wave_lr_2p:nlr, flux_lr_2p:nlr, err_lr_2p:nlr, order_lr_2p:nlr, bit_lr_2p:nlr, module_lr_2p:nlr, $
;		wave_sh_2p:nsh, flux_sh_2p:nsh, err_sh_2p:nsh, order_sh_2p:nsh, bit_sh_2p:nsh, module_sh_2p:nsh, $
;		wave_lh_2p:nlh, flux_lh_2p:nlh, err_lh_2p:nlh, order_lh_2p:nlh, bit_lh_2p:nlh, module_lh_2p:nlh}
	
save,sedsize,filename='~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/sedsize.sav'

end
