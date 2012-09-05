pro ohmmake
;+
; NAME: 
;       OHMMAKE 
;
; PURPOSE:
;
;	Create an IDL structure for each OHM containing full Spitzer spectra plus ancillary data
;
; CATEGORY:
;	ASTRONOMY; DATABASE; BATCH FILE
;
; INPUTS:
;	
; OUTPUTS:
;
; KEYWORDS:
;
; REQUIRES:
;
;	The data must already exist for this program to write it to a structure; this can be in a variety of forms, 
;		including spreadsheets, IDL .sav files, and ASCII files. 
;
; EXAMPLE:
;
;	IDL> ohmmake
;
; NOTES:
;
;	The first time this is run with new objects added, the sections for ancillary data (dust temperatures, H2 temps/masses, 
;		photometry, L_OH, spectral indices, silicates and PAH fluxes) must be commented out and the program run to create
;		a structure with just the Spitzer spectra. The individual programs are then run one at a time (ex. DUSTTEMP.pro, 
;		SPINDEX.pro) to generate the data and save it as IDL .sav files. The relevant lines are then uncommented here
;		and run again to make a structure that can be read by OHMDAT.pro. This isn't the best way of doing things, but it
;		will work until the amount of data becomes overly cumbersome. 
;
; MODIFICATION HISTORY:
;
;	Renamed from STRMAKE - Aug 08
;	PAH data now read from both spline fit and PAHFIT files - Sep 08
;	Spectra are now read from BADPIXELSREMOVED directory, if present - Oct 08
;	Removed mega025, mega033 - Jun 10
;-

fname = ['mega001', 'mega002', 'mega004', 'mega005', 'mega006', $
	'mega007', 'mega008', 'mega009', 'mega010', 'mega013', $
	'mega014', 'mega016', 'mega017', 'mega018', 'mega020', $
	'mega022', 'mega023',            'mega026', 'mega027', $
	'mega028', 'mega029', 'mega031', 'mega032',            $
	'mega034']

;fname = ['mega001', 'mega002', 'mega004', 'mega005', 'mega006', $
;	'mega007', 'mega008', 'mega009', 'mega010', 'mega013', $
;	'mega014', 'mega016', 'mega017', 'mega018', 'mega020', $
;	'mega022', 'mega023', 'mega025', 'mega026', 'mega027', $
;	'mega028', 'mega029', 'mega031', 'mega032', 'mega033', $
;	'mega034']

nf = n_elements(fname)

; Peakup and IRAS photometry

;	readcol,'~/Astronomy/Research/Spitzer/OHM/fluxes.csv',$
;		obj, tag, pu16, pu16err, pu22, pu22err, $
;		skipline = 4, format='a,a,f,f,f,f',/silent, numline=26
	
;	nobj = n_elements(obj)
;	pu16err = fltarr(nobj)+0.2	; Assume 20% error on all peakups (conservative)
;	pu22err = fltarr(nobj)+0.2
	
;	pu16 = pu16 * 1d-3 & pu16err = pu16err * pu16
;	pu22 = pu22 * 1d-3 & pu22err = pu22err * pu22

	restore,'~/Astronomy/Research/Spitzer/ohm/data/idl_sav/pu_sur_ohm.sav'
	restore,'~/Astronomy/Research/Spitzer/ohm/data/idl_sav/ohm_iras.sav'

; OH luminosity (Darling et al)

	restore,'~/Astronomy/Research/Spitzer/OHM/ohlum.sav'

; Dust temperatures

	restore,'~/Astronomy/Research/Spitzer/OHM/dusttemp.sav'

; Spectral indices

	restore,'~/Astronomy/Research/Spitzer/OHM/spindex.sav'

; Silicate strengths and PAH EWs

	restore,'~/Astronomy/Research/Spitzer/ohm/data/idl_sav/sil_pah_ohm.sav'

; H2 temperatures and masses

	restore,'~/Astronomy/Research/Spitzer/OHM/h2.sav'

zz=0
if zz eq 1 then begin

; Loop over each target in the list

endif
for i = 0, nf - 1 do begin

	; PAHFIT results

	restore,'~/Astronomy/Research/Spitzer/ohm/pahfit/'+fname[i]+'_fit.sav'
	restore,'~/Astronomy/Research/Spitzer/ohm/data/idl_sav/silicate/'+fname[i]+'_sil.sav'

	; Read in the full LR spectrum
	
	tag,fname[i],dirtag
	targets,fname[i], z, obj, dl
	
	specdir_or = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
	specdir_bp = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/badpixelsremoved/'
	a = file_search(specdir_bp+'*'+fname[i]+'*')
	if a[0] ne '' then specdir = specdir_bp else specdir = specdir_or

	readcol, specdir+fname[i]+'_sl1_cal.tbl', $
		sl1_order, sl1_wave, sl1_flux, sl1_err, sl1_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname[i]+'_sl2_cal.tbl', $
		sl2_order, sl2_wave, sl2_flux, sl2_err, sl2_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname[i]+'_sl3_cal.tbl', $
		sl3_order, sl3_wave, sl3_flux, sl3_err, sl3_bit, format = 'i,f,f,f,i', /silent

	readcol, specdir+fname[i]+'_ll1_cal.tbl', $
		ll1_order, ll1_wave, ll1_flux, ll1_err, ll1_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname[i]+'_ll2_cal.tbl', $
		ll2_order, ll2_wave, ll2_flux, ll2_err, ll2_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname[i]+'_ll3_cal.tbl', $
		ll3_order, ll3_wave, ll3_flux, ll3_err, ll3_bit, format = 'i,f,f,f,i', /silent

	readcol, specdir+fname[i]+'_sh_cal.tbl', $
		sh_order, sh_wave, sh_flux, sh_err, sh_bit, format = 'i,f,f,f,i', /silent
	
	readcol, specdir+fname[i]+'_lh_cal.tbl', $
		lh_order, lh_wave, lh_flux, lh_err, lh_bit, format = 'i,f,f,f,i', /silent

	wave1 = [sl2_wave, sl1_wave, ll2_wave, ll1_wave]
	flux1 = [sl2_flux, sl1_flux, ll2_flux, ll1_flux]
	err1  = [sl2_err, sl1_err, ll2_err, ll1_err]
	order1  = [sl2_order, sl1_order, ll2_order, ll1_order]
	bit1  = [sl2_bit, sl1_bit, ll2_bit, ll1_bit]
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
	
	sed = {obj:'', tag:'', date:'', inst:'', frame:'', redshift:0d, redshift_err:0d, dl:0d, $
		logoh:0d, f1667:0d, f1420:0d, czoh:0d, czoh_err:0d, cz_opt:0L, cz_opt_err:0L, rh:0d, $
		dtemp:fltarr(2), angsize:fltarr(2), physize:fltarr(2), $
		h2temp:fltarr(2), h2mass:fltarr(2), $
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
	

	irasind = where(fname[i] eq atag_ir, icount) & if icount ne 1 then $
		message,'Did not find correct IRAS fluxes in OHMMAKE for '+fname[i],/info
	puind = where(fname[i] eq pufilelist, pcount)& if pcount ne 1 then $
		message,'Did not find correct peakup fluxes in OHMMAKE for '+fname[i],/info
	ohind = where(fname[i] eq meganames, ocount)& if ocount ne 1 then $
		message,'Did not find correct OH data in OHMMAKE for '+fname[i],/info
	silpahind = where(fname[i] eq ohmfiles_silpah, scount)& if scount ne 1 then $
		message,'Did not find correct silicate/PAH data in OHMMAKE for '+fname[i],/info
	dtempind = where(fname[i] eq dtemp_fnames, dtcount)& if dtcount ne 1 then $
		message,'Did not find correct dust temperature data in OHMMAKE for '+fname[i],/info
	spindexind = where(fname[i] eq spindex_fname, spcount)& if spcount ne 1 then $
		message,'Did not find correct spectral index data in OHMMAKE for '+fname[i],/info
	h2ind = where(fname[i] eq h2.tag, h2count)

	; Shift to rest-frame wavelengths
	
	sed.obj = obj
	sed.tag = fname[i]
	sed.date = systime()
	sed.inst = 'IRS'
	sed.frame = 'Rest frame'
	sed.redshift = z
	sed.redshift_err = zerr_opt_ohm[ohind]
	sed.dl = dl
	sed.logoh = logoh[ohind]
	sed.f1667 = f1667_ohm[ohind]
	sed.f1420 = f1420_ohm[ohind]
	sed.rh = mega_rh[ohind]

	sed.czoh = czoh_ohm[ohind]
	sed.czoh_err = czoh_err_ohm[ohind]
	sed.cz_opt = cz_opt_ohm[ohind]
	sed.cz_opt_err = cz_opt_err_ohm[ohind]

	sed.dtemp = [dtemp(dtempind,0),dtemp(dtempind,1)]
	sed.angsize = [angsize(dtempind,0),angsize(dtempind,1)]
	sed.physize = [physize(dtempind,0),physize(dtempind,1)]

	sed.spindex  = [alpha1[spindexind], alpha2[spindexind]]		; alpha1 = 15-6, alpha2 = 30-20
	sed.spindexerr  = alphaerr(spindexind,*)

	sed.sil    = [sil10[silpahind],sil18[silpahind]]
	sed.silerr = [sil10err[silpahind],sil18err[silpahind]]
	sed.ice = [ice[silpahind],iceerr[silpahind]]

	if h2count ne 1 then begin
  		sed.h2temp = [0,0]
  		sed.h2mass = [0,0]
	endif else begin
	  	sed.h2temp = [h2.temp_warm[h2ind],h2.temp_hot[h2ind]]
	  	sed.h2mass = [h2.mass_warm[h2ind],h2.mass_hot[h2ind]]
	endelse

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
	sed.pahfit62lum = alog10(sed.pahfit62flux * (dl * 3.086d24)^2 * 4d * !dpi * 1d7 / 3.862d33)			; log L_sun
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

	save,sed,filename='~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname[i]+'.sav'

endfor

; Save a file with the sizes of each array so I don't have to enter it manually in OHMDAT

sedsize = {obj:1, tag:1, date:1, inst:1, frame:1, redshift:1, redshift_err:1, dl:1, $
		logoh:1, f1667:1, f1420:1, czoh:1, czoh_err:1, cz_opt:1, cz_opt_err:1, rh:1, $
		dtemp:2, angsize:2, physize:2, $
		h2temp:2, h2mass:2, $
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
	
save,sedsize,filename='~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/sedsize.sav'

end
