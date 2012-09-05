
;+
; NAME:
;       
;	PKS1413
;
; PURPOSE:
;
;	Plot zoomed-in nods of possible features in PKS 1413+135 IRS spectrum	
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
;	IDL> .r pks1413
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                
;-

set_plot,'ps'
;device, filename = '~/Astronomy/Research/Spitzer/plots/pks1413_zoom.ps', /landscape
device, filename = '~/Astronomy/Research/Spitzer/plots/pks1413_diffspec.ps', /landscape

!p.multi=[0,3,3]



; Difference spectra


		tag, 'cso005', dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+'cso005'+'.sav'

		specdir_nod = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'
		specdir_coadd = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/badpixelsremoved/'
 	
 		readcol, specdir_nod+'cso005'+'_sl1_1p_cal.tbl', $
 			sl1_order_1p, sl1_wave_1p, sl1_flux_1p, sl1_err_1p, sl1_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_sl2_1p_cal.tbl', $
 			sl2_order_1p, sl2_wave_1p, sl2_flux_1p, sl2_err_1p, sl2_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_sl3_1p_cal.tbl', $
 			sl3_order_1p, sl3_wave_1p, sl3_flux_1p, sl3_err_1p, sl3_bit_1p, format = 'i,f,f,f,i', /silent
 	
 		readcol, specdir_nod+'cso005'+'_ll1_1p_cal.tbl', $
 			ll1_order_1p, ll1_wave_1p, ll1_flux_1p, ll1_err_1p, ll1_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_ll2_1p_cal.tbl', $
 			ll2_order_1p, ll2_wave_1p, ll2_flux_1p, ll2_err_1p, ll2_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_ll3_1p_cal.tbl', $
 			ll3_order_1p, ll3_wave_1p, ll3_flux_1p, ll3_err_1p, ll3_bit_1p, format = 'i,f,f,f,i', /silent
 	
 	
 		readcol, specdir_nod+'cso005'+'_sh_1p_cal.tbl', $
 			sh_order_1p, sh_wave_1p, sh_flux_1p, sh_err_1p, sh_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_lh_1p_cal.tbl', $
 			lh_order_1p, lh_wave_1p, lh_flux_1p, lh_err_1p, lh_bit_1p, format = 'i,f,f,f,i', /silent


  		readcol, specdir_nod+'cso005'+'_sl1_2p_cal.tbl', $
 			sl1_order_2p, sl1_wave_2p, sl1_flux_2p, sl1_err_2p, sl1_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_sl2_2p_cal.tbl', $
 			sl2_order_2p, sl2_wave_2p, sl2_flux_2p, sl2_err_2p, sl2_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_sl3_2p_cal.tbl', $
 			sl3_order_2p, sl3_wave_2p, sl3_flux_2p, sl3_err_2p, sl3_bit_2p, format = 'i,f,f,f,i', /silent
 	
 		readcol, specdir_nod+'cso005'+'_ll1_2p_cal.tbl', $
 			ll1_order_2p, ll1_wave_2p, ll1_flux_2p, ll1_err_2p, ll1_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_ll2_2p_cal.tbl', $
 			ll2_order_2p, ll2_wave_2p, ll2_flux_2p, ll2_err_2p, ll2_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+'cso005'+'_ll3_2p_cal.tbl', $
 			ll3_order_2p, ll3_wave_2p, ll3_flux_2p, ll3_err_2p, ll3_bit_2p, format = 'i,f,f,f,i', /silent
 	
 	
 		readcol, specdir_coadd+'cso005'+'_sh_cal.tbl', $
 			sh_order, sh_wave, sh_flux, sh_err, sh_bit, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_coadd+'cso005'+'_lh_cal.tbl', $
 			lh_order, lh_wave, lh_flux, lh_err, lh_bit, format = 'i,f,f,f,i', /silent
 	

  		readcol, specdir_coadd+'cso005'+'_sl1_cal.tbl', $
 			sl1_order, sl1_wave, sl1_flux, sl1_err, sl1_bit, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_coadd+'cso005'+'_sl2_cal.tbl', $
 			sl2_order, sl2_wave, sl2_flux, sl2_err, sl2_bit, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_coadd+'cso005'+'_sl3_cal.tbl', $
 			sl3_order, sl3_wave, sl3_flux, sl3_err, sl3_bit, format = 'i,f,f,f,i', /silent
 	
 		readcol, specdir_coadd+'cso005'+'_ll1_cal.tbl', $
 			ll1_order, ll1_wave, ll1_flux, ll1_err, ll1_bit, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_coadd+'cso005'+'_ll2_cal.tbl', $
 			ll2_order, ll2_wave, ll2_flux, ll2_err, ll2_bit, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_coadd+'cso005'+'_ll3_cal.tbl', $
 			ll3_order, ll3_wave, ll3_flux, ll3_err, ll3_bit, format = 'i,f,f,f,i', /silent
 	
 	
 		readcol, specdir_coadd+'cso005'+'_sh_cal.tbl', $
 			sh_order, sh_wave, sh_flux, sh_err, sh_bit, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_coadd+'cso005'+'_lh_cal.tbl', $
 			lh_order, lh_wave, lh_flux, lh_err, lh_bit, format = 'i,f,f,f,i', /silent
 	

		flux_sh = sh_flux_1p - sh_flux_2p
		wave_sh = sh_wave_2p / (1d + sed.redshift)

		wave_sl = [sl2_wave_1p, sl1_wave_1p] / (1d + sed.redshift)
		
		ymin = min(flux_sh)
		ymax = max(flux_sh) * 1.1

		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/cso005.sav'
		
		flux_lh = sed.flux_lh
		wave_lh = sed.wave_lh 
		order_lh = sed.order_lh
		
		noneg_lh = where(sed.flux_lh gt 0)
		flux_lh = flux_lh(noneg_lh)
		wave_lh = wave_lh(noneg_lh)
		order_lh = order_lh(noneg_lh)

; Nods

plot, wave_sh, sh_flux_1p, psym=10, xr=[10,11], yr=[0,0.3], thick=2, xthick=2, ythick=2, charthick=2, charsize=2, xtitle = '!7k!3!Iobs!N [!7l!3m]', ytitle = 'Flux density [Jy]', title=sed.obj, /xstyle, /ystyle
oplot, wave_sh, sh_flux_2p + 0.2, psym=10, thick=2 
ver, 10.511, linestyle=1
xyouts,10.1,0.15,'SIV',/data, charsize=1.0

; NeII

plot, wave_sh, sh_flux_1p, psym=10, xr=[12.3,13.3], yr=[0,0.3], thick=2, xthick=2, ythick=2, charthick=2, charsize=2, xtitle = '!7k!3!Iobs!N [!7l!3m]', ytitle = 'Flux density [Jy]', title=sed.obj, /xstyle, /ystyle
oplot, wave_sh, sh_flux_2p + 0.2, psym=10, thick=2 
ver, 12.814, linestyle=1
xyouts,12.5,0.15,'NeII',/data, charsize=1.0

; 6.2 PAH

plot, wave_sl, [sl2_flux_1p,sl1_flux_1p], psym=10, xr=[4,10], yr=[1d-3,1d0], thick=2, xthick=2, ythick=2, charthick=2, charsize=2, xtitle = '!7k!3!Iobs!N [!7l!3m]', ytitle = 'Flux density [Jy]', title=sed.obj, /xlog, /ylog, /xstyle, /ystyle
oplot, wave_sl, [sl2_flux_2p,sl1_flux_2p] *4, psym=10, thick=2 
ver, 6.2, linestyle=1
xyouts, 4.3, 3d-1, 'PAH 6.2', /data, charsize=1.0

; Diff spec. 

plot, wave_sh, flux_sh, $
	xr = [10,11], /xstyle, $
	yrange = [-1 * ymax,ymax], /ystyle, $
	xtitle = '!7k!3!Iobs!N [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = 2, $
	xthick = 2, $
	ythick = 2, $
	thick = lthick, $
	charthick = lthick, $
	/nodata

oplot, wave_sh,flux_sh, psym = 10, thick = lthick
	
ver, 10.511, linestyle=1
xyouts,10.1,0.015,'SIV',/data, charsize=1.0


plot, wave_sh, flux_sh, $
	xr = [12.3,13.3], /xstyle, $
	yrange = [-1 * ymax,ymax], /ystyle, $
	xtitle = '!7k!3!Iobs!N [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = 2, $
	xthick = 2, $
	ythick = 2, $
	thick = lthick, $
	charthick = lthick, $
	/nodata

oplot, wave_sh,flux_sh, psym = 10, thick = lthick
	
ver, 12.814, linestyle=1
xyouts,12.5,0.015,'NeII',/data, charsize=1.0


plot, wave_sh, flux_sh, $
	xr = [4,10], /xstyle, $
	yrange = [-1 * ymax,ymax], /ystyle, $
	xtitle = '!7k!3!Iobs!N [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = sed.obj, $
	charsize = 2, $
	xthick = 2, $
	ythick = 2, $
	thick = lthick, $
	charthick = lthick, $
	/xlog, $
	/nodata

oplot, wave_sl,[sl2_flux_1p - sl2_flux_2p, sl1_flux_1p - sl1_flux_2p], psym = 10, thick = lthick
	
ver, 6.2, linestyle=1
xyouts, 4.3, 0.01, 'PAH 6.2', /data, charsize=1.0


; Co-added

plot, wave_sh, sh_flux, psym=10, xr=[10,11], yr=[0,0.05], thick=2, xthick=2, ythick=2, charthick=2, charsize=2, xtitle = '!7k!3!Iobs!N [!7l!3m]', ytitle = 'Flux density [Jy]', title=sed.obj, /xstyle, /ystyle
;fullql_hr,'cso005',xr=[10,11], /panel, /bw, /nolines, yr=[0,0.1]
ver, 10.511, linestyle=1
xyouts,10.1,0.015,'SIV',/data, charsize=1.0

plot, wave_sh, sh_flux, psym=10, xr=[12.3,13.3], yr=[0,0.1], thick=2, xthick=2, ythick=2, charthick=2, charsize=2, xtitle = '!7k!3!Iobs!N [!7l!3m]', ytitle = 'Flux density [Jy]', title=sed.obj, /xstyle, /ystyle
;fullql_hr,'cso005',xr=[12.3,13.3], /panel, /bw, /nolines, yr=[0,0.1]
ver, 12.814, linestyle=1
xyouts,12.5,0.04,'NeII',/data, charsize=1.0

plot, wave_sl, [sl2_flux,sl1_flux], psym=10, xr=[4,10], yr=[1d-3,1d-1], thick=2, xthick=2, ythick=2, charthick=2, charsize=2, xtitle = '!7k!3!Iobs!N [!7l!3m]', ytitle = 'Flux density [Jy]', title=sed.obj, /xlog, /ylog, /xstyle, /ystyle
;spec2pahfit,'cso005',xr=[4,10], /panel, /bw, /nolines, yr=[1d-3,1d-1], /nolabel, /nobonus
ver, 6.2, linestyle=1
xyouts, 4.3, 3d-2, 'PAH 6.2', /data, charsize=1.0

device, /close
set_plot,'x'
end
