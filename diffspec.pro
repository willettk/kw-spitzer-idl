; pro diffspec
;+
; NAME:
;       
;	DIFFSPEC
;
; PURPOSE:
;
;	Plot a difference spectrum of two nods for IRS data
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
;       Written by K. Willett                Jan 09
;-

; Load SH data


csotag = csodat('tag')
ncso = n_elements(csotag)
ncso = 1

!p.multi=[0,1,ncso]

;for m = 0,1 do begin

;	for j = 3*m,3*m + 2 do begin
	for j = 0,ncso-1 do begin
	
		tag, csotag[j], dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+csotag[j]+'.sav'

		specdir_nod = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/nods/'
		specdir_coadd = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/calibrated/coadd/'
 	
 		readcol, specdir_nod+csotag[j]+'_sl1_1p_cal.tbl', $
 			sl1_order_1p, sl1_wave_1p, sl1_flux_1p, sl1_err_1p, sl1_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_sl2_1p_cal.tbl', $
 			sl2_order_1p, sl2_wave_1p, sl2_flux_1p, sl2_err_1p, sl2_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_sl3_1p_cal.tbl', $
 			sl3_order_1p, sl3_wave_1p, sl3_flux_1p, sl3_err_1p, sl3_bit_1p, format = 'i,f,f,f,i', /silent
 	
 		readcol, specdir_nod+csotag[j]+'_ll1_1p_cal.tbl', $
 			ll1_order_1p, ll1_wave_1p, ll1_flux_1p, ll1_err_1p, ll1_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_ll2_1p_cal.tbl', $
 			ll2_order_1p, ll2_wave_1p, ll2_flux_1p, ll2_err_1p, ll2_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_ll3_1p_cal.tbl', $
 			ll3_order_1p, ll3_wave_1p, ll3_flux_1p, ll3_err_1p, ll3_bit_1p, format = 'i,f,f,f,i', /silent
 	
 	
 		readcol, specdir_nod+csotag[j]+'_sh_1p_cal.tbl', $
 			sh_order_1p, sh_wave_1p, sh_flux_1p, sh_err_1p, sh_bit_1p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_lh_1p_cal.tbl', $
 			lh_order_1p, lh_wave_1p, lh_flux_1p, lh_err_1p, lh_bit_1p, format = 'i,f,f,f,i', /silent



  		readcol, specdir_nod+csotag[j]+'_sl1_2p_cal.tbl', $
 			sl1_order_2p, sl1_wave_2p, sl1_flux_2p, sl1_err_2p, sl1_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_sl2_2p_cal.tbl', $
 			sl2_order_2p, sl2_wave_2p, sl2_flux_2p, sl2_err_2p, sl2_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_sl3_2p_cal.tbl', $
 			sl3_order_2p, sl3_wave_2p, sl3_flux_2p, sl3_err_2p, sl3_bit_2p, format = 'i,f,f,f,i', /silent
 	
 		readcol, specdir_nod+csotag[j]+'_ll1_2p_cal.tbl', $
 			ll1_order_2p, ll1_wave_2p, ll1_flux_2p, ll1_err_2p, ll1_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_ll2_2p_cal.tbl', $
 			ll2_order_2p, ll2_wave_2p, ll2_flux_2p, ll2_err_2p, ll2_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_ll3_2p_cal.tbl', $
 			ll3_order_2p, ll3_wave_2p, ll3_flux_2p, ll3_err_2p, ll3_bit_2p, format = 'i,f,f,f,i', /silent
 	
 	
 		readcol, specdir_nod+csotag[j]+'_sh_2p_cal.tbl', $
 			sh_order_2p, sh_wave_2p, sh_flux_2p, sh_err_2p, sh_bit_2p, format = 'i,f,f,f,i', /silent
 		
 		readcol, specdir_nod+csotag[j]+'_lh_2p_cal.tbl', $
 			lh_order_2p, lh_wave_2p, lh_flux_2p, lh_err_2p, lh_bit_2p, format = 'i,f,f,f,i', /silent
 	

		flux_sh = sh_flux_1p - sh_flux_2p
		wave_sh = sh_wave_2p
		
;		noneg_sh = where(flux_sh gt 0)
;		flux_sh = flux_sh[noneg_sh]
;		wave_sh = wave_sh[noneg_sh]
		
		; Identified IR lines for overplotting
		
		templines = ir_lines(/hr)
		lines = templines(*,0) & line_id = templines(*,1)
		
		; Plot data
		
		; Hard copy option
		
		lthick = 1
		cs = 2
		
		ymin = min(flux_sh)
		ymax = max(flux_sh) * 1.1

		plot, wave_sh, flux_sh, $
;			xrange = [floor(min(wave_sh)),ceil(max(wave_sh))], /xstyle, $
			yrange = [-1 * ymax,ymax], /ystyle, $
			xtitle = 'Wavelength (observed frame) [!7l!3m]', $
			ytitle = 'Flux density [Jy]', $
			title = sed.tag+' '+sed.obj, $
			charsize = 2, $
			xthick = 2, $
			ythick = 2, $
			thick = lthick, $
			charthick = lthick, $
			/nodata
	
		oplot, wave_sh,flux_sh, psym = 10, thick = lthick
		
		; Instead of vertical lines, create bars with fine-structure positions
	
;		fslines = [8.99138d,10.5110d,12.368d,12.8140d,14.3220d,14.368d,15.555d,17.9340d,18.7130d]
;		fsid = ['ArIII','SIV','Hu-!7a!3','NeII','NeV','ClII','NeIII','FeII','SIII']
;		fsind = where(fslines lt max(wave_sh) and fslines gt min(wave_sh))
;		fslines = fslines[fsind] & fsid = fsid[fsind]
;		bheight = (!y.crange[0] + 0.15*(!y.crange[1] - !y.crange[0]))
;		blevel = (!y.crange[0]+0.1*(!y.crange[1] - !y.crange[0]))
;		blabel = (!y.crange[0]+0.05*(!y.crange[1] - !y.crange[0]))
;		blabellev = (!y.crange[0]+0.12*(!y.crange[1] - !y.crange[0]))
;	
;		plots, [fslines(0),fslines(n_elements(fslines)-1)], [blevel,blevel], thick = lthick
;		for k=0,n_elements(fslines)-1 do begin
;			plots, [fslines(k),fslines(k)], [blevel,bheight], thick = lthick
;			if fsid[k] eq 'Hu-!7a!3' then xyouts,fslines[k]-0.4,blabellev,fsid[k],charsize=2, charthick=lthick $
;				else if fsid[k] eq 'NeV' or fsid[k] eq 'FeII' $
;				then xyouts,fslines[k]-0.3,blabellev,fsid[k],charsize=2, charthick=lthick $
;				else xyouts,fslines[k]+0.1,blabellev,fsid[k],charsize=2, charthick=lthick
;		endfor
;;		xyouts,12.0,blabel,'Atomic transitions',charsize=0.5
;	
;		h2lines = [9.665d,12.279d,17.035d,28.218d]
;		h2id = ['S(3)','S(2)','S(1)','S(0)']
;		h2ind = where(h2lines lt max(wave_sh))
;		h2lines = h2lines[h2ind] & h2id = h2id[h2ind]
;		bheight = (!y.crange[1]-0.15*(!y.crange[1] - !y.crange[0]))
;		blevel = (!y.crange[0]+0.90*(!y.crange[1] - !y.crange[0]))
;		blabel = (!y.crange[0]+0.92*(!y.crange[1] - !y.crange[0]))
;		blabellev = (!y.crange[0]+0.83*(!y.crange[1] - !y.crange[0]))
;	
;		plots, [h2lines(0),h2lines(n_elements(h2lines)-1)], [blevel,blevel], thick = lthick
;		for k=0,n_elements(h2lines)-1 do begin
;			plots, [h2lines(k),h2lines(k)], [blevel, bheight], thick = lthick
;			xyouts,h2lines[k]+0.1,blabellev,h2id[k],charsize=2, charthick=lthick
;		endfor
;		xyouts,mean([h2lines[0],h2lines[n_elements(h2lines)-1]]),blabel,'H!I2!N',charsize=2, charthick=lthick
	
	endfor	
	


end
