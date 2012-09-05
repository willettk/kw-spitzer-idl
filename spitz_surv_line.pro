pro spitz_surv_line, line, lh = lh
;+
; NAME:
;       
;	SPITZ_SURV
;
; PURPOSE:
;
;	
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
;       Written by K. Willett                Feb 09
;-

	; HR line luminosities
	
	linelist, line, linearr, /quiet
	linelimlist, line, linelimarr, /quiet, lh = lh
	
	linelist, line, linearr_con, /quiet, /con
	linelimlist, line, linelimarr_con, /quiet, /con, lh = lh
	
	a = size(linelimarr)
	b = size(linelimarr_con)
	
	if a[1] or b[1] eq 1 then message, 'No limits found for line in this module'
	
	ndet = n_elements(linearr[0,*])
	nlim = n_elements(linelimarr[0,*])
	
	ndet_con = n_elements(linearr_con[0,*])
	nlim_con = n_elements(linelimarr_con[0,*])
	
	tags_det = transpose(linearr[0,*])
	tags_lim = transpose(linelimarr[0,*])
	
	tags_det_con = transpose(linearr_con[0,*])
	tags_lim_con = transpose(linelimarr_con[0,*])
	
	dl_lim = dblarr(nlim) & dl_lim_con = dblarr(nlim_con)
	dl_det = dblarr(ndet) & dl_det_con = dblarr(ndet_con)
	
	for i = 0, ndet-1 do begin
		targets, tags_det[i], z, obj, dl
		dl_det[i] = dl * 3.0826d24
	endfor
	
	for i = 0, nlim-1 do begin
		targets, tags_lim[i], z, obj, dl
		dl_lim[i] = dl * 3.0826d24
	endfor
	
	for i = 0, ndet_con-1 do begin
		targets, tags_det_con[i], z, obj, dl
		dl_det_con[i] = dl * 3.0826d24
	endfor
	
	for i = 0, nlim_con-1 do begin
		targets, tags_lim_con[i], z, obj, dl
		dl_lim_con[i] = dl * 3.0826d24
	endfor


censor = [replicate(' 0',ndet),replicate('-1',nlim),$
	replicate(' 0',ndet_con),replicate('-1',nlim_con)]

data = string($
	alog10(([transpose(linearr[1,*] * dl_det^2), $
	transpose(linelimarr[1,*] * dl_lim^2), $
	transpose(linearr_con[1,*] * dl_det_con^2), $
	transpose(linelimarr_con[1,*] * dl_lim_con^2)]) * 1d-21 * 1d7 / 3.826d33),$
	format='(f10.2)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

; Col 1 - censor
; Col 2 - data
; Col 3 - group

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_'+line+'.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

end
