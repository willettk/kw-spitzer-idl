pro bisurv_line, line, lh = lh
;+
; NAME:
;       
;	BISURV_LINE
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

tags = [transpose(ohmdat('tag')),transpose(archdat('tag'))]
objs = [transpose(ohmdat('obj')),transpose(archdat('obj'))]
dl   = [transpose(ohmdat('dl')),transpose(archdat('dl'))]

tags_con = [transpose(condat('tag'))]
objs_con = [transpose(condat('obj'))]
dl_con = [transpose(condat('dl'))]

; Load OH data

logoh_mega = transpose(float(ohmdat('logoh')))
logoh_arch = transpose(float(archdat('logoh')))

logoh_ohm = [logoh_mega, logoh_arch]

logoh_con = transpose(float(condat('logoh')))

f1667_mega = transpose(float(ohmdat('f1667')))
f1667_arch = transpose(float(archdat('f1667')))

f1667_ohm = [f1667_mega, f1667_arch]

f1667_con = transpose(float(condat('f1667')))

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

	ohm_line     = transpose(alog10(linearr[1,*] * dl_det^2 * 1d-21 * 1d7 / 3.826d33))
	ohm_line_lim = transpose(alog10(linelimarr[1,*] * dl_lim^2 * 1d-21 * 1d7 / 3.826d33))
	con_line     = transpose(alog10(linearr_con[1,*] * dl_det_con^2 * 1d-21 * 1d7 / 3.826d33))
	con_line_lim = transpose(alog10(linelimarr_con[1,*] * dl_lim_con^2 * 1d-21 * 1d7 / 3.826d33))

	match, tags_det, tags, od1, od2
	match, tags_lim, tags, ol1, ol2

	match, tags_det_con, tags_con, cd1, cd2
	match, tags_lim_con, tags_con, cl1, cl2

censor = [replicate('0',ndet),$
	replicate('-1', nlim),$
	replicate('-2', ndet_con),$
	replicate('-3', nlim_con)]

ind_data = string($
	logoh_ohm[od2], $
	logoh_ohm[ol2], $
	logoh_con[cd2], $
	logoh_con[cl2], $
	format='(f10.2)')

dep_data = string($
	ohm_line, $
	ohm_line_lim, $
	con_line, $
	con_line_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent data
; Col 3 - dependent data

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_'+line+'.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

end
