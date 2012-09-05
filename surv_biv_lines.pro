
;+
; NAME:
;       
;	SURV_BIV_LINES
;
; PURPOSE:
;
;	Batch file for line luminosities for bivariate survival analysis
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
;       Written by K. Willett                
;-


lines_sh = ['sIV','hI76','neII','neV','neIII','h2s3','h2s2','arIII']
lines_lh = ['neV24','oIV','feII26','h2s0']

both_lines = ['feII','sIII','h2s1']
outerlim_lines = ['sIII33','siII34']


for i = 0, n_elements(lines_sh) - 1 do begin

	bisurv_line, lines_sh[i]
	print,'Did bivariate SH '+lines_sh[i]

endfor

for i = 0, n_elements(lines_lh) - 1 do begin

	bisurv_line, lines_lh[i], /lh
	print,'Did bivariate LH '+lines_lh[i]

endfor

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


; FeII (18 um)

line = 'feII'

linelist, line, linearr, /quiet
linelimlist, line, linelimarr, /quiet, /lh
linelimlist, line, linelimarr2, /quiet

linelist, line, linearr_con, /quiet, /con
linelimlist, line, linelimarr_con, /quiet, /con, /lh
linelimlist, line, linelimarr_con2, /quiet, /con

ndet = n_elements(linearr[0,*])
nlim = n_elements(linelimarr[0,*])
nlim2 = n_elements(linelimarr2[0,*])

ndet_con = n_elements(linearr_con[0,*])
nlim_con = n_elements(linelimarr_con[0,*])
nlim_con2 = n_elements(linelimarr_con2[0,*])

match, transpose(linelimarr[0,*]), transpose(linelimarr2[0,*]), ind1, ind2

linelimarr[1,ind1] = float(linelimarr[1,ind1]) > float(linelimarr2[1,ind2])

linelimarr_con = [[linelimarr_con],[linelimarr_con2]]
nlim_con = nlim_con + nlim_con2

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

print,'Did bivariate feII'

;;;;;;;;;;;;;
;
; SIII (18.7 um)

line = 'sIII'

linelist, line, linearr, /quiet
linelimlist, line, linelimarr, /quiet, /lh
linelimlist, line, linelimarr2, /quiet

linelist, line, linearr_con, /quiet, /con
linelimlist, line, linelimarr_con, /quiet, /con, /lh

a = size(linelimarr)
b = size(linelimarr_con)

if a[1] or b[1] eq 1 then message, 'No limits found for line in this module'

ndet = n_elements(linearr[0,*])
nlim = n_elements(linelimarr[0,*])
nlim2 = n_elements(linelimarr2[0,*])

ndet_con = n_elements(linearr_con[0,*])
nlim_con = n_elements(linelimarr_con[0,*])

match, transpose(linelimarr[0,*]), transpose(linelimarr2[0,*]), ind1, ind2

linelimarr[1,ind1] = float(linelimarr[1,ind1]) > float(linelimarr2[1,ind2])

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

print,'Did bivariate SIII'

;;;;;;;;;;;;

; H2 S1 (18.7 um)

line = 'h2s1'

linelist, line, linearr, /quiet
linelimlist, line, linelimarr, /quiet, /lh

linelist, line, linearr_con, /quiet, /con
linelimlist, line, linelimarr_con, /quiet, /con, /lh
linelimlist, line, linelimarr_con2, /quiet, /con

a = size(linelimarr)
b = size(linelimarr_con)

if a[1] or b[1] eq 1 then message, 'No limits found for line in this module'

ndet = n_elements(linearr[0,*])
nlim = n_elements(linelimarr[0,*])

ndet_con = n_elements(linearr_con[0,*])
nlim_con = n_elements(linelimarr_con[0,*])

linelimarr_con = [[linelimarr_con],[linelimarr_con2]]
nlim_con = 2

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

print,'Did bivariate h2s1'

;;;;;;;;;;;;

; SIII33 and SiII34

for j = 0,1 do begin

	line = outerlim_lines[j]

	linelist, line, linearr, /quiet
	linelimlist, line, linelimarr, /quiet, /lh
	
	linelist, line, linearr_con, /quiet, /con
	
	a = size(linelimarr)
	
	if a[1] eq 1 then message, 'No limits found for line in this module'
	
	ndet = n_elements(linearr[0,*])
	nlim = n_elements(linelimarr[0,*])
	
	ndet_con = n_elements(linearr_con[0,*])
	
	tags_det = transpose(linearr[0,*])
	tags_lim = transpose(linelimarr[0,*])
	
	tags_det_con = transpose(linearr_con[0,*])
	
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
	
	ohm_line     = transpose(alog10(linearr[1,*] * dl_det^2 * 1d-21 * 1d7 / 3.826d33))
	ohm_line_lim = transpose(alog10(linelimarr[1,*] * dl_lim^2 * 1d-21 * 1d7 / 3.826d33))
	con_line     = transpose(alog10(linearr_con[1,*] * dl_det_con^2 * 1d-21 * 1d7 / 3.826d33))

	match, tags_det, tags, od1, od2
	match, tags_lim, tags, ol1, ol2

	match, tags_det_con, tags_con, cd1, cd2

	censor = [replicate('0',ndet),$
		replicate('-1', nlim),$
		replicate('-2', ndet_con)]
	
	ind_data = string($
		logoh_ohm[od2], $
		logoh_ohm[ol2], $
		logoh_con[cd2], $
		format='(f10.2)')
	
	dep_data = string($
		ohm_line, $
		ohm_line_lim, $
		con_line, $
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
	
	
	print,'Did bivariate '+line
	
endfor

;;;;;;;;;;;;

; Absorption features

for j = 0,2 do begin

	line = absorption[j]

	linelist, line, linearr, /quiet
	linelimlist, line, linelimarr, /quiet
	
	linelimlist, line, linelimarr_con, /quiet, /con
	
	a = size(linelimarr)
	
	if a[1] eq 1 then message, 'No limits found for line in this module'
	
	ndet = n_elements(linearr[0,*])
	nlim = n_elements(linelimarr[0,*])
	
	nlim_con = n_elements(linelimarr_con[0,*])
	
	tags_det = transpose(linearr[0,*])
	tags_lim = transpose(linelimarr[0,*])
	
	tags_lim_con = transpose(linelimarr_con[0,*])
	
	dl_det = dblarr(ndet) 
	dl_lim = dblarr(nlim) & dl_lim_con = dblarr(nlim_con)
	
	for i = 0, ndet-1 do begin
		targets, tags_det[i], z, obj, dl
		dl_det[i] = dl * 3.0826d24
	endfor
	
	for i = 0, nlim-1 do begin
		targets, tags_lim[i], z, obj, dl
		dl_lim[i] = dl * 3.0826d24
	endfor
	
	for i = 0, nlim_con-1 do begin
		targets, tags_lim_con[i], z, obj, dl
		dl_lim_con[i] = dl * 3.0826d24
	endfor
	
	ohm_line     = transpose(alog10(-1d * linearr[1,*] * dl_det^2 * 1d-21 * 1d7 / 3.826d33))
	ohm_line_lim = transpose(alog10(linelimarr[1,*] * dl_lim^2 * 1d-21 * 1d7 / 3.826d33))
	con_line_lim = transpose(alog10(linelimarr_con[1,*] * dl_lim_con^2 * 1d-21 * 1d7 / 3.826d33))

	match, tags_det, tags, od1, od2
	match, tags_lim, tags, ol1, ol2

	match, tags_lim_con, tags_con, cl1, cl2

	censor = [replicate('0',ndet),$
		replicate('-1', nlim),$
		replicate('-3', nlim_con)]
	
	ind_data = string($
		logoh_ohm[od2], $
		logoh_ohm[ol2], $
		logoh_con[cl2], $
		format='(f10.2)')
	
	dep_data = string($
		ohm_line, $
		ohm_line_lim, $
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
	
	
	print,'Did bivariate absorption '+line
	

endfor

;;;;;;;;;;;;

;;;;;;;;;;;;


end
