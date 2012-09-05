;+
; NAME:
;       
;	SURV_BIVARIATE
;
; PURPOSE:
;
;	Batch file to make data for bivariate survival analysis for Spitzer OHMs/non-masers
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
;       Written by K. Willett                Feb 08
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

;;;;;;;;;;;;;;;
; Silicate
;;;;;;;;;;;;;;;

silicate_mega = transpose(-1d * float(ohmdat('sil')))
silicate_arch = transpose(-1d * float(archdat('sil')))

silicate_con = transpose(-1d * float(condat('sil')))

silicate_ohm = [silicate_mega[*,0], silicate_arch[*,0]]
silicate_con = silicate_con[*,0]

ohm_index = where(silicate_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(silicate_con ne 0. and logoh_con ne 0.)

ohm_dep = silicate_ohm[ohm_index]
con_dep = silicate_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-1',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(dep_data), transpose(ind_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_regression/bireg_data_silicate.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; Spectral slopes
;;;;;;;;;;;;;;;

slope_mega = transpose(-1d * float(ohmdat('spindex')))
slope_arch = transpose(-1d * float(archdat('spindex')))

slope_con = transpose(-1d * float(condat('spindex')))

slope30_ohm = [slope_mega[*,1], slope_arch[*,1]]
slope30_con = slope_con[*,1]

ohm_index = where(slope30_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(slope30_con ne 0. and logoh_con ne 0.)

ohm_dep = slope30_ohm[ohm_index]
con_dep = slope30_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-1',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(dep_data), transpose(ind_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_regression/bireg_data_slope30.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

end
