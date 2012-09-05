
;+
; NAME:
;       
;	SURV_BIVARIATE
;
; PURPOSE:
;
;	Batch file to make data for bivarate survival analysis for Spitzer OHMs/non-masers
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
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_silicate.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; log L_FIR
;;;;;;;;;;;;;;;

lfir_mega = transpose(float(ohmdat('lfir')))
lfir_arch = transpose(float(archdat('lfir')))

lfir_ohm = [lfir_mega, lfir_arch]

lfir_con = transpose(float(condat('lfir')))

ohm_index = where(lfir_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(lfir_con ne 0. and logoh_con ne 0.)

ohm_dep = lfir_ohm[ohm_index]
con_dep = lfir_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_lfir.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; Dust temperature
;;;;;;;;;;;;;;;

dtemp_mega = transpose(float(ohmdat('dtemp')))
dtemp_arch = transpose(float(archdat('dtemp')))

dtemp_ohm = [dtemp_mega, dtemp_arch]

dtemp_con = transpose(float(condat('dtemp')))

ohm_index = where(dtemp_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(dtemp_con ne 0. and logoh_con ne 0.)

ohm_dep = dtemp_ohm[ohm_index]
con_dep = dtemp_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_dtemp.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; H2 mass and temperature
;;;;;;;;;;;;;;;

restore,'~/Astronomy/Research/Spitzer/ohm/h2.sav'

h2_temp_warm_ohm = h2.temp_warm

restore,'~/Astronomy/Research/Spitzer/control/h2_con.sav'

h2_temp_warm_con = h2.temp_warm

ohm_index = where(h2_temp_warm_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(h2_temp_warm_con ne 0. and logoh_con ne 0.)

ohm_dep = h2_temp_warm_ohm[ohm_index]
con_dep = h2_temp_warm_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_h2_temp_warm.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

restore,'~/Astronomy/Research/Spitzer/ohm/h2.sav'

h2_mass_warm_ohm = h2.mass_warm

restore,'~/Astronomy/Research/Spitzer/control/h2_con.sav'

h2_mass_warm_con = h2.mass_warm

ohm_index = where(h2_mass_warm_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(h2_mass_warm_con ne 0. and logoh_con ne 0.)

ohm_dep = h2_mass_warm_ohm[ohm_index]
con_dep = h2_mass_warm_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_h2_mass_warm.txt'

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

slope15_ohm = [slope_mega[*,0], slope_arch[*,0]]
slope15_con = slope_con[*,0]

slope30_ohm = [slope_mega[*,1], slope_arch[*,1]]
slope30_con = slope_con[*,1]

ohm_index = where(slope15_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(slope15_con ne 0. and logoh_con ne 0.)

ohm_dep = slope15_ohm[ohm_index]
con_dep = slope15_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_slope15.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

ohm_index = where(slope30_ohm ne 0. and logoh_ohm ne 0.)
con_index = where(slope30_con ne 0. and logoh_con ne 0.)

ohm_dep = slope30_ohm[ohm_index]
con_dep = slope30_con[con_index]

ohm_ind = logoh_ohm[ohm_index]
con_ind = logoh_con[con_index]

ndep = n_elements(ohm_dep) & ndep_con = n_elements(con_dep)

censor = [replicate(' 0',ndep),$				
	replicate('-2',ndep_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	format='(f6.2)')

dep_data = string($
	ohm_dep, $
	con_dep, $
	format='(f6.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_slope30.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; f60/f100
;;;;;;;;;;;;;;;

iras_mega = ohmdat('iras')
f60_mega  = transpose(float(iras_mega[2,*]))
f100_mega = transpose(float(iras_mega[3,*]))

iras_arch = archdat('iras')
f60_arch  = transpose(float(iras_arch[2,*]))
f100_arch = transpose(float(iras_arch[3,*]))

f60_ohm = [f60_mega, f60_arch]
f100_ohm = [f100_mega, f100_arch]

iras_con = condat('iras')
f60_con  = transpose(float(iras_con[2,*]))
f100_con = transpose(float(iras_con[3,*]))

; Upper limits on data from NED

;lim100_ohm = [1,2,21,23,25,50]		<--- Indices of data w/100 um flux limits
;lim100_con = [1,8]

lim100_ohm = [4.789,2.145,1.547,1.033,1.373,1.633]
lim100_con = [0.9896,1.453]

ohm_detind_1667 = where(f100_ohm ne 0 and f1667_ohm ne 0)
con_detind_1667 = where(f100_con ne 0 and f1667_con ne 0)

ohm_limind_1667 = where(f100_ohm eq 0 and f1667_ohm ne 0)
con_limind_1667 = where(f100_con eq 0 and f1667_con ne 0)

ohm_det = f60_ohm[ohm_detind_1667] / f100_ohm[ohm_detind_1667]
con_det = f60_con[con_detind_1667] / f100_con[con_detind_1667]

ohm_lim = f60_ohm[ohm_limind_1667] / lim100_ohm
con_lim = f60_con[con_limind_1667] / lim100_con

f1667_ohm_fir     = f1667_ohm[ohm_detind_1667]
f1667_con_fir     = f1667_con[con_detind_1667]
f1667_ohm_fir_lim = f1667_ohm[ohm_limind_1667]
f1667_con_fir_lim = f1667_con[con_limind_1667]

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate('0',ndet),$
	replicate('1',nlim),$
	replicate('-2',ndet_con),$
	replicate('4',nlim_con)]

ind_data = string($
	f1667_ohm_fir, $
	f1667_ohm_fir_lim, $
	f1667_con_fir, $
	f1667_con_fir_lim, $
	format='(f10.2)')

dep_data = string($
	ohm_det, $
	ohm_lim, $
	con_det, $
	con_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent variable
; Col 3 - dependent variable

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_f60_f100.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; L_PAH and EW PAH
;;;;;;;;;;;;;;;

pah62_mega = transpose(float(ohmdat('pah62lum')))
pah62_arch = transpose(float(archdat('pah62lum')))

pah62_ohm = [pah62_mega, pah62_arch]

pah62_con = transpose(float(condat('pah62lum')))

pah11_mega = transpose(float(ohmdat('pah11lum')))
pah11_arch = transpose(float(archdat('pah11lum')))

pah11_ohm = [pah11_mega, pah11_arch]

pah11_con = transpose(float(condat('pah11lum')))

; 

ewpah62_mega = transpose(float(ohmdat('pah62ew')))
ewpah62_arch = transpose(float(archdat('pah62ew')))

ewpah62_ohm = [ewpah62_mega, ewpah62_arch]

ewpah62_con = transpose(float(condat('pah62ew')))

ewpah11_mega = transpose(float(ohmdat('pah11ew')))
ewpah11_arch = transpose(float(archdat('pah11ew')))

ewpah11_ohm = [ewpah11_mega, ewpah11_arch]

ewpah11_con = transpose(float(condat('pah11ew')))

; 

index_ohm = setdifference(indgen(n_elements(pah62_ohm)),[5,40,25])
index_con = setdifference(indgen(n_elements(pah62_con)),[1,8])

; Upper limits on data from NED

;lim62_ohm = [5,40]		;      <--- Indices of data w/limits
;lim62_con = [1]
;lim11_ohm = []
;lim11_con = [1]

pah62_ohm_lim = [8.95, 9.79]
pah62_con_lim = [9.39]

ewpah62_ohm_lim = [0.10, 0.01]
ewpah62_con_lim = [0.01]

pah62_ohm_det = pah62_ohm[index_ohm] 
pah62_con_det = pah62_con[index_con]

ewpah62_ohm_det = ewpah62_ohm[index_ohm] 
ewpah62_con_det = ewpah62_con[index_con]

ohm_ind = logoh_ohm[index_ohm]
con_ind = logoh_con[index_con]

ohm_ind_lim = logoh_ohm[[5,40]]
con_ind_lim = logoh_con[1]

ndet = n_elements(ohm_ind) & ndet_con = n_elements(con_ind)
nlim = n_elements(ohm_ind_lim) & nlim_con = n_elements(con_ind_lim)

censor = [replicate('0',ndet),$
	replicate('-1', nlim),$
	replicate('-2', ndet_con),$
	replicate('-3', nlim_con)]

ind_data = string($
	ohm_ind, $
	ohm_ind_lim, $
	con_ind, $
	con_ind_lim, $
	format='(f10.2)')

dep_data = string($
	pah62_ohm_det, $
	pah62_ohm_lim, $
	pah62_con_det, $
	pah62_con_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent variable
; Col 3 - dependent variable

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_pah62lum.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

dep_data = string($
	ewpah62_ohm_det, $
	ewpah62_ohm_lim, $
	ewpah62_con_det, $
	ewpah62_con_lim, $
	format='(f10.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_pah62ew.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

; Upper limits on data from NED

;lim11_ohm = []
;lim11_con = [1]

index_ohm = indgen(n_elements(pah11_ohm))
index_con = setdifference(indgen(n_elements(pah11_con)),[1,8])

pah11_con_lim = [9.10]

ewpah11_con_lim = [0.01]

pah11_ohm_det = pah11_ohm[index_ohm] 
pah11_con_det = pah11_con[index_con]

ewpah11_ohm_det = ewpah11_ohm[index_ohm] 
ewpah11_con_det = ewpah11_con[index_con]

ohm_ind = logoh_ohm[index_ohm]
con_ind = logoh_con[index_con]

con_ind_lim = logoh_con[1]

ndet = n_elements(ohm_ind) & ndet_con = n_elements(con_ind)
                             nlim_con = n_elements(con_ind_lim)

censor = [replicate('0',ndet),$
	replicate('-2', ndet_con),$
	replicate('-3', nlim_con)]

ind_data = string($
	ohm_ind, $
	con_ind, $
	con_ind_lim, $
	format='(f10.2)')

dep_data = string($
	pah11_ohm_det, $
	pah11_con_det, $
	pah11_con_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent variable
; Col 3 - dependent variable

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_pah11lum.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

dep_data = string($
	ewpah11_ohm_det, $
	ewpah11_con_det, $
	ewpah11_con_lim, $
	format='(f10.2)')

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_pah11ew.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; Water ice depth
;;;;;;;;;;;;;;;

ice_mega = transpose(-1d * float(ohmdat('ice')))
ice_arch = transpose(-1d * float(archdat('ice')))

ice_ohm = [ice_mega[*,0], ice_arch[*,0]]

ice_con = transpose(-1d * float(condat('ice')))
ice_con = ice_con[*,0]

restore, '~/Astronomy/Research/Spitzer/icelist62.sav'

match, tags, icelist62, oa, ob
match, tags_con, icelist62, ca, cb

ohm_det = ice_ohm[oa] 
con_det = ice_con[ca]

ohm_lim = ice_ohm[setdifference(indgen(n_elements(tags)),oa)]
con_lim = ice_con[setdifference(indgen(n_elements(tags_con)),ca)]

ohm_ind = logoh_ohm[oa]
con_ind = logoh_con[ca]

ohm_ind_lim = logoh_ohm[setdifference(indgen(n_elements(logoh_ohm)),oa)]
con_ind_lim = logoh_con[setdifference(indgen(n_elements(logoh_con)),ca)]

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate('0',ndet),$
	replicate('-1', nlim),$
	replicate('-2', ndet_con),$
	replicate('-3', nlim_con)]

ind_data = string($
	ohm_ind, $
	ohm_ind_lim, $
	con_ind, $
	con_ind_lim, $
	format='(f10.2)')

dep_data = string($
	ohm_det, $
	ohm_lim, $
	con_det, $
	con_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent variable
; Col 3 - dependent variable

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_ice.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; HAC absorption
;;;;;;;;;;;;;;;

; Entered in optical depths manually from table (did not save exact spline fits
; 	and features widths, it seems (unfortunately)

tau685_detection = [0.49, 0.20, 0.13, 0.35, 0.23, 0.37, 0.36, 0.15, 0.24, 0.23, 0.18, 0.25, $
	0.27, 0.16, 0.55, 0.41, 0.40, 0.24, 0.35, 0.45, 0.38, 0.23, 0.23, 0.21, 0.04, 0.19]
tau725_detection = [ 0.10, 0.22, 0.2 , 0.30, 0.15, 0.24, 0.21, 0.23]

restore, filename='~/Astronomy/Research/Spitzer/hac_limits.sav'

; 6.85 um feature

ohm_det = tau685_detection

ohm_lim = arr_685[1,where(strmid(arr_685[0,*],0,3) ne 'con')]
con_lim = arr_685[1,where(strmid(arr_685[0,*],0,3) eq 'con')]

match, tags,     transpose(arr_685[0,where(strmid(arr_685[0,*],0,3) ne 'con')]), oa, ob
match, tags_con, transpose(arr_685[0,where(strmid(arr_685[0,*],0,3) eq 'con')]), ca, cb

ohm_ind = logoh_ohm[setdifference(indgen(n_elements(tags)),oa)]

ohm_ind_lim = logoh_ohm[oa]
con_ind_lim = logoh_con[ca]

ndet = n_elements(ohm_det) 
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate('0',ndet),$
	replicate('-1', nlim),$
	replicate('-3', nlim_con)]

ind_data = string($
	ohm_ind, $
	ohm_ind_lim, $
	con_ind_lim, $
	format='(f10.2)')

dep_data = string($
	ohm_det, $
	ohm_lim, $
	con_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent variable
; Col 3 - dependent variable

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_hac685.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

; 7.25 um feature

ohm_det = tau725_detection

ohm_lim = arr_725[1,where(strmid(arr_725[0,*],0,3) ne 'con')]
con_lim = arr_725[1,where(strmid(arr_725[0,*],0,3) eq 'con')]

match, tags,     transpose(arr_725[0,where(strmid(arr_725[0,*],0,3) ne 'con')]), oa, ob
match, tags_con, transpose(arr_725[0,where(strmid(arr_725[0,*],0,3) eq 'con')]), ca, cb

ohm_ind = logoh_ohm[setdifference(indgen(n_elements(tags)),oa)]

ohm_ind_lim = logoh_ohm[oa]
con_ind_lim = logoh_con[ca]

ndet = n_elements(ohm_det) 
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate('0',ndet),$
	replicate('-1', nlim),$
	replicate('-3', nlim_con)]

ind_data = string($
	ohm_ind, $
	ohm_ind_lim, $
	con_ind_lim, $
	format='(f10.2)')

dep_data = string($
	ohm_det, $
	ohm_lim, $
	con_lim, $
	format='(f10.2)')

; Col 1 - censor
; Col 2 - independent variable
; Col 3 - dependent variable

array = [transpose(censor), transpose(ind_data), transpose(dep_data)]

filename = '~/Astronomy/Research/Spitzer/surv/bisurv_data/bisurv_hac725.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

end
