
;+
; NAME:
;       
;	SURV_BATCH
;
; PURPOSE:
;
;	Create files readable by IRAF's survival analysis package (ASURV) for various mid-IR properties
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
; REVISION HISTORY
;       Written by K. Willett
;	Added crystalline silicates - Feb 10
;-

tags = [transpose(ohmdat('tag')),transpose(archdat('tag'))]
objs = [transpose(ohmdat('obj')),transpose(archdat('obj'))]
dl   = [transpose(ohmdat('dl')),transpose(archdat('dl'))]

tags_con = [transpose(condat('tag'))]
objs_con = [transpose(condat('obj'))]
dl_con = [transpose(condat('dl'))]

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

ohm_det = f60_ohm[where(f100_ohm ne 0)] / f100_ohm[where(f100_ohm ne 0)]
con_det = f60_con[where(f100_con ne 0)] / f100_con[where(f100_con ne 0)]

ohm_lim = f60_ohm[where(f100_ohm eq 0)] / lim100_ohm
con_lim = f60_con[where(f100_con eq 0)] / lim100_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate('0',ndet),replicate('1',nlim),$			; Note that this is a LOWER limit
	replicate('0',ndet_con),replicate('1',nlim_con)]

data = string($
	ohm_det, $
	ohm_lim, $
	con_det, $
	con_lim, $
	format='(f10.2)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

; Col 1 - censor
; Col 2 - data
; Col 3 - group

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_f60_f100.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; Flux and luminosity of 1420 MHz continuum
;;;;;;;;;;;;;;;

f1420_mega = transpose(float(ohmdat('f1420')))
f1420_arch = transpose(float(archdat('f1420')))

f1420_ohm = [f1420_mega, f1420_arch]

f1420_con = transpose(float(condat('f1420')))

; Upper limits on data from NED

;lim1420_ohm = [28,48,54]		<--- Indices of data w/limits; 5 mJy for all targets (limits of NVSS)
;lim1420_con = [0,3,5,14]

ohm_det = f1420_ohm[where(f1420_ohm ne 0)]
con_det = f1420_con[where(f1420_con ne 0)]

ohm_lim = replicate(5., 3)
con_lim = replicate(5., 4)

dl_det     = dl[where(f1420_ohm ne 0)]
dl_lim     = dl[where(f1420_ohm eq 0)]
dl_det_con = dl[where(f1420_con ne 0)]
dl_lim_con = dl[where(f1420_con eq 0)]

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet),replicate('-1',nlim),$			; UPPER limit
	replicate(' 0',ndet_con),replicate('-1',nlim_con)]

data = string($
	ohm_det, $
	ohm_lim, $
	con_det, $
	con_lim, $
	format='(f10.2)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_f1420.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

data = string(alog10([$
	ohm_det * dl_det^2, $
	ohm_lim * dl_lim^2, $
	con_det * dl_det_con^2, $
	con_lim * dl_lim_con^2]), $
	format='(f10.2)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_lum1420.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; Luminosity of PAH features
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

ind_ohm = setdifference(indgen(n_elements(pah62_ohm)),[5,40,25])
ind_con = setdifference(indgen(n_elements(pah62_con)),[1,8])

ohm_det = pah62_ohm[ind_ohm] & con_det = pah62_con[ind_con]
ewohm_det = ewpah62_ohm[ind_ohm] & ewcon_det = ewpah62_con[ind_con]

; Upper limits on data from NED

;lim62_ohm = [5,40]		;      <--- Indices of data w/limits
;lim62_con = [1]
;lim11_ohm = []
;lim11_con = [1]

ohm_lim = [8.95, 9.79]
con_lim = [9.39]

ew_ohm_lim = [0.10, 0.01]
ew_con_lim = [0.01]

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet),replicate('-1',nlim),$				; UPPER limit
	replicate(' 0',ndet_con),replicate('-1',nlim_con)]

data = string($
	ohm_det, $
	ohm_lim, $
	con_det, $
	con_lim, $
	format='(f10.2)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_pah62lum.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

data = string($
	ewohm_det, $
	ew_ohm_lim, $
	ewcon_det, $
	ew_con_lim, $
	format='(f10.2)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_pah62ew.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

ind_ohm = setdifference(indgen(n_elements(pah11_ohm)),[25])
ind_con = setdifference(indgen(n_elements(pah11_con)),[1,8])

ohm_det = pah11_ohm[ind_ohm] & con_det = pah11_con[ind_con]
ewohm_det = ewpah11_ohm[ind_ohm] & ewcon_det = ewpah11_con[ind_con]

con_lim = [9.10]
ew_con_lim = [0.01]

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet), 			$			; UPPER limit
	replicate(' 0',ndet_con),replicate('-1',nlim_con)]

data = string([$
	ewohm_det, $
	ewcon_det, $
	ew_con_lim], $
	format='(f10.2)')

group = [replicate('1', ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_pah11lum.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

data = string([$
	ohm_det, $
	con_det, $
	con_lim], $
	format='(f10.2)')

group = [replicate('1', ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_pah11ew.txt'

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

ndet = n_elements(ohm_det) 
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet),replicate('-1',nlim),$				; UPPER limit
	replicate('-1',nlim_con)]

data = string($
	ohm_det, $
	ohm_lim, $
	con_lim, $
	format='(f10.3)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_hac685.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

; 7.25 um feature

ohm_det = tau725_detection

ohm_lim = arr_725[1,where(strmid(arr_725[0,*],0,3) ne 'con')]
con_lim = arr_725[1,where(strmid(arr_725[0,*],0,3) eq 'con')]

ndet = n_elements(ohm_det) 
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet),replicate('-1',nlim),$				; UPPER limit
	replicate('-1',nlim_con)]

data = string($
	ohm_det, $
	ohm_lim, $
	con_lim, $
	format='(f10.3)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_hac725.txt'

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

ohm_det = ice_ohm[oa] & con_det = ice_con[ca]

ohm_lim = ice_ohm[setdifference(indgen(n_elements(tags)),oa)]
con_lim = ice_con[setdifference(indgen(n_elements(tags_con)),ca)]

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)
nlim = n_elements(ohm_lim) & nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet),replicate('-1',nlim),$				; UPPER limit
	replicate(' 0',ndet_con),replicate('-1',nlim_con)]

data = string($
	ohm_det, $
	ohm_lim, $
	con_det, $
	con_lim, $
	format='(f10.3)')

group = [replicate('1', nlim + ndet), replicate('2', nlim_con + ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_ice.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun


; Perform same tests on objects with no limits, for completeness' sake


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

ohm_det = slope15_ohm & con_det = slope15_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f10.3)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_slope15.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

ohm_det = slope30_ohm & con_det = slope30_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f10.3)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_slope30.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun


;;;;;;;;;;;;;;;
; H2 mass and temperature
;;;;;;;;;;;;;;;

restore,'~/Astronomy/Research/Spitzer/ohm/h2.sav'

h2_temp_warm_ohm = h2.temp_warm
h2_temp_hot_ohm = h2.temp_hot
h2_mass_warm_ohm = h2.mass_warm
h2_mass_hot_ohm = h2.mass_hot

h2_temp_warm_ohm = h2_temp_warm_ohm[where(h2_temp_warm_ohm ne 0.)]
h2_temp_hot_ohm = h2_temp_hot_ohm[where(h2_temp_hot_ohm ne 0.)]
h2_mass_warm_ohm = h2_mass_warm_ohm[where(h2_mass_warm_ohm ne 0.)]
h2_mass_hot_ohm = h2_mass_hot_ohm[where(h2_mass_hot_ohm ne 0.)]

restore,'~/Astronomy/Research/Spitzer/control/h2_con.sav'

h2_temp_warm_con = h2.temp_warm
h2_temp_hot_con = h2.temp_hot
h2_mass_warm_con = h2.mass_warm
h2_mass_hot_con = h2.mass_hot

h2_temp_warm_con = h2_temp_warm_con[where(h2_temp_warm_con ne 0.)]
h2_temp_hot_con = h2_temp_hot_con[where(h2_temp_hot_con ne 0.)]
h2_mass_warm_con = h2_mass_warm_con[where(h2_mass_warm_con ne 0.)]
h2_mass_hot_con = h2_mass_hot_con[where(h2_mass_hot_con ne 0.)]

ohm_det = h2_temp_warm_ohm & con_det = h2_temp_warm_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				; No limits
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(i4)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_h2_temp_warm.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

ohm_det = h2_temp_hot_ohm & con_det = h2_temp_hot_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				; No limits
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(i4)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_h2_temp_hot.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

ohm_det = h2_mass_warm_ohm & con_det = h2_mass_warm_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				; No limits
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f10.4)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_h2_mass_warm.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

ohm_det = h2_mass_hot_ohm & con_det = h2_mass_hot_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				; No limits
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f10.4)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_h2_mass_hot.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun


;;;;;;;;;;;;;;;
; Dust temperature
;;;;;;;;;;;;;;;

dtemp_mega = transpose(float(ohmdat('dtemp')))
dtemp_arch = transpose(float(archdat('dtemp')))

dtemp_con = transpose(float(condat('dtemp')))

dtemp_ohm = [dtemp_mega[*,0], dtemp_arch[*,0]]
dtemp_con = dtemp_con[*,0]

ohm_det = dtemp_ohm & con_det = dtemp_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				; UPPER limit
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f10.1)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_dtemp.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun


;;;;;;;;;;;;;;;
; Silicate depth
;;;;;;;;;;;;;;;

silicate_mega = transpose(-1d * float(ohmdat('sil')))
silicate_arch = transpose(-1d * float(archdat('sil')))

silicate_con = transpose(-1d * float(condat('sil')))

silicate_ohm = [silicate_mega[*,0], silicate_arch[*,0]]
silicate_con = silicate_con[*,0]

silicate_ohm = silicate_ohm[where(silicate_ohm ne 0.)]
silicate_con = silicate_con[where(silicate_con ne 0.)]

ohm_det = silicate_ohm & con_det = silicate_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f6.2)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_silicate.txt'

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

ohm_det = lfir_ohm & con_det = lfir_con

ndet = n_elements(ohm_det) & ndet_con = n_elements(con_det)

censor = [replicate(' 0',ndet),$				
	replicate(' 0',ndet_con)]

data = string($
	ohm_det, $
	con_det, $
	format='(f6.2)')

group = [replicate('1', ndet), replicate('2', ndet_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_lfir.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

;;;;;;;;;;;;;;;
; Crystalline silicates
;;;;;;;;;;;;;;;

csil_tags = ['mega'+string([8,9,25,26,27,32],format='(i03)'),$
	'arch'+string([4,5,9,10,12,18,29,30,31,32,33,36,45], format='(i03)')]

csil_tags_con = ['control'+string([35],format='(i03)')]

alltags = [transpose(ohmdat('tag')),transpose(archdat('tag'))]
contags = condat('tag')

match, alltags, csil_tags, a, b
nocrysind = setdifference(indgen(n_elements(alltags)),a)

; Empty arrays 

csil16 = fltarr(n_elements(csil_tags))
csil23 = fltarr(n_elements(csil_tags))
csil16_con = fltarr(n_elements(csil_tags_con))
csil23_con = fltarr(n_elements(csil_tags_con))
csil16lim = fltarr(n_elements(nocrysind))
csil23lim = fltarr(n_elements(nocrysind))
csil16lim_con = fltarr(n_elements(contags))
csil23lim_con = fltarr(n_elements(contags))

match, contags, csil_tags_con, c, d
nocrysind_con = setdifference(indgen(n_elements(contags)),c)

; Run detections

for i = 0, n_elements(csil_tags) - 1 do begin
	crystalline, csil_tags[i], tau16, tau16err, tau23, tau23err,/quiet, /noplot
	csil16[i] = -1 * tau16
	csil23[i] = -1 * tau23
endfor

for i = 0, n_elements(csil_tags_con) - 1 do begin
	crystalline, csil_tags_con[i], tau16, tau16err, tau23, tau23err,/quiet, /noplot
	csil16_con[i] = -1 * tau16
	csil23_con[i] = -1 * tau23
endfor

; Run non-detections

for i = 0, n_elements(nocrysind) - 1 do begin
	crystalline, alltags[nocrysind[i]], tau16, tau16err, tau23, tau23err,/quiet, /noplot
	csil16lim[i] = -1 * tau16
	csil23lim[i] = -1 * tau23
endfor

for i = 0, n_elements(nocrysind_con) - 1 do begin
	crystalline, contags[nocrysind_con[i]], tau16, tau16err, tau23, tau23err, /quiet, /noplot
	csil16lim_con[i] = -1 * tau16
	csil23lim_con[i] = -1 * tau23
endfor

ohm_det = csil16
con_det = csil16_con
ohm_lim = csil16lim
con_lim = csil16lim_con

ndet     = n_elements(ohm_det)
ndet_con = n_elements(con_det)
nlim     = n_elements(ohm_lim)
nlim_con = n_elements(con_lim)

censor = [replicate(' 0',ndet), $
	replicate(' 0',ndet_con), $
	replicate('-1',nlim), $				; UPPER limit
	replicate('-1',nlim_con)]

data = string($
	ohm_det, $
	con_det, $
	ohm_lim, $
	con_lim, $
	format='(f6.2)')

group = [replicate('1',ndet), $
	replicate('2',ndet_con), $
	replicate('1',nlim), $
	replicate('2',nlim_con)]

array = [transpose(censor), transpose(data), transpose(group)]

filename = '~/Astronomy/Research/Spitzer/surv/univ_data/univ_crystalline.txt'

openw, lun, filename, /get_lun
printf, lun, '# Censor	Data	Group'
printf, lun, array
free_lun, lun

; Create columns for table

;a = [csil_tags,alltags[nocrysind]]
;na = n_elements(a)
;b = strarr(na)
;for i = 0, na - 1 do begin
;	targets,a[i],r,o,z
;	b[i] = o
;endfor
;
;c=string([ohm_det,ohm_lim],format='(f4.2)')
;c[n_elements(ohm_det):n_elements(c)-1] = '--  '
;print,[transpose(b[sort(b)]),transpose(c[sort(b)])]
;
;print,''
;
;aa = [csil_tags_con,contags[nocrysind_con]]
;naa = n_elements(aa)
;bb = strarr(naa)
;for i = 0, naa - 1 do begin
;	targets,aa[i],r,o,z
;	bb[i] = o
;endfor
;
;cc=string([con_det,con_lim],format='(f4.2)')
;cc[n_elements(con_det):n_elements(cc)-1] = '--  '
;print,[transpose(bb[sort(bb)]),transpose(cc[sort(bb)])]

end
