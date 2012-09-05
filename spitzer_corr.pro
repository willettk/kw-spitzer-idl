pro spitzer_corr, stop=stop, verbose=verbose
;+
; NAME:
;       
;	SPITZER_CORR
;
; PURPOSE:
;
;	Return correlation results for two parameters
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
;	VERBOSE - print correlation matrices
;
; EXAMPLE:
;
;	IDL> .r corrmatrix4
;
; NOTES:
;
;	A shortened version of CORRMATRIX3, removing all parameters that did not reveal any significant rank correlations. 
;
; REVISION HISTORY
;       Written by K. Willett                Mar 2010
;-

otag = ohmdat('tag')
atag = archdat('tag')
ctag = condat('tag')

; List of properties to use

; Lo-res IRS data

nir = 11

; L_FIR

otemp = ohmdat('lfir')
olfir = float(transpose(otemp))

atemp = archdat('lfir')
alfir = float(transpose(atemp))

ctemp = condat('lfir')
clfir = float(transpose(ctemp))

; NIR index

ospindex = ohmdat('spindex')
aspindex = archdat('spindex')
cspindex = condat('spindex')

onir = float(transpose(ospindex[0,*]))
anir = float(transpose(aspindex[0,*]))
cnir = float(transpose(cspindex[0,*]))

; MIR index

omir = float(transpose(ospindex[1,*]))
amir = float(transpose(aspindex[1,*]))
cmir = float(transpose(cspindex[1,*]))

; PAH 6.2 um EW

otemp = ohmdat('pahfit62ew')
opah62ew = float(transpose(otemp[0,*]))

atemp = archdat('pahfit62ew')
apah62ew = float(transpose(atemp[0,*]))

ctemp = condat('pahfit62ew')
cpah62ew = float(transpose(ctemp[0,*]))

; PAH 6.2 um luminosity (PAHFIT)

otemp = ohmdat('pahfit62lum')
opah62lum = float(transpose(otemp[0,*]))

atemp = archdat('pahfit62lum')
apah62lum = float(transpose(atemp[0,*]))

ctemp = condat('pahfit62lum')
cpah62lum = float(transpose(ctemp[0,*]))

; NeV 14 / NeII ratio

one5ne2 = fltarr(n_elements(otag))

for i=0, n_elements(otag) - 1 do begin
	otemp2 = getlineflux(/quiet,otag[i], 'neII')
	otemp5 = getlineflux(/quiet,otag[i], 'neV')
	if otemp2[0] gt 0. and otemp5[0] gt 0. then one5ne2[i] = float(otemp5[0])/float(otemp2[0])
endfor

ane5ne2 = fltarr(n_elements(atag))

for i=0, n_elements(atag) - 1 do begin
	atemp2 = getlineflux(/quiet,atag[i], 'neII')
	atemp5 = getlineflux(/quiet,atag[i], 'neV')
	if atemp2[0] gt 0. and atemp5[0] gt 0. then ane5ne2[i] = float(atemp5[0])/float(atemp2[0])
endfor

cne5ne2 = fltarr(n_elements(ctag))

for i=0, n_elements(ctag) - 1 do begin
	ctemp2 = getlineflux(/quiet,ctag[i], 'neII')
	ctemp5 = getlineflux(/quiet,ctag[i], 'neV')
	if ctemp2[0] gt 0. and ctemp5[0] gt 0. then cne5ne2[i] = float(ctemp5[0])/float(ctemp2[0])
endfor

; 9.7 um silicate depth

osiltemp = ohmdat('sil')
asiltemp = archdat('sil')
csiltemp = condat('sil')

osil10 = transpose(float(osiltemp[0,*]))
asil10 = transpose(float(asiltemp[0,*]))
csil10 = transpose(float(csiltemp[0,*]))

; 18 um silicate depth

osil18 = transpose(float(osiltemp[1,*]))
asil18 = transpose(float(asiltemp[1,*]))
csil18 = transpose(float(csiltemp[1,*]))

; Blackbody dust temperature

otemp = ohmdat('dtemp')
odtemp_bb = float(transpose(otemp(0,*)))

atemp = archdat('dtemp')
adtemp_bb = float(transpose(atemp(0,*)))

ctemp = condat('dtemp')
cdtemp_bb = float(transpose(ctemp(0,*)))

; DUSTY best fits (Y, tau_V, T_dust)

oy = fltarr(n_elements(otag))
otauv = fltarr(n_elements(otag))
odtemp_dusty = fltarr(n_elements(otag))

for i=0, n_elements(otag) - 1 do begin
	restore,'~/Astronomy/Research/Spitzer/ohm/DUSTY/sphere_sav/'+otag[i]+'_sphere_pahfitgrid.sav'
	oy[i] = y_min
	otauv[i] = tauv_min
	odtemp_dusty[i] = tdust
endfor

ay = fltarr(n_elements(atag))
atauv = fltarr(n_elements(atag))
adtemp_dusty = fltarr(n_elements(atag))

for i=0, n_elements(atag) - 1 do begin
	restore,'~/Astronomy/Research/Spitzer/archived/DUSTY/sphere_sav/'+atag[i]+'_sphere_pahfitgrid.sav'
	ay[i] = y_min
	atauv[i] = tauv_min
	adtemp_dusty[i] = tdust
endfor

cy = fltarr(n_elements(ctag))
ctauv = fltarr(n_elements(ctag))
cdtemp_dusty = fltarr(n_elements(ctag))

for i=0, n_elements(ctag) - 1 do begin
	restore,'~/Astronomy/Research/Spitzer/control/DUSTY/sphere_sav/'+ctag[i]+'_sphere_pahfitgrid.sav'
	cy[i] = y_min
	ctauv[i] = tauv_min
	cdtemp_dusty[i] = tdust
endfor

; Radio data

nrad = 3

o_dl = float(ohmdat('dl'))
a_dl = float(archdat('dl'))
c_dl = float(condat('dl'))

; P_1667 - proportional to power at peak of OHM

o1667_flux = float(ohmdat('f1667'))
a1667_flux = float(archdat('f1667'))
c1667_flux = float(condat('f1667'))

o1667_power = transpose(alog10(o1667_flux * o_dl^2))
a1667_power = transpose(alog10(a1667_flux * a_dl^2))
c1667_power = transpose(alog10(c1667_flux * c_dl^2))

; P_1420 - proportional to power at 1.4 GHz

o1420_flux = float(ohmdat('f1420'))
a1420_flux = float(archdat('f1420'))
c1420_flux = float(condat('f1420'))

o1420_power = transpose(alog10(o1420_flux * o_dl^2))
a1420_power = transpose(alog10(a1420_flux * a_dl^2))
c1420_power = transpose(alog10(c1420_flux * c_dl^2))

; log OH

ologoh = float(transpose(ohmdat('logoh')))
alogoh = float(transpose(archdat('logoh')))
clogoh = float(transpose(condat('logoh')))

; Loop, eliminating bad data and making sure all object samples are identical

nvar = nrad + nir

ohm_matrix = fltarr(nvar,nvar)
con_matrix = fltarr(nvar,nvar)

om_string = strarr(nvar,nvar)+' \ldots '
cm_string = strarr(nvar,nvar)+' \ldots '

for j = 0, nvar-1 do begin
	case j of
	0: begin
		ohvar1 = [olfir,alfir]
		convar1 = [clfir]
		whatis1 = 'log L_FIR'
	   end
	1: begin
		ohvar1 = [onir, anir]
		convar1 = [cnir]
		whatis1 = '15-6 slope'
	   end
	2: begin
		ohvar1 = [omir, amir]
		convar1 = [cmir]
		whatis1 = '30-20 slope'
	   end
	3: begin
		ohvar1 = [opah62ew,apah62ew]
		convar1 = [cpah62ew]
		whatis1 = 'PAH 6.2 EW'
	   end
	4: begin
		ohvar1 = [opah62lum,apah62lum]
		convar1 = [cpah62lum]
		whatis1 = 'PAH 6.2 luminosity'
	   end
	5: begin
		ohvar1 = [osil10,asil10]
		convar1 = [csil10]
		whatis1 = 'S_9.7'
	   end
	6: begin
		ohvar1 = [osil18,asil18]
		convar1 = [csil18]
		whatis1 = 'S_18'
	   end
	7: begin
		ohvar1 = [odtemp_bb,adtemp_bb]
		convar1 = [cdtemp_bb]
		whatis1 = 'Blackbody dust temp.'
	   end
	8: begin
		ohvar1 = [oy,ay]
		convar1 = [cy]
		whatis1 = 'Y from DUSTY'
	   end
	9: begin
		ohvar1 = [otauv,atauv]
		convar1 = [ctauv]
		whatis1 = 'tau_V from DUSTY'
	   end
	10: begin
		ohvar1 = [odtemp_dusty,adtemp_dusty]
		convar1 = [cdtemp_dusty]
		whatis1 = 'Dust temp. from DUSTY'
	   end
	11: begin
		ohvar1 = [o1667_power,a1667_power]
		convar1 = [c1667_power]
		whatis1 = 'log P_1667'
	   end
	12: begin
		ohvar1 = [o1420_power,a1420_power]
		convar1 = [c1420_power]
		whatis1 = 'log P_1420'
	   end
	13: begin
		ohvar1 = [ologoh,alogoh]
		convar1 = [clogoh]
		whatis1 = 'log L_OH'
	   end
	endcase

	for i = j, nvar-1 do begin

		case i of
			0: begin
				ohvar2 = [olfir,alfir]
				convar2 = [clfir]
				whatis2 = 'log L_FIR'
			   end
			1: begin
				ohvar2 = [onir, anir]
				convar2 = [cnir]
				whatis2 = '15-6 slope'
			   end
			2: begin
				ohvar2 = [omir, amir]
				convar2 = [cmir]
				whatis2 = '30-20 slope'
			   end
			3: begin
				ohvar2 = [opah62ew,apah62ew]
				convar2 = [cpah62ew]
				whatis2 = 'PAH 6.2 EW'
			   end
			4: begin
				ohvar2 = [opah62lum,apah62lum]
				convar2 = [cpah62lum]
				whatis2 = 'PAH 6.2 luminosity'
			   end
			5: begin
				ohvar2 = [osil10,asil10]
				convar2 = [csil10]
				whatis2 = 'S_9.7'
			   end
			6: begin
				ohvar2 = [osil18,asil18]
				convar2 = [csil18]
				whatis2 = 'S_18'
			   end
			7: begin
				ohvar2 = [odtemp_bb,adtemp_bb]
				convar2 = [cdtemp_bb]
				whatis2 = 'Blackbody dust temp.'
			   end
			8: begin
				ohvar2 = [oy,ay]
				convar2 = [cy]
				whatis2 = 'Y from DUSTY'
			   end
			9: begin
				ohvar2 = [otauv,atauv]
				convar2 = [ctauv]
				whatis2 = 'tau_V from DUSTY'
			   end
			10: begin
				ohvar2 = [odtemp_dusty,adtemp_dusty]
				convar2 = [cdtemp_dusty]
				whatis2 = 'Dust temp. from DUSTY'
			   end
			11: begin
				ohvar2 = [o1667_power,a1667_power]
				convar2 = [c1667_power]
				whatis2 = 'log P_1667'
			   end
			12: begin
				ohvar2 = [o1420_power,a1420_power]
				convar2 = [c1420_power]
				whatis2 = 'log P_1420'
			   end
			13: begin
				ohvar2 = [ologoh,alogoh]
				convar2 = [clogoh]
				whatis2 = 'log L_OH'
			   end
		endcase

		; Find indices w/missing data

			badohm1 = where(ohvar1 eq 0 or finite(ohvar1) eq 0)
			badcon1 = where(convar1 eq 0 or finite(convar1) eq 0)
			badohm2 = where(ohvar2 eq 0 or finite(ohvar2) eq 0)
			badcon2 = where(convar2 eq 0 or finite(convar2) eq 0)

		; Check to make sure sizes of arrays are the same

		if n_elements(ohvar1) ne n_elements(ohvar2) then begin
			print,'Number of targets in correlating parameters for OHMs are not the same'
			print,whatis1
			print,whatis2
			stop
		endif

		if n_elements(convar1) ne n_elements(convar2) then begin
			print,'Number of targets in correlating parameters for control sample are not the same'
			print,whatis1
			print,whatis2
			stop
		endif

		; Find indices of good data

		goodindices_ohm = setdifference(indgen(n_elements(ohvar1)),[badohm1,badohm2])
		goodindices_con = setdifference(indgen(n_elements(convar1)),[badcon1,badcon2])

		; Define variables for correlation

		ohvar1_new = ohvar1(goodindices_ohm)
		ohvar2_new = ohvar2(goodindices_ohm)

		convar1_new = convar1(goodindices_con)
		convar2_new = convar2(goodindices_con)

		; Using the Spearman's rank coefficient

			ohmcorr = r_correlate(ohvar1_new, ohvar2_new, zd=zd,     probd=probd)
			concorr = r_correlate(convar1_new,convar2_new,zd=zd_con, probd=probd_con)

			ohm_matrix[i,j] = zd[0]
			con_matrix[i,j] = zd_con[0]
	
			om_string[i,j] = string(-1 * zd[0],format='(f8.1)')
			cm_string[i,j] = string(-1 * zd_con[0],format='(f8.1)')

			if i eq j then begin
				om_string[i,j] = '   --    '
				cm_string[i,j] = '   --    '
			endif

			; Print strong correlations from both OHMs and non-masing

			if abs(zd) gt 4.0  and i ne j then begin
				print,'Variable 1: '+whatis1
				print,'Variable 2: '+whatis2
				print,'r_s: ',ohmcorr
				print,'zd: ',zd
				print,'probd: ',probd
				print,''
			endif

			; Print differences between correlations in OHMs and control samples

			if abs(zd - zd_con) gt 4.0 and i ne j then begin
				print,'Variable 1: '+whatis1
				print,'Variable 2: '+whatis2
				print,'r_OHM: ',ohmcorr
				print,'r_con: ',concorr
				print,'zd_OHM: ',zd
				print,'zd_con: ',zd_con
				print,''
			endif

			filetag='spearman'

	endfor
endfor

; Write matrices to file

writeohm = '~/Astronomy/Research/Spitzer/ohm/ohm_corrmatrix4_'+filetag+'.txt'
writecon = '~/Astronomy/Research/Spitzer/control/con_corrmatrix4_'+filetag+'.txt'

openw,lun1,writeohm,/get_lun
printf,lun1,string(ohm_matrix,format='(f8.2)')
close,lun1

openw,lun2,writecon,/get_lun
printf,lun2,string(con_matrix,format='(f8.2)')
close,lun2

close,/all

; Option to print correlation matrices to screen

if keyword_set(verbose) then begin
	print,''
	print,'OHMs'
	print,om_string
	print,''
	print,'Control'
	print,cm_string
	print,''
endif

if keyword_set(stop) then stop
end
