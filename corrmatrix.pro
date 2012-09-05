pro corrmatrix, stop=stop, spearman = spearman, verbose = verbose
;+
; NAME:
;       
;	CORRMATRIX
;
; PURPOSE:
;
;	Create a correlation matrix of various properties measured for OHM galaxies using the IRS
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
;	SPEARMAN - use Spearman's rank correlation coefficient (non-parametric) instead of Pearson's (default parametric)
;
;	VERBOSE - print correlation matrices
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
;       Written by K. Willett                Mar 08
; 	Added SPEARMAN keyword, converted f1420, f1667 to log luminosities - Mar 09
;-

; List of properties to use

; Lo-res IRS data

nirs = 8

; NIR index

ospindex = ohmdat('spindex')
onir = float(transpose(ospindex(0,*)))

aspindex = archdat('spindex')
anir = float(transpose(aspindex(0,*)))

cspindex = condat('spindex')
cnir = float(transpose(cspindex(0,*)))

; MIR index

omir = float(transpose(ospindex(1,*)))
amir = float(transpose(aspindex(1,*)))
cmir = float(transpose(cspindex(1,*)))

; PAH 6.2 um EW

otemp = ohmdat('pahfit62ew')
opah62ew = float(transpose(otemp(0,*)))

atemp = archdat('pahfit62ew')
apah62ew = float(transpose(atemp(0,*)))

ctemp = condat('pahfit62ew')
cpah62ew = float(transpose(ctemp(0,*)))

; PAH 11.3 um EW

otemp = ohmdat('pahfit11ew')
opah11ew = float(transpose(otemp(0,*)))

atemp = archdat('pahfit11ew')
apah11ew = float(transpose(atemp(0,*)))

ctemp = condat('pahfit11ew')
cpah11ew = float(transpose(ctemp(0,*)))

; PAH 6.2 um luminosity (PAHFIT)

otemp = ohmdat('pahfit62lum')
opah62lum = float(transpose(otemp(0,*)))

atemp = archdat('pahfit62lum')
apah62lum = float(transpose(atemp(0,*)))

ctemp = condat('pahfit62lum')
cpah62lum = float(transpose(ctemp(0,*)))

; PAH 11.3 um luminosity (PAHFIT)

otemp = ohmdat('pahfit11lum')
opah11lum = float(transpose(otemp(0,*)))

atemp = archdat('pahfit11lum')
apah11lum = float(transpose(atemp(0,*)))

ctemp = condat('pahfit11lum')
cpah11lum = float(transpose(ctemp(0,*)))

; Silicate strength

otemp = ohmdat('sil')
osil = float(transpose(otemp(0,*)))

atemp = archdat('sil')
asil = float(transpose(atemp(0,*)))

ctemp = condat('sil')
csil = float(transpose(ctemp(0,*)))

; Dust temperature

otemp = ohmdat('dtemp')
odtemp = float(transpose(otemp(0,*)))

atemp = archdat('dtemp')
adtemp = float(transpose(atemp(0,*)))

ctemp = condat('dtemp')
cdtemp = float(transpose(ctemp(0,*)))

; IR data - 2

nirphot = 2

; L_FIR

otemp = ohmdat('iras')
o60 = float(transpose(otemp(2,*)))
o100 = float(transpose(otemp(3,*)))
odl = float(transpose(ohmdat('dl')))
olfir = lfir(o60,o100,odl)

atemp = archdat('iras')
a60 = float(transpose(atemp(2,*)))
a100 = float(transpose(atemp(3,*)))
adl = float(transpose(archdat('dl')))
alfir = lfir(a60,a100,adl)

ctemp = condat('iras')
c60 = float(transpose(ctemp(2,*)))
c100 = float(transpose(ctemp(3,*)))
cdl = float(transpose(condat('dl')))
clfir = lfir(c60,c100,cdl)

; f60/f100

ofluxratio = alog10(o60/o100)
afluxratio = alog10(a60/a100)
cfluxratio = alog10(c60/c100)

; Radio data

nrad = 3

o_dl = float(ohmdat('dl'))
a_dl = float(archdat('dl'))
c_dl = float(condat('dl'))

; L_1667 - proportional to luminosity (needs factor of 4pi, units to work out)

o1667_flux = float(ohmdat('f1667'))
a1667_flux = float(archdat('f1667'))
c1667_flux = float(condat('f1667'))

o1667_lum = transpose(alog10(o1667_flux * o_dl^2))
a1667_lum = transpose(alog10(a1667_flux * a_dl^2))
c1667_lum = transpose(alog10(c1667_flux * c_dl^2))

; L_1420 

o1420_flux = float(ohmdat('f1420'))
a1420_flux = float(archdat('f1420'))
c1420_flux = float(condat('f1420'))

o1420_lum = transpose(alog10(o1420_flux * o_dl^2))
a1420_lum = transpose(alog10(a1420_flux * a_dl^2))
c1420_lum = transpose(alog10(c1420_flux * c_dl^2))

; log OH

ologoh = float(transpose(ohmdat('logoh')))
alogoh = float(transpose(archdat('logoh')))
clogoh = float(transpose(condat('logoh')))

; Loop, eliminating bad data and making sure all object samples are identical

nvar = nrad + nirphot + nirs

ohm_matrix = fltarr(nvar,nvar)
con_matrix = fltarr(nvar,nvar)

om_string = strarr(nvar,nvar)+'######'
cm_string = strarr(nvar,nvar)+'######'

for j = 0, nvar-1 do begin
	case j of
	0: begin
		ohvar1 = [onir, anir]
		convar1 = [cnir]
		whatis1 = '15-6 slope'
	   end
	1: begin
		ohvar1 = [omir, amir]
		convar1 = [cmir]
		whatis1 = '30-20 slope'
	   end
	2: begin
		ohvar1 = [opah62ew,apah62ew]
		convar1 = [cpah62ew]
		whatis1 = 'PAH 6.2 EW'
	   end
	3: begin
		ohvar1 = [opah11ew,apah11ew]
		convar1 = [cpah11ew]
		whatis1 = 'PAH 11.3 EW'
	   end
	4: begin
		ohvar1 = [opah62lum,apah62lum]
		convar1 = [cpah62lum]
		whatis1 = 'PAH 6.2 luminosity'
	   end
	5: begin
		ohvar1 = [opah11lum,apah11lum]
		convar1 = [cpah11lum]
		whatis1 = 'PAH 11.3 luminosity'
	   end
	6: begin
		ohvar1 = [osil,asil]
		convar1 = [csil]
		whatis1 = 'Silicate strength'
	   end
	7: begin
		ohvar1 = [odtemp,adtemp]
		convar1 = [cdtemp]
		whatis1 = 'Dust temperature'
	   end
	8: begin
		ohvar1 = [olfir,alfir]
		convar1 = [clfir]
		whatis1 = 'log L_FIR'
	   end
	9: begin
		ohvar1 = [ofluxratio,afluxratio]
		convar1 = [cfluxratio]
		whatis1 = 'f60/f100'
	   end
	10: begin
		ohvar1 = [o1667_lum,a1667_lum]
		convar1 = [c1667_lum]
		whatis1 = 'log L_1667'
	   end
	11: begin
		ohvar1 = [o1420_lum,a1420_lum]
		convar1 = [c1420_lum]
		whatis1 = 'log L_1420'
	   end
	12: begin
		ohvar1 = [ologoh,alogoh]
		convar1 = [clogoh]
		whatis1 = 'log L_OH'
	   end
	endcase

	for i = j, nvar-1 do begin

		case i of
		0: begin
			ohvar2 = [onir, anir]
			convar2 = [cnir]
			whatis2 = 'NIR'
		   end
		1: begin
			ohvar2 = [omir, amir]
			convar2 = [cmir]
			whatis2 = 'MIR'
		   end
		2: begin
			ohvar2 = [opah62ew,apah62ew]
			convar2 = [cpah62ew]
			whatis2 = 'PAH 6.2 EW'
		   end
		3: begin
			ohvar2 = [opah11ew,apah11ew]
			convar2 = [cpah11ew]
			whatis2 = 'PAH 11.3 EW'
		   end
		4: begin
			ohvar2 = [opah62lum,apah62lum]
			convar2 = [cpah62lum]
			whatis2 = 'PAH 6.2 luminosity'
		   end
		5: begin
			ohvar2 = [opah11lum,apah11lum]
			convar2 = [cpah11lum]
			whatis2 = 'PAH 11.3 luminosity'
		   end
		6: begin
			ohvar2 = [osil,asil]
			convar2 = [csil]
			whatis2 = 'Silicate strength'
		   end
		7: begin
			ohvar2 = [odtemp,adtemp]
			convar2 = [cdtemp]
			whatis2 = 'Dust temperature'
		   end
		8: begin
			ohvar2 = [olfir,alfir]
			convar2 = [clfir]
			whatis2 = 'L_FIR'
		   end
		9: begin
			ohvar2 = [ofluxratio,afluxratio]
			convar2 = [cfluxratio]
			whatis2 = 'f60/f100'
		   end
		10: begin
			ohvar2 = [o1667_lum,a1667_lum]
			convar2 = [c1667_lum]
			whatis2 = 'log L_1667'
		   end
		11: begin
			ohvar2 = [o1420_lum,a1420_lum]
			convar2 = [c1420_lum]
			whatis2 = 'log L_1420'
		   end
		12: begin
			ohvar2 = [ologoh,alogoh]
			convar2 = [clogoh]
			whatis2 = 'log L_OH'
		   end
		endcase

		; Find indices w/missing data

		badohm1 = where(ohvar1 eq 0 or finite(ohvar1) eq 0)
		badohm2 = where(ohvar2 eq 0 or finite(ohvar2) eq 0)
		badcon1 = where(convar1 eq 0 or finite(convar1) eq 0)
		badcon2 = where(convar2 eq 0 or finite(convar2) eq 0)

		; Check to make sure sizes of arrats are the same

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

		; Option for using the Spearman's rank coefficient

		if keyword_set(spearman) then begin

			ohmcorr = r_correlate(ohvar1_new, ohvar2_new, zd=zd,     probd=probd)
			concorr = r_correlate(convar1_new,convar2_new,zd=zd_con, probd=probd_con)

			ohm_matrix[i,j] = ohmcorr[0]
			con_matrix[i,j] = concorr[0]
	
			om_string[i,j] = string(ohmcorr[0],format='(f6.2)')
			cm_string[i,j] = string(concorr[0],format='(f6.2)')

;			if abs(zd) gt 3.0  and i ne j then begin
;				print,'Variable 1: '+whatis1
;				print,'Variable 2: '+whatis2
;				print,'r_s: ',ohmcorr
;				print,'zd: ',zd,' probd: ',probd
;				print,''
;			endif

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

		; Option for using the Pearson linear correlation coefficient

		endif else begin

			ohmcorr = correlate(ohvar1_new,ohvar2_new)
			concorr = correlate(convar1_new,convar2_new)

			ohm_matrix[i,j] = ohmcorr
			con_matrix[i,j] = concorr
	
			om_string[i,j] = string(ohmcorr,format='(f6.2)')
			cm_string[i,j] = string(concorr,format='(f6.2)')
	
			N = n_elements(ohvar1_new)
			r = ohmcorr
			t = r * sqrt(N - 2) / sqrt(1d - r^2)		; Wall & Jenkins, p. 62 - Student's t-test

;			if i ne j and abs(ohmcorr) gt 0.7 then begin
			if i ne j and t gt 2.7 then begin		; 2-sigma significance for nu = 50 is roughly 2.7 
				print,'Variable 1: '+whatis1
				print,'Variable 2: '+whatis2
				print,'r = ',ohmcorr
				print,'t = ',t
				print,''
			endif

			filetag='pearson'

		endelse
	endfor
endfor

; Write matrices to file

writeohm = '~/Astronomy/Research/Spitzer/ohm/ohm_corrmatrix_'+filetag+'.txt'
writecon = '~/Astronomy/Research/Spitzer/control/con_corrmatrix_'+filetag+'.txt'

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
