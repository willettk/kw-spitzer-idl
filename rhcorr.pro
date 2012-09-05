pro rhcorr, stop=stop
;+
; NAME:
;       
;	Look at correlations between IR properties and hyperfine ratios
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
;       Written by K. Willett                
;-


; List of properties to use

; Lo-res IRS data

nirs = 8

; NIR index

ospindex = ohmdat('spindex')
onir = float(transpose(ospindex(0,*)))

aspindex = archdat('spindex')
anir = float(transpose(aspindex(0,*)))


; MIR index

omir = float(transpose(ospindex(1,*)))
amir = float(transpose(aspindex(1,*)))

; PAH 6.2 um EW

otemp = ohmdat('pahfit62ew')
opah62ew = float(transpose(otemp(0,*)))

atemp = archdat('pahfit62ew')
apah62ew = float(transpose(atemp(0,*)))


; PAH 11.3 um EW

otemp = ohmdat('pahfit11ew')
opah11ew = float(transpose(otemp(0,*)))

atemp = archdat('pahfit11ew')
apah11ew = float(transpose(atemp(0,*)))


; PAH 6.2 um luminosity (PAHFIT)

otemp = ohmdat('pahfit62lum')
opah62lum = float(transpose(otemp(0,*)))

atemp = archdat('pahfit62lum')
apah62lum = float(transpose(atemp(0,*)))


; PAH 11.3 um luminosity (PAHFIT)

otemp = ohmdat('pahfit11lum')
opah11lum = float(transpose(otemp(0,*)))

atemp = archdat('pahfit11lum')
apah11lum = float(transpose(atemp(0,*)))


; Silicate strength

otemp = ohmdat('sil')
osil = float(transpose(otemp(0,*)))

atemp = archdat('sil')
asil = float(transpose(atemp(0,*)))


; Dust temperature

otemp = ohmdat('dtemp')
odtemp = float(transpose(otemp(0,*)))

atemp = archdat('dtemp')
adtemp = float(transpose(atemp(0,*)))


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

; f60/f100

ofluxratio = alog10(o60/o100)
afluxratio = alog10(a60/a100)

; Radio data

nrad = 4

o_dl = float(ohmdat('dl'))
a_dl = float(archdat('dl'))

; L_1667 - proportional to luminosity (needs factor of 4pi, units to work out)

o1667_flux = float(ohmdat('f1667'))
a1667_flux = float(archdat('f1667'))

o1667_lum = transpose(alog10(o1667_flux * o_dl^2))
a1667_lum = transpose(alog10(a1667_flux * a_dl^2))

; L_1420 

o1420_flux = float(ohmdat('f1420'))
a1420_flux = float(archdat('f1420'))

o1420_lum = transpose(alog10(o1420_flux * o_dl^2))
a1420_lum = transpose(alog10(a1420_flux * a_dl^2))

; log OH

ologoh = float(transpose(ohmdat('logoh')))
alogoh = float(transpose(archdat('logoh')))

; Hyperfine ratio

o_rh = float(transpose(ohmdat('rh')))
a_rh = float(transpose(archdat('rh')))

; Loop, eliminating bad data and making sure all object samples are identical

nvar = nrad + nirphot + nirs

ohm_matrix = fltarr(nvar,nvar)

om_string = strarr(nvar,nvar)+'######'

for j = 0, nvar-1 do begin

	case j of
	0: begin
		ohvar1 = [onir, anir]
		whatis1 = '15-6 slope'
	   end
	1: begin
		ohvar1 = [omir, amir]
		whatis1 = '30-20 slope'
	   end
	2: begin
		ohvar1 = [opah62ew,apah62ew]
		whatis1 = 'PAH 6.2 EW'
	   end
	3: begin
		ohvar1 = [opah11ew,apah11ew]
		whatis1 = 'PAH 11.3 EW'
	   end
	4: begin
		ohvar1 = [opah62lum,apah62lum]
		whatis1 = 'PAH 6.2 luminosity'
	   end
	5: begin
		ohvar1 = [opah11lum,apah11lum]
		whatis1 = 'PAH 11.3 luminosity'
	   end
	6: begin
		ohvar1 = [osil,asil]
		whatis1 = 'Silicate strength'
	   end
	7: begin
		ohvar1 = [odtemp,adtemp]
		whatis1 = 'Dust temperature'
	   end
	8: begin
		ohvar1 = [olfir,alfir]
		whatis1 = 'log L_FIR'
	   end
	9: begin
		ohvar1 = [ofluxratio,afluxratio]
		whatis1 = 'f60/f100'
	   end
	10: begin
		ohvar1 = [o1667_lum,a1667_lum]
		whatis1 = 'log L_1667'
	   end
	11: begin
		ohvar1 = [o1420_lum,a1420_lum]
		whatis1 = 'log L_1420'
	   end
	12: begin
		ohvar1 = [ologoh,alogoh]
		whatis1 = 'log L_OH'
	   end
	13: begin
		ohvar1 = [o_rh,a_rh]
		whatis1 = 'R_H'
	   end
	endcase

		ohvar2 = [o_rh,a_rh]
		whatis2 = 'R_H'

		; Find indices w/missing data

		badohm1 = where(ohvar1 eq 0 or finite(ohvar1) eq 0)
		badohm2 = where(ohvar2 eq 0 or finite(ohvar2) eq 0)

		; Check to make sure sizes of arrays are the same

		if n_elements(ohvar1) ne n_elements(ohvar2) then begin
			print,'Number of targets in correlating parameters for OHMs are not the same'
			print,whatis1
			print,whatis2
			stop
		endif

		; Find indices of good data

		goodindices_ohm = setdifference(indgen(n_elements(ohvar1)),[badohm1,badohm2])

		; Define variables for correlation

		ohvar1_new = ohvar1[goodindices_ohm]
		ohvar2_new = ohvar2[goodindices_ohm]

		; Option for using the Spearman's rank coefficient

			ohmcorr = r_correlate(ohvar1_new, ohvar2_new, zd=zd, probd=probd)

				print,'Variable 1: '+whatis1
				print,'Variable 2: '+whatis2
				print,'r_s: ',ohmcorr
				print,'zd: ',zd,' probd: ',probd
				print,''

			filetag='spearman'

endfor

if keyword_set(stop) then stop
end


