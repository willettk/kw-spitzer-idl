
;+
; NAME:
;       
;	NED_IRAS
;
; PURPOSE:
;
;	Extract IRAS fluxes from a formatted NED batch file
;
; INPUTS:
;
;	AFILE_IR -	text file containing names and photometry from the NED batch system
;
; OUTPUTS:
;
;	IDL sav file with fluxes, errors, and names
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> .r ned_iras
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;-

; Files on which to operate

spitzdir = '~/Astronomy/Research/Spitzer/'
filelist = spitzdir+['ohm/ned_mega.lst','archived/ned_arch.lst','cso/ned_cso.lst','control/ned_con.lst']
writefile = spitzdir+['ohm/data/idl_sav/ohm_iras.sav', $
	'archived/data/idl_sav/arch_iras.sav', $
	'cso/data/idl_sav/cso_iras.sav', $
	'control/data/idl_sav/con_iras.sav']

for j = 0, n_elements(filelist) do begin

	openr,lun,filelist(j),/get_lun
	nobj = numlines(filelist(j))  & nasource = nobj/5
	if nobj mod 5 ne 0 then print,'Wrong number of files in list.'
	irflux = strarr(nobj)
	readf,lun,irflux
	close,lun
	
	; Extract the object names
	
	aname_ir = irflux(indgen(nasource)*5)
	aname_ir = strmid(aname_ir,10,100)
	
	; Format the names
	
	for i = 0, nasource - 1 do begin
		temp2 = strsplit(aname_ir(i),/extract)
		if temp2(0) eq 'IRAS' then temp2(0) = 'IRAS '
		aname_ir(i) = strupcase(strjoin(temp2))
	endfor
	
	;stop
	
	; Extract the fluxes
	
	a12 = irflux(indgen(nasource)*5 + 1)
	a12 = strmid(a12,27,11)
	a25 = irflux(indgen(nasource)*5 + 2)
	a25 = strmid(a25,27,11)
	a60 = irflux(indgen(nasource)*5 + 3)
	a60 = strmid(a60,27,11)
	a100 = irflux(indgen(nasource)*5 + 4)
	a100 = strmid(a100,27,11)
	
	; Extract the errors
	
	a12_err_raw = irflux(indgen(nasource)*5 + 1)
	a12_err = strmid(a12_err_raw,42,5)
	a25_err_raw = irflux(indgen(nasource)*5 + 2)
	a25_err = strmid(a25_err_raw,42,5)
	a60_err_raw = irflux(indgen(nasource)*5 + 3)
	a60_err = strmid(a60_err_raw,42,5)
	a100_err_raw = irflux(indgen(nasource)*5 + 4)
	a100_err = strmid(a100_err_raw,42,5)
	
	for i = 0, nasource - 1 do begin
		temp3_12 = strmid(a12_err(i),strlen(a12_err(i))-1,1)			; Find if error is a percentage 
		temp3_25 = strmid(a25_err(i),strlen(a25_err(i))-1,1)
		temp3_60 = strmid(a60_err(i),strlen(a60_err(i))-1,1)
		temp3_100 = strmid(a100_err(i),strlen(a100_err(i))-1,1)
		if strmid(a12_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
			a12(i) = 0
			a12_err(i) = 0 
		endif else begin
			if temp3_12 eq '%' then begin
				temp4_12 = strmid(a12_err(i),0,strlen(a12_err(i))-1)	; Convert percentage to absolute
				a12_err(i) = string(1d-2 * temp4_12 * a12(i),format='(f5.3)')
			endif
		endelse
		if strmid(a25_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
			a25(i) = 0
			a25_err(i) = 0 
		endif else begin
			if temp3_25 eq '%' then begin
				temp4_25 = strmid(a25_err(i),0,strlen(a25_err(i))-1)
				a25_err(i) = string(1d-2 * temp4_25 * a25(i),format='(f5.3)')
			endif
		endelse
		if strmid(a60_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
			a60(i) = 0
			a60_err(i) = 0 
		endif else begin
			if temp3_60 eq '%' then begin
				temp4_60 = strmid(a60_err(i),0,strlen(a60_err(i))-1)
				a60_err(i) = string(1d-2 * temp4_60 * a60(i),format='(f5.3)')
			endif
		endelse
		if strmid(a100_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
			a100(i) = 0
			a100_err(i) = 0 
		endif else begin
			if temp3_100 eq '%' then begin
				temp4_100 = strmid(a100_err(i),0,strlen(a100_err(i))-1)
				a100_err(i) = string(1d-2 * temp4_100 * a100(i),format='(f5.3)')
			endif
		endelse
	endfor
	
	
	;a_lir = lir(a12,a25,a60,a100,a_dl)
	;a_lfir = lfir(a60,a100,a_dl)
	
	save,aname_ir,a12,a25,a60,a100,a12_err,a25_err,a60_err,a100_err,filename=writefile(j)

endfor	; j

end
