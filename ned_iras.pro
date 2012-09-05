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
;	AFILE_IR -	formatted text file containing names and photometry from the NED batch system
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
;	Added 11 new archived sources		Jul 08
;-

; Files on which to operate

spitzdir = '~/Astronomy/Research/Spitzer/'
filelist = spitzdir+['ohm/ned_mega.lst','archived/ned_arch.lst','cso/ned_cso.lst','control/ned_con.lst']
writefile = spitzdir+['ohm/data/idl_sav/ohm_iras.sav', $
	'archived/data/idl_sav/archived_iras.sav', $
	'cso/data/idl_sav/cso_iras.sav', $
	'control/data/idl_sav/control_iras.sav']

for j = 0, n_elements(filelist)-1 do begin

	openr,lun,filelist(j),/get_lun
	nobj = numlines(filelist(j))  & nasource = nobj/5
	if nobj mod 5 ne 0 then print,'Wrong number of files in list.'
	irflux = strarr(nobj)
	readf,lun,irflux
	close,lun
	
	; Extract the object names
	
	aname_ir = irflux(indgen(nasource)*5)
	aname_ir = strmid(aname_ir,9,100)
	atag_ir = strarr(n_elements(aname_ir))
	
	; Format the names
	
	objlist = [transpose(archdat('obj')),transpose(condat('obj')),$
		transpose(csodat('obj')),transpose(ohmdat('obj'))]
	objlist = strtrim(objlist,2)
	taglist = [transpose(archdat('tag')),transpose(condat('tag')),$
		transpose(csodat('tag')),transpose(ohmdat('tag'))]

	for i = 0, nasource - 1 do begin
		temp2 = strsplit(aname_ir(i),/extract)
		if n_elements(temp2) eq 2 then temp2[0] = temp2[0]+' '
		aname_ir(i) = strupcase(strjoin(temp2))

		; Kluge for a couple of misnamed targets

		if aname_ir[i] eq 'VIIZW485' then aname_ir[i] = 'VII Zw 485'
		if aname_ir[i] eq 'IRAS 17539+2935' then aname_ir[i] = 'IRAS 17540+2935'
		if aname_ir[i] eq 'IRAS 17208-0014' then aname_ir[i] = 'IRAS 17207-0014'
		
		objind = where(aname_ir[i] eq objlist, ocount)
		if ocount eq 1 then atag_ir[i] = taglist[objind] else $
			message,'Did not find match for '+aname_ir[i],/info
	endfor

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
	
	save,aname_ir,atag_ir,a12,a25,a60,a100,a12_err,a25_err,a60_err,a100_err,filename=writefile(j)

endfor	; j

end
