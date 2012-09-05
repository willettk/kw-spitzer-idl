;pro photo_control

; PHOTO_CONTROL

; Extract the IRAS fluxes from results of a NED batch file

; Adapted from NED_IRAS.pro

confile_ir='~/Astronomy/Research/Spitzer/control/ned_control.lst'
openr,lun,confile_ir,/get_lun
ncon = numlines(confile_ir)  & n_targets = ncon/5
if ncon mod 5 ne 0 then print,'Wrong number of files in list.'
conir = strarr(ncon)
readf,lun,conir
close,lun

; Extract the object names

cname_ir = conir(indgen(n_targets)*5)
cname_ir = strmid(cname_ir,10,100)

; Format the names

for i = 0, n_targets - 1 do begin
	temp2 = strsplit(cname_ir(i),/extract)
	if temp2(0) eq 'IRAS' then temp2(0) = 'IRAS '
	cname_ir(i) = strupcase(strjoin(temp2))
endfor


; Extract the fluxes

c12 = conir(indgen(n_targets)*5 + 1)
c12 = strmid(c12,27,11)
c25 = conir(indgen(n_targets)*5 + 2)
c25 = strmid(c25,27,11)
c60 = conir(indgen(n_targets)*5 + 3)
c60 = strmid(c60,27,11)
c100 = conir(indgen(n_targets)*5 + 4)
c100 = strmid(c100,27,11)

; Extract the errors

c12_err_raw = conir(indgen(n_targets)*5 + 1)
c12_err = strmid(c12_err_raw,42,5)
c25_err_raw = conir(indgen(n_targets)*5 + 2)
c25_err = strmid(c25_err_raw,42,5)
c60_err_raw = conir(indgen(n_targets)*5 + 3)
c60_err = strmid(c60_err_raw,42,5)
c100_err_raw = conir(indgen(n_targets)*5 + 4)
c100_err = strmid(c100_err_raw,42,5)

for i = 0, n_targets - 1 do begin
	temp3_12 = strmid(c12_err(i),strlen(c12_err(i))-1,1)			; Find if error is a percentage 
	temp3_25 = strmid(c25_err(i),strlen(c25_err(i))-1,1)
	temp3_60 = strmid(c60_err(i),strlen(c60_err(i))-1,1)
	temp3_100 = strmid(c100_err(i),strlen(c100_err(i))-1,1)
	if strmid(c12_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
		c12(i) = 0
		c12_err(i) = 0 
	endif else begin
		if temp3_12 eq '%' then begin
			temp4_12 = strmid(c12_err(i),0,strlen(c12_err(i))-1)	; Convert percentage to absolute
			c12_err(i) = string(1d-2 * temp4_12 * c12(i),format='(f5.3)')
		endif
	endelse
	if strmid(c25_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
		c25(i) = 0
		c25_err(i) = 0 
	endif else begin
		if temp3_25 eq '%' then begin
			temp4_25 = strmid(c25_err(i),0,strlen(c25_err(i))-1)
			c25_err(i) = string(1d-2 * temp4_25 * c25(i),format='(f5.3)')
		endif
	endelse
	if strmid(c60_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
		c60(i) = 0
		c60_err(i) = 0 
	endif else begin
		if temp3_60 eq '%' then begin
			temp4_60 = strmid(c60_err(i),0,strlen(c60_err(i))-1)
			c60_err(i) = string(1d-2 * temp4_60 * c60(i),format='(f5.3)')
		endif
	endelse
	if strmid(c100_err(i),0,1) eq ' ' then begin				; Set flux + error to zero for upper limit
		c100(i) = 0
		c100_err(i) = 0 
	endif else begin
		if temp3_100 eq '%' then begin
			temp4_100 = strmid(c100_err(i),0,strlen(c100_err(i))-1)
			c100_err(i) = string(1d-2 * temp4_100 * c100(i),format='(f5.3)')
		endif
	endelse
endfor


;c_lir = lir(c12,c25,c60,c100,c_dl)
;c_lfir = lfir(c60,c100,c_dl)

save,cname_ir,c12,c25,c60,c100,c12_err,c25_err,c60_err,c100_err,filename='~/Astronomy/Research/Spitzer/control/data/idl_sav/con_iras.sav'

;stop
end

