
;+
; NAME:
;       
;	H2MASS_BATCH
;
; PURPOSE:
;
;	Print list of warm, hot H2 masses derived from H2 S(1) lines and T_ex for all OHMs and non-masing galaxies
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
;       Written by K. Willett                Mar 10
;-

anames = archdat('tag')
onames = ohmdat('tag')

aobj = archdat('obj')
oobj = ohmdat('obj')

ah2temp = archdat('h2temp')
oh2temp = ohmdat('h2temp')

ohmtags = [transpose(anames),transpose(onames)]
ohmobjs = [transpose(aobj),transpose(oobj)]
ohmh2temp     = [transpose(ah2temp[0,*]),transpose(oh2temp[0,*])]
ohmh2temp_hot = [transpose(ah2temp[1,*]),transpose(oh2temp[1,*])]

ohmtags = ohmtags[sort(ohmobjs)]
ohmh2temp     = ohmh2temp[sort(ohmobjs)]
ohmh2temp_hot = ohmh2temp_hot[sort(ohmobjs)]

for i=0, n_elements(ohmtags) - 1 do begin
	
	targets, ohmtags[i], r, o, dl
	h2s1_flux = getlineflux(ohmtags[i], 'h2s1',/quiet)	
	h2s3_flux = getlineflux(ohmtags[i], 'h2s3',/quiet)	
	tex = ohmh2temp[i]
	tex_hot = ohmh2temp_hot[i]
	
	if tex gt 0. then begin
		mass_temp      = h2mass(h2s1_flux[0], tex,      dl, 1,/msun,/quiet)
	endif else begin
		mass_temp_hot = 0
	endelse
	
	if tex_hot gt 0. then begin
		mass_temp_hot  = h2mass(h2s3_flux[0], tex_hot,  dl, 3,/msun,/quiet)
	endif else begin
		mass_temp_hot = 0
	endelse
	
;	if mass_temp eq 0 then print,o,'  --' else print,o,' ',string(mass_temp, format='(f5.2)')
	if mass_temp_hot eq 0 then print,o,'  --' else print,o,' ',string(mass_temp_hot, format='(f6.2)')
	
endfor

print,''

contags = condat('tag',/obj)
conh2 = condat('h2temp',/obj)

contags = contags[1,*]
conh2temp = conh2[1,*]
conh2temp_hot = conh2[2,*]

for i=0, n_elements(contags) - 1 do begin
	
	targets, contags[i], r, o, dl
	h2s1_flux = getlineflux(contags[i], 'h2s1',/quiet)	
	h2s3_flux = getlineflux(contags[i], 'h2s3',/quiet)	
	tex = conh2temp[i]
	tex_hot = conh2temp_hot[i]
	
	if tex gt 0. then begin
		mass_temp      = h2mass(h2s1_flux[0], tex,      dl, 1,/msun,/quiet)
	endif else begin
		mass_temp_hot = 0
	endelse
	
	if tex_hot gt 0. then begin
		mass_temp_hot  = h2mass(h2s3_flux[0], tex_hot,  dl, 3,/msun,/quiet)
	endif else begin
		mass_temp_hot = 0
	endelse
	
;	if mass_temp eq 0 then print,o,'  --' else print,o,' ',string(mass_temp, format='(f5.2)')
	if mass_temp_hot eq 0 then print,o,'  --' else print,o,' ',string(mass_temp_hot, format='(f6.2)')
	
endfor


end

