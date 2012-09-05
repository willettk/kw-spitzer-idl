
;+
; NAME:
;       
;	DUSTY_TABLES
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

restore,'~/Astronomy/Research/Spitzer/dusty_pahfit_results.sav'

nohm = n_elements(ohmtags)
ohmobj = strarr(nohm)

for i=0,nohm - 1 do begin
	targets, ohmtags[i], r, o, d
	ohmobj[i] = o
endfor
	

ncon = n_elements(contags)
conobj = strarr(ncon)

for i=0,ncon - 1 do begin
	targets, contags[i], r, o, d
	conobj[i] = o
endfor
	
ohmind = sort(ohmobj)
ohm_yarr = string(ohm_yarr,format='(i4)')
ohm_qarr = string(ohm_qarr,format='(f3.1)')
ohm_tauvarr = string(ohm_tauvarr,format='(f5.1)')
ohm_tdustarr = string(ohm_tdustarr,format='(i3)')

print,[transpose(ohmobj[ohmind]),transpose(ohm_yarr[ohmind]),transpose(ohm_qarr[ohmind]),transpose(ohm_tauvarr[ohmind]),transpose(ohm_tdustarr[ohmind])]

conind = sort(conobj)
con_yarr = string(con_yarr,format='(i4)')
con_qarr = string(con_qarr,format='(f3.1)')
con_tauvarr = string(con_tauvarr,format='(f5.1)')
con_tdustarr = string(con_tdustarr,format='(i3)')

print,''
print,[transpose(conobj[conind]),transpose(con_yarr[conind]),transpose(con_qarr[conind]),transpose(con_tauvarr[conind]),transpose(con_tdustarr[conind])]

end
