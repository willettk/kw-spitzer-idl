;+
; NAME:
;       
;	CRYSTALLINE_STATS
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
;       Written by K. Willett                Oct 09
;-

; OHMs with detected crystalline silicates

csil_tags = ['mega'+string([8,9,26,27,32],format='(i03)'),'arch'+string([4,5,7,9,10,12,18,29,30,32,45],format='(i03)')]

csil16 = [0.16,0.14,0.25,0.27,0.24,0.14,0.30,0.19,0.37,0.41,0.25,0.18,0.13,0.20,0.19,0.25]
csil23 = [0.09,0.13,0.15,0.03,0.20,0.12,0.32,0.09,0.19,0.25,0.16,0.12,0.10,0.14,0.12,0.19]

n = n_elements(csil_tags)

sil10 = fltarr(n)
obj = strarr(n)

for i = 0, n-1 do begin
	
	megasil = ohmdat('sil',/ver)
	archsil = archdat('sil',/ver)
	allsil = [transpose(megasil[1,*]), transpose(archsil[1,*])]
	alltag = [transpose(megasil[0,*]), transpose(archsil[0,*])]
	allobj = [transpose(ohmdat('obj')),transpose(archdat('obj'))]

	ind = where(alltag eq csil_tags[i], count)

	if count eq 1 then begin
		sil10[i] = -1d * float(allsil[ind]) 
		obj[i] = allobj[ind]
	endif else print,'Did not find 10 um silicate depth'

endfor

!p.multi=[0,1,2]

plot, sil10, csil16, $
	xtitle='10 um amorphous', $
	ytitle='16 um crystalline', $
	psym = symcat(16)

;xyouts, sil10, csil16, obj, /data

plot, sil10, csil23, $
	xtitle='10 um amorphous', $
	ytitle='23 um crystalline', $
	psym = symcat(16)

print,csil16 / sil10
range, (csil16/sil10)

print,1.5 * csil16 / sil10
range, 1.5 * (csil16/sil10)

end
