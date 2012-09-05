
;+
; NAME:
;       
;	DUSTY_BATCH
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
;       Written by K. Willett                Nov 09
;-

;Y_arr = [1.25, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
;power_arr = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

Y_arr = [1.25, 10, 100]
power_arr = [0.0, 1.0]

dustydir = '~/Astronomy/Research/Spitzer/cso/DUSTY/'

for i = 0, n_elements(Y_arr) - 1 do begin
	
	for j=0, n_elements(power_arr) - 1 do begin

		; Retrieve template input file
		
		file = dustydir + 'kw1.inp'
		
		nlines = file_lines(file)
		emptyarr = strarr(nlines)
		
		openr, lun, file, /get_lun
		readf, lun, emptyarr
		close, lun & free_lun, lun
		
		; Set radial power (r^-p or r^-q), geometric thickness (Y)
		
		lineI_phrase = '3) Density Distribution'
		lineI = where(strmid(strtrim(emptyarr,2),0,strlen(lineI_phrase)) eq lineI_phrase)
		
		emptyarr[lineI + 3] = " - shell's relative thickness = "+string(Y_arr[i],format='(i3)')
		emptyarr[lineI + 4] = " - power = "+string(power_arr[j],format='(f3.1)')
		
		; Write new input file
		
		filename = 'kw_sphere_Y'+string(Y_arr[i],format='(i03)')+'_q'+string(powerarr[j],format='(f3.1)')
		wfile = dustydir+filename+'.inp'
		
		newarr = emptyarr
		
		openw, lun2, wfile, /get_lun
		printf, lun2, newarr
		close, lun2 & free_lun, lun2
		
		; Append new file into dusty.inp
		
		dustyinp_file = dustydir + 'dusty.inp'
		
		nd = file_lines(dustyinp_file)
		dustyarr = strarr(nd)
		
		openr, lun3, dustyinp_file, /get_lun
		readf, lun3, dustyarr
		close, lun3 & free_lun, lun3
		
		new_dustyarr = strarr(nd + 1)
		new_dustyarr[0:nd-1] = dustyarr+string(13B)
		new_dustyarr[nd] = new_dustyarr[nd - 1]
		new_dustyarr[nd-1] = filename+string(13B)
		
		openw, lun4, dustyinp_file, /get_lun
		printf, lun4, new_dustyarr
		close, lun4 & free_lun, lun4

	endfor

endfor

; Call DUSTY program

;cd, dustydir
;spawn, './dusty'

end
