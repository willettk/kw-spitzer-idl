
;+
; NAME:
;       
;	DUSTYBATCH_HIGHTAU4
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
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Feb 10
;-

starttime = systime(1)

Y_arr = [400, 500, 750, 1000]
power_arr = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
nmodels=90

dustydir = '~/Astronomy/Research/Spitzer/cso/DUSTY4/'
cd, dustydir

;if keyword_set(hightau) then begin
	tag = '_hightau' 
	specdir = 'hightau/'
;endif else begin
;	tag = ''
;	specdir = 'regtau/'
;endelse

file_copy, 'dusty_copy.inp','dusty.inp', /overwrite

for i = 0, n_elements(Y_arr) - 1 do begin
	
	for j=0, n_elements(power_arr) - 1 do begin

		; Retrieve template input file
		
		file = dustydir + 'kw1'+tag+'.inp'
		
		nlines = file_lines(file)
		emptyarr = strarr(nlines)
		
		openr, lun, file, /get_lun
		readf, lun, emptyarr
		close, lun & free_lun, lun
		
		; Set radial power (r^-p or r^-q), geometric thickness (Y)
		
		lineI_phrase = '3) Density Distribution'
		lineI = where(strmid(strtrim(emptyarr,2),0,strlen(lineI_phrase)) eq lineI_phrase)

		emptyarr[lineI + 3] = " - shell's relative thickness = "+string(Y_arr[i],format='(i4)')
		emptyarr[lineI + 4] = " - power = "+string(power_arr[j],format='(f3.1)')
		
		lineII_phrase = '- number of models'
		lineII = where(strmid(strtrim(emptyarr,2),0,strlen(lineII_phrase)) eq lineII_phrase)
		emptyarr[lineII] = strjoin(replicate(' ',8))+lineII_phrase+' = '+string(nmodels,format='(i3)')

		; Write new input file
		
		filename = 'kw_sphere_Y'+string(Y_arr[i],format='(i04)')+'_q'+string(power_arr[j],format='(f3.1)')+tag
		inpfile = dustydir+filename+'.inp'
		outfile = dustydir+filename+'.out'
		sppfile = dustydir+filename+'.spp'
		specfile_start = dustydir+filename+'.s*'
		
		newarr = emptyarr
		
		openw, lun2, inpfile, /get_lun
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

		spawn, './dusty'

		file_copy, inpfile, 'inp/'+specdir, /overwrite
		file_copy, outfile, 'out/'+specdir, /overwrite
		file_copy, sppfile, 'spp/'+specdir, /overwrite
		file_copy, specfile_start, 'spectra/'+specdir, /overwrite

		file_delete, inpfile
		file_delete, outfile
		file_delete, sppfile
		file_delete, dustydir+filename+'.s0'+string(indgen(nmodels)+1,format='(i02)')

	endfor

endfor

endtime = systime(1)

print,''
print, 'Time elapsed [min] = ',string((endtime - starttime) / 60., format='(f8.1)')
print, 'Time elapsed [hr]  = ',string((endtime - starttime) / 3600., format='(f8.1)')
print,''

end
