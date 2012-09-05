
;+
; NAME:
;       
;	SURV_BATCH2
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

lines = ['arIII','sIV','hI76','neII','neV','clII','neIII','h2s3','h2s2','neV24','oIV','feII26','h2s0','feII','sIII','h2s1']
;,'sIII33','siII34'] 

for i = 0, n_elements(lines) - 1 do begin

	file = '~/Astronomy/Research/Spitzer/surv/surv_res_'+lines[i]+'.txt'

	nlines = file_lines(file)
	emptyarr = strarr(nlines)

	openr, lun, file, /get_lun
	readf, lun, emptyarr
	close, lun & free_lun, lun
	

	gehan_lines = where(strmid(emptyarr, 8, 5) eq 'Gehan')

	peto_lines = where(strmid(emptyarr, 8, 4) eq 'Peto')

	logrank_line = where(strmid(emptyarr, 8, 7) eq 'Logrank')
	
	gehan_perm_stat = float(strmid(emptyarr(gehan_lines[0]+2),33,21))
	gehan_hyper_stat = float(strmid(emptyarr(gehan_lines[1]+2),33,21))
	peto_peto_stat = float(strmid(emptyarr(peto_lines[0]+2),33,21))
	peto_prentice_stat = float(strmid(emptyarr(peto_lines[1]+2),33,21))
	logrank_stat = float(strmid(emptyarr(logrank_line+2),33,21))
	
	gehan_perm_prob = float(strmid(emptyarr(gehan_lines[0]+3),33,21))
	gehan_hyper_prob = float(strmid(emptyarr(gehan_lines[1]+3),33,21))
	peto_peto_prob = float(strmid(emptyarr(peto_lines[0]+3),33,21))
	peto_prentice_prob = float(strmid(emptyarr(peto_lines[1]+3),33,21))
	logrank_prob = float(strmid(emptyarr(logrank_line+3),33,21))
	
endfor	




end
