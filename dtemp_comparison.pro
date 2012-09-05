
;+
; NAME:
;       
;	DTEMP_COMPARISON
;
; PURPOSE:
;
;	Compare the dust temperatures from Yun & Carilli analytical fit to PAHFIT primary component
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
;       Written by K. Willett                Mar 09
;-


; Load data from PAHFIT output

megafiles = findfile('~/Astronomy/Research/Spitzer/ohm/pahfit/mega*.txt')
archfiles = findfile('~/Astronomy/Research/Spitzer/archived/pahfit/arch*.txt')
confiles  = findfile('~/Astronomy/Research/Spitzer/archived/pahfit/arch*.txt')

ohmfiles = [megafiles, archfiles]

for k = 0, 1 do begin

	if k eq 0 then begin
		nfiles = n_elements(ohmfiles)
		files = ohmfiles
		dtemparr = fltarr(nfiles)
		dtauarr  = fltarr(nfiles)
		tagarr   = strarr(nfiles)
		galtype = 'OHM'
	endif else begin
		nfiles = n_elements(confiles)
		files = confiles
		dtemparr = fltarr(nfiles)
		dtauarr  = fltarr(nfiles)
		tagarr   = strarr(nfiles)
		galtype = 'control'
	endelse

	for i = 0, nfiles - 1 do begin
	
		; Open file and read in data
	
		file = files[i]
		emptyarr = strarr(file_lines(file))
	
		openr, lun, file, /get_lun
		readf, lun, emptyarr
		close, lun 
		free_lun, lun
	
		; Find tag of object involved
	
		j1 = strsplit(file,'/',/ex)
		j2 = j1[n_elements(j1) - 1]
		j3 = strsplit(j2,/ex,'_')
		j4 = j3[0]
	
		tag = j4
	
		; Extract dust temp information
	
		dustlines = where(strmid(emptyarr, 0, 8) eq '  T_dust', count)
	
		dtemp = strarr(count)
		dtau  = strarr(count)
	
		; Extract dust temperature and tau
	
		for j = 0, count - 1 do begin
	
			dtemp[j] = strmid(emptyarr[dustlines[j]], 10, 10)
			dtau[j]  = strmid(emptyarr[dustlines[j]], 30, 9)
		endfor
	
		dtemp = float(dtemp) & dtau = float(dtau)
	
		; Find dominant component
	
		mostdust = where(dtau eq max(dtau))
	
		dtemparr[i] = dtemp[mostdust]
		dtauarr[i]  = dtau[mostdust]
		tagarr[i]   = tag
	
	endfor

	print,''
	print,'Avg. ',galtype,' maximum dust tau is at: ',mean(dtemparr),' +- ', stddev(dtemparr),' K'
	print,''

endfor

end
