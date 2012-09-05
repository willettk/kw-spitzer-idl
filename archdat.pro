function archdat, field, names = names, exp=exp, float = float, int = int, verbose = verbose, obj = obj, table = table
;+
; NAME: 
;       ARCHDAT 
;
; PURPOSE:
;
;	Call data from IDL structures created by ARCHMAKE
;
; CATEGORY:
;	ASTRONOMY; DATABASE
;
; INPUTS:
;
;	FIELD - 		string giving the tag for the requested data
;
; OUTPUTS:
;
;	SZ x N array of data, where SZ is the specified size of the field and N is the number of targets
;
; REQUIRES:
;
;	ARCHMAKE must have been run and the IDL structures exist in ~/Astronomy/Research/Spitzer/archived/data/structures/
;
; EXAMPLE:
;
;	IDL> print, archdat('sil')	; Outputs a 2 x N array where the first col. is the sil strength and the second the error
;
; NOTES:
;	Will not work for retrieving spectral data, since each spectrum has a different number of elements and so cannot
;		all be placed in a single array.
;
; MODIFICATION HISTORY:
;
;	Adapted from OHMDAT.pro - KW, 11 Mar 08
; 	Added format keywords - Sep 08
;-

; Define data location
dir = '~/Astronomy/Research/Spitzer/archived/data/structures/'

; Locate files
cd,dir
files = file_search(dir+'arch*.sav')

; Set size of array
count = n_elements(files)
restore,'~/Astronomy/Research/Spitzer/archived/data/structures/sedsize.sav'

if keyword_set(names) then print,transpose(tag_names(sedsize))

if n_elements(field) gt 0 then begin
	szindex = where(tag_names(sedsize) eq strupcase(field))
	sz = sedsize.(szindex[0])
	temparr = strarr(count,sz)
	
	; Read in data to array
	for i = 0, count - 1 do begin
		restore, files[i]
		index = where(tag_names(sed) eq strupcase(field), jj)
		if keyword_set(exp) then entry = string(sed.(index[0]),format='(e9.2)') else $
			if keyword_set(float) then entry = string(sed.(index[0]),format='(f9.2)') else $
			if keyword_set(table) then entry = string(sed.(index[0]),format='(f9.1)') else $
			if keyword_set(int) then entry = string(sed.(index[0]),format='(i5)') else $
			entry = sed.(index[0])
		if jj ne 0 then temparr[i,*] = entry else print,'Field not found'
	endfor
	
	; Return array

	sidearr = transpose(temparr)
	if keyword_set(verbose) then begin
		a = archdat('tag')
		s = [a,sidearr]
		return, s
	endif else if keyword_set(obj) then begin
		b = archdat('obj')
		bind = sort(b)
		s = [transpose(b[bind]),sidearr[*,bind]]
		return, s
	endif else return,sidearr

endif

end
