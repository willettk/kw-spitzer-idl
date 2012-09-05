function csodat, field, names = names, exp=exp, float = float, verbose = verbose, int = int, obj = obj, rasort = rasort
;+
; NAME: 
;       CSODAT 
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
;	ARCHMAKE must have been run and the IDL structures exist in ~/Astronomy/Research/Spitzer/cso/data/structures/
;
; EXAMPLE:
;
;	IDL> print, csodat('sil')	; Outputs a 2 x N array where the first col. is the sil strength and the second the error
;
; NOTES:
;	Will not work for retrieving spectral data, since each spectrum has a different number of elements and so cannot
;		all be placed in a single array.
;
; MODIFICATION HISTORY:
;
;	Adapted from OHMDAT.pro - KW, 11 Mar 08
; 	Added format keywords - Sep 08
;	Added RASORT keyword, which sorts results by RA of target (as they appear in paper) - Sep 09
;	Padded RASORT target names so that results all are in the same column - Sep 09
;-

; Define data location
dir = '~/Astronomy/Research/Spitzer/cso/data/structures/'

; Locate files
cd,dir
files = file_search(dir+'cso*.sav')

; Set size of array
count = n_elements(files)
restore,'~/Astronomy/Research/Spitzer/cso/data/structures/sedsize.sav'

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
			if keyword_set(int) then entry = string(sed.(index[0]),format='(i5)') else $
			entry = sed.(index[0])
		if jj ne 0 then temparr[i,*] = entry else print,'Field not found'
	endfor
	
	; Return array

	sidearr = transpose(temparr)
	if keyword_set(verbose) then begin
		a = csodat('tag')
		s = [a,sidearr]
	endif else if keyword_set(obj) then begin
		b = csodat('obj')
		bind = sort(b)
		s = [transpose(b[bind]),sidearr[*,bind]]
	endif else if keyword_set(rasort) then begin
		b = csodat('obj')
		bind = [1,7,0,9,5,3,4,2,8,6]
		ra_names = transpose(b[bind])
		maxlen = max(strlen(ra_names))
		for j=0, n_elements(ra_names) - 1 do $
			if strlen(ra_names[j]) lt maxlen then $
				ra_names[j] = ra_names[j] + string(replicate(32B, maxlen - strlen(ra_names[j])))
		s = [ra_names,sidearr[*,bind]]
	endif else s = sidearr

	return, s


endif

end
