
;+
; NAME:
;       
;	FARRAH
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

; Read in data from Farrah Tbl. 3

readcol,'~/Astronomy/Research/Spitzer/OHM/farrah_tbl3.txt',$
	gal, arIII, h2s3, sIV, h2s2, neII, neV14, neIII, h2s1, sIII, neV24, oIV, h2s0, sIII33, siII34, $
	format = 'a,a,a,a,a,a,a,a,a,a,a,a,a,a,a',/silent

all = [[gal], [arIII], [h2s3], [sIV], [h2s2], [neII], [neV14], [neIII], [h2s1], [sIII], $
	[neV24], [oIV], [h2s0], [sIII33], [siII34]]

sz = size(all)
lim=intarr(sz(1),sz(2))

for j = 0, sz(2)-1 do begin
	for i = 0, sz(1)-1 do begin
		cell = all(i,j)
		if strmid(cell,0,1) eq '<'  then begin
			;cell = strmid(cell,1,strlen(cell)-1)
			cell = '0.'
			lim(i,j) = 1
		endif
		if strmid(cell,0,2) eq '\l' then cell = '0.'
		;if strmid(cell,strlen(cell)-2,2) eq '::' then cell = strmid(cell,0,strlen(cell)-2)
		;if strmid(cell,strlen(cell)-1,1) eq ':'  then cell = strmid(cell,0,strlen(cell)-1)
		if strmid(cell,strlen(cell)-2,2) eq '::' then cell = '0.'				; 30% error
		if strmid(cell,strlen(cell)-1,1) eq ':'  then cell = strmid(cell,0,strlen(cell)-1)	; 20% error
		all(i,j) = cell
	endfor
endfor


allobj = [transpose(ohmdat('obj')),transpose(archdat('obj')),transpose(condat('obj'))]
alltag = [transpose(ohmdat('tag')),transpose(archdat('tag')),transpose(condat('tag'))]

all[45:52] = ['01572+0009','09320+6134','3C273','12540+5708','13428+5608','Mrk463','15327+2340','NGC6240']
farrobj = 'IRAS '+ all[*,0]

match, allobj, farrobj, a, f

matched_targets = all[f,*]
matched_tags = alltag[a]
	
; Example on NeII flux

nf = n_elements(f)
ionlist = ['ariii','h2s3','siv','h2s2','neii','nev','neiii','h2s1','siii','nev24','oiv','h2s0','siii33','siii34']

for j = 0, n_elements(ionlist) - 1 do begin

	arr = fltarr(2,nf)
	ion = ionlist[j]
	ionno = j+1

	for i = 0, nf - 1 do begin
	
		tempflux = getlineflux(matched_tags[i],ion,/quiet)
	
		arr[0,i] = tempflux[0] * 1d21
		arr[1,i] = all[f[i],ionno]
	
	endfor
	
	noline = where(arr[0,*] le 0 or arr[1,*] le 0, count)
	if count gt 0 then begin
		ind = setdifference(indgen(nf),noline)
		newarr = arr[*,ind]
	endif else newarr = arr
	
	noline_me = where(arr[0,*] le 0 and arr[1,*] gt 0, count_me)

	!p.multi=[0,2,1]
	plot, newarr[0,*], newarr[1,*], xtitle='Willett',ytitle='Farrah',title=strupcase(ion), /nodata
	oplot, newarr[0,*], newarr[1,*], psym = symcat(14), color=fsc_color("Red")
	oplot,indgen(100),linestyle=2
	plot, newarr[0,*]/newarr[1,*], psym = symcat(14),title=strupcase(ion)+' ratio',xr=[-1,n_elements(newarr[0,*])+1]
	oplot, newarr[0,*]/newarr[1,*], psym = symcat(14), color=fsc_color("Red")
	hor, 1.0, linestyle=1
	hor, mean(newarr[0,*]/newarr[1,*]), linestyle=2
	xyouts, 0.75, 0.8, string(mean(newarr[0,*]/newarr[1,*]),format='(f6.2)'), /normal

	print,ion+'    '+string(mean(newarr[0,*]/newarr[1,*]),format='(f6.2)')+'    '$
		+string(n_elements(newarr[0,*]),format='(i2)')+' objects'$
		+'     '+string(count_me)

	if count_me gt 0 then print,matched_targets[noline_me]

endfor

end
