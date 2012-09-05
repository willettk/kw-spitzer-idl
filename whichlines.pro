function whichlines, obj, hrsky = hrsky, absorb = absorb
;+
; NAME:
;       
;	WHICHLINES
;
; PURPOSE:
;
;	Return list of detected HR lines for each object
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
;	IDL> result = whichlines('mega001')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Jan 09
;-

if keyword_set(hrsky) then begin
	linehr, /quiet, /hrsky 
	skydir = 'hrsky/'
endif else begin
	linehr, /quiet
	skydir = 'nosky/'
endelse

; ID the directory

tag, obj, dirtag

if dirtag eq 'CSO' then skydir = 'hrsky/'		; All CSOs only have the HRSKY option

; Option for detecting absorption features

if keyword_set(absorb) then absdir = 'abs/' else absdir = ''

dir = '~/Astronomy/Research/Spitzer/'+dirtag+'/lines/'+skydir+'hires/'+absdir

	cd,dir
	files = file_search('*'+obj+'*lines.txt',count=fcount)
	if fcount gt 0 then begin
		thelines = strarr(fcount)
		for i = 0, fcount - 1 do begin
			t1 = strsplit(files[i],'_',/extract)
			thelines[i] = t1[1]
		endfor
	
	return, thelines
	endif else return, 0

end
