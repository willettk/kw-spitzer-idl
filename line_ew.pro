;+
; NAME:
;       
;	LINE_EW
;
; PURPOSE:
;
;	Return the equivalent width of IRS HR lines measured with SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	EW - 		equivalent width of line [microns]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_ew('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;-

function line_ew, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

ew_line = where(strmid(emptyarr, 0, 16) eq 'Equivalent Width')

if n_elements(ew_line) ne 1 then begin
	message,'ERROR - more than one value found for EW'
	stop
endif

ew_string = strmid(emptyarr(ew_line),30,12)
ew = float(ew_string)

return, ew

end
