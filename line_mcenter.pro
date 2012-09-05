;+
; NAME:
;       
;	LINE_MCENTER
;
; PURPOSE:
;
;	Return the measured center of the best fit Gaussian for HR IRS lines
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	CENTER - 	center of best Gaussian fit to line in observed frame [microns]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_mcenter('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;-

function line_mcenter, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

mcenter_line = where(strmid(emptyarr, 0, 12) eq ' Line Center')

if n_elements(mcenter_line) ne 1 then begin
	message,'ERROR - more than one value found for line center'
	stop
endif

mcenter_string = strmid(emptyarr(mcenter_line),30,12)
mcenter = float(mcenter_string)

return, mcenter

end
