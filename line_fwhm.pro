;+
; NAME:
;       
;	LINE_FWHM
;
; PURPOSE:
;
;	Return the full-width half-max of IRS HR lines measured with SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	FWHM - 		full-width half-max of line [microns]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_fwhm('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;-

function line_fwhm, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

fwhm_line = where(strmid(emptyarr, 0, 24) eq ' Line FWHM (2.354*sigma)')

if n_elements(fwhm_line) ne 1 then begin
	message,'ERROR - more than one value found for FWHM'
	stop
endif

fwhm_string = strmid(emptyarr(fwhm_line),33,11)
fwhm = float(fwhm_string)

return, fwhm

end
