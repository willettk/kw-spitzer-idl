;+
; NAME:
;       
;	LINE_FLUX
;
; PURPOSE:
;
;	Return the line flux of IRS HR lines measured with SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	FLUX - 		line flux of line [W/cm^2]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_flux('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;-

function line_flux, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

flux_line = where(strmid(emptyarr, 0, 10) eq ' Line Flux')

if n_elements(flux_line) ne 1 then begin
	message,'ERROR - more than one value found for line flux'
	stop
endif

flux_string = strmid(emptyarr(flux_line),30,14)
flux = float(flux_string)

return, flux

end
