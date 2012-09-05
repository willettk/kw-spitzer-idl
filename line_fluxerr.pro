;+
; NAME:
;       
;	LINE_FLUXERR
;
; PURPOSE:
;
;	Return the flux error of IRS HR lines measured with SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	FLUX - 		flux of line [W/cm^2]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_fluxerr('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;-

function line_fluxerr, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

flux_line = where(strmid(emptyarr, 0, 10) eq ' Line Flux')

if n_elements(flux_line) ne 1 then begin
	message,'ERROR - more than one value found for line flux error'
	stop
endif

flux_string = strmid(emptyarr(flux_line),48,10)
flux = float(flux_string)

return, flux

end
