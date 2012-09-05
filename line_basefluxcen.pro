;+
; NAME:
;       
;	LINE_BASEFLUXCEN
;
; PURPOSE:
;
;	Return the expected flux density at the center of the line for HR fits from SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	BASEFLUXCEN - 	expected flux density at center of line from baseline fit [W/cm^2/um]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_basefluxcen('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Feb 09 (never used)
;-

function line_basefluxcen, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

basefluxcen_line = where(strmid(emptyarr, 0, 22) eq ' Baseline Flux Density')

if n_elements(basefluxcen_line) ne 1 then begin
	message,'ERROR - more than one value found for line center'
	stop
endif

basefluxcen_string = strmid(emptyarr(basefluxcen_line),38,14)
basefluxcen = float(basefluxcen_string)

return, basefluxcen

end
