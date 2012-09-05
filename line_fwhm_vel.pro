;+
; NAME:
;       
;	LINE_FWHM_VEL
;
; PURPOSE:
;
;	Return the full-width half-max (in km/s) of IRS HR lines measured with SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	FWHM - 		full-width half-max of line [km/s]
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_fwhm_vel('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 09
;-

function line_fwhm_vel, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

fwhm_vel_line = where(strmid(emptyarr, 0, 11) eq ' Line FWHM:')
fwhm_line = where(strmid(emptyarr, 0, 24) eq ' Line FWHM (2.354*sigma)')

if n_elements(fwhm_vel_line) ne 1 then begin
	message,'ERROR - more than one value found for FWHM'
	stop
endif

fwhm_string = strmid(emptyarr(fwhm_vel_line),33,11)
fwhm_vel = float(fwhm_string)

fwhm_um = float(strmid(emptyarr(fwhm_line),33,11))
fwhm_um_err = float(strmid(emptyarr(fwhm_line),48,11))
fwhm_vel_err = fwhm_um_err / fwhm_um * fwhm_vel

return, [fwhm_vel, fwhm_vel_err]

end
