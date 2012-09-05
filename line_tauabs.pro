;+
; NAME:
;       
;	LINE_TAUABS
;
; PURPOSE:
;
;	Return the optical depth for hi-res absorption lines measured in SMART
;
; INPUTS:
;
;	FILE - 		string giving path of text file containing line measurement data
;
; OUTPUTS:
;
;	TAU - 		optical depth of the absorption feature
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> result = line_tauabs('mega001_neII_lines.txt')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;-

function line_tauabs, file, float=float

tag, file, dirtag
dir = '~/Astronomy/Research/Spitzer/'+dirtag+'/lines/nosky/hires/'
fullfile = dir+file

nlines = file_lines(fullfile)
emptyarr = strarr(nlines)

openr, lun, fullfile, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

; Baseline flux

baseflux = where(strmid(emptyarr, 0, 8) eq 'Midpoint')

if n_elements(baseflux) ne 1 then begin
	message,'ERROR - more than one value found for baseline flux'
	stop
endif

baseflux_string = strmid(emptyarr(baseflux),32,14)
baseflux = float(baseflux_string)

; Line depth below baseline

absdepth = where(strmid(emptyarr, 0, 12) eq ' Line Height')

if n_elements(absdepth) ne 1 then begin
	message,'ERROR - more than one value found for line depth'
	stop
endif

absdepth_string = strmid(emptyarr(absdepth),32,11)
absdepth = float(absdepth_string)

tau = -1d * alog((baseflux + absdepth) / baseflux)
tau_str = string(tau,format='(f5.2)')

if keyword_set(float) then return, tau else return, tau_str

end
