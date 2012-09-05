function hrwavelength, ion
;+
; NAME:
;       
;	HRWAVELENGTH
;
; PURPOSE:
;
;	Function returns double-precision value of rest wavelength for common mid-IR transitions
;
; INPUTS:
;
;	ION - 		string giving name of ion (eg, 'neII')
;
; OUTPUTS:
;
;	WAVELENGTH - 	float giving rest wavelength of transition in microns
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> wavelength = hrwavelength('neII')
;
; NOTES:
;
;	Rest wavelengths are taken from the ISO line list at
;		http://www.ipac.caltech.edu/iso/lws/ir_lines.html 
;
;	HI 7-6 taken from MPE-Garching line list: http://www.mpe-garching.mpg.de/iso/linelists/Hydrogenic.html 
;
;	OIV discrepancy - 25.913 (ISO) vs. 25.8903 (SMART, Farrah07, MPE)
;
; REVISION HISTORY
;       Written by K. Willett                May 08
;	Added kluges for LR transitions - Aug 08
;-


; Remove suffix for LR transitions

if strmid(ion,strlen(ion)-3,3) eq '_lr' then ion_short = strmid(ion,0,strlen(ion)-3) else ion_short = ion

; Fix line center for the ArII / H2 S(5) complex

if ion eq 'arIIh2s5_lr' then ion_short = 'arII'

; Lookup arrays

lines = double([5.51118, 6.10857, 6.90952, 6.985274, $
	8.02505, 8.99103, 9.66492, 10.5105, 12.27861, 12.371898, 12.81355, $
	13.7, 14.0, 14.322, 14.3678, 15.0, 15.555, 17.03484, 17.9363, 18.7129, 24.318, 25.8903, $
	25.9882, 28.21883, 33.481, 34.815])

line_id = ['H2S7', 'H2S6', 'H2S5', 'ArII', $
	'H2S4', 'ArIII', 'H2S3', 'SIV', 'H2S2', 'HI76', 'NeII', $
	'C2H2', 'HCN', 'NeV', 'ClII', 'CO2', 'NeIII', 'H2S1', 'FeII', 'SIII', 'NeV24', 'OIV', $
	'FeII26', 'H2S0', 'SIII33', 'SiII34']

if n_elements(lines) ne n_elements(line_id) then begin
	message,'List of wavelengths and line IDs have different numbers of elements.'
	stop
endif

; Find match to input string

matchind = where(strupcase(line_id) eq strupcase(ion_short))

; Error messages

if n_elements(matchind) ne 1 then begin
	message,'Multiple transitions found'
	stop
endif else if matchind(0) eq -1 then begin
	message,'Ion not found'
	stop
endif else return, lines(matchind)

end
