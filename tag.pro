pro tag, specname, dirtag
;+
; NAME:
;       TAG
;
; PURPOSE:
; 	ID the directory for different Spitzer spectra
;
; OUTPUTS:
;	
;	DIRTAG - string telling calling program where to look for data (eg, 'CSO')
;
; EXAMPLE:
;	IDL> tag, 'mega001', dirtag
;
; REVISION HISTORY
;       Written by K. Willett                May 2007
;	Added ULIRG tag		- KW, Jun 07
;-

tagtemp = strmid(specname,0,3)
if tagtemp eq 'meg' then dirtag = 'ohm' else $
if tagtemp eq 'uli' then dirtag = 'ulirg' else $
if tagtemp eq 'cso'  then dirtag = 'cso' else $
if tagtemp eq 'con'  then dirtag = 'control' else $
if tagtemp eq 'arc'  then dirtag = 'archived' else begin
	on_error,2
	message,'Did not find directory for Spitzer data'
endelse

end
