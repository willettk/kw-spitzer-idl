function getlineflux, tag, ion, quiet = quiet
;+
; NAME:
;       
;	GETLINEFLUX
;
; PURPOSE:
;
;	Return the flux for a specific object and transition
;
; INPUTS:
;
;	TAG - 		identifier of target (eg, 'mega001')
;
;	ION - 		line transition (eg, 'neII')
;
; OUTPUTS:
;
;	Two-element array with [flux, fluxerror]
;
; KEYWORDS:
;
;	QUIET - 	does not print results to screen
;
; EXAMPLE:
;
;	IDL> result = getlineflux('mega001','neII')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;	Added CSO keyword - Aug 08
;-

; Restore the data

if strmid(tag,0,3) eq 'cso' then begin
	linehr, /quiet, /cso
	dir = '~/Astronomy/Research/Spitzer/cso/linedata/'
endif else begin
	linehr, /quiet
	dir = '~/Astronomy/Research/Spitzer/linedata/'
endelse

restore,dir+ion+'.sav'

taglist = line.tag
fluxlist = line.flux
fluxerrlist = line.fluxerr

ind = where(taglist eq tag)
if ind[0] ne -1 then begin
	flux = dblarr(2)
	flux[0] = fluxlist[ind]
	flux[1] = fluxerrlist[ind]
endif else begin
	if not keyword_set(quiet) then begin
		print,''
		message,'No flux found for '+ion,/info
	endif
	flux = -1
endelse

return, flux

end
