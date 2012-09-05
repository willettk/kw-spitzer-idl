function getlinelist, ion, ohm=ohm, arch=arch, control=control, cso = cso, stop = stop, hrsky = hrsky
;+
; NAME:
;       
;	GETLINELIST
;
; PURPOSE:
;
;	Print list of line fluxes from LINEHR.pro
;
; INPUTS:
;
;
;
; OUTPUTS:
;
;
;
; KEYWORDS:
;
;	CONTROL - 	displays data for control sample galaxies
;
;	CSO - 		displays data for CSOs
;
; EXAMPLE:
;
;	IDL> list = getlinelist('neII')
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Dec 08
; 	Eliminated ARCH keyword (all OHMs displayed as the default) - Jan 10
;-

; Restore the data

if keyword_set(hrsky) then begin
	linehr, /quiet, /hrsky 
	skydir = 'hrsky/'
endif else begin
	linehr, /quiet
	skydir = ''
endelse

;if keyword_set(arch) then begin
;	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
;	restore,dir+ion+'.sav'
;	tags = line.tag
;	flux = line.flux
;	archind = where(strmid(tags,0,4) eq 'arch')
;	if archind[0] ne -1 then begin
;		tags = tags[archind]
;		flux = flux[archind]
;	endif else message,'No data found for '+ion+' in archived files'
;endif else if keyword_set(control) then begin
if keyword_set(control) then begin
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
	conind = where(strmid(tags,0,3) eq 'con')
	if conind[0] ne -1 then begin
		tags = tags[conind]
		flux = flux[conind]
	endif else message,'No data found for '+ion+' in control files'
endif else if keyword_set(cso) then begin
	dir = '~/Astronomy/Research/Spitzer/cso/linedata/'
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
endif else begin
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
	ohmind = where(strmid(tags,0,4) eq 'mega' or strmid(tags,0,4) eq 'arch')
	if ohmind[0] ne -1 then begin
		tags = tags[ohmind]
		flux = flux[ohmind]
	endif else message,'No data found for '+ion+' in mega files'
endelse


tagarr = transpose(tags)
fluxarr = transpose([string(flux * 1d21,format='(f6.2)')])
;print,''
;print,'Object   Flux                              '+string(n_elements(tagarr))+' objects'
;print,'       [10^-21 W/cm^2]'
;print,''
;print,[tagarr,fluxarr]

if keyword_set(stop) then stop

thing = [tagarr,fluxarr]
return, thing

end
