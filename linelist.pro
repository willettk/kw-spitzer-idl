pro linelist, ion, linearr, control=control, cso = cso, stop = stop, hrsky = hrsky, quiet = quiet, error = error
;+
; NAME:
;       
;	LINELIST
;
; PURPOSE:
;
;	Print list of line fluxes from LINEHR.pro
;
; INPUTS:
;
;	ION - 		string of ion to list (eg, 'neII')
;
; OUTPUTS:
;
;	LINEARR - 	2 x N string array with the tags and fluxes
;
; KEYWORDS:
;
;	CONTROL - 	displays data for control sample galaxies
;
;	CSO - 		displays data for CSOs
;
;	ERROR - 	display error in flux measurement from SMART
;
; EXAMPLE:
;
;	IDL> linelist,'neII'
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Jul 08
; 	Merged ARCH/MEGA, added output array LINEARR - Feb 09
;	Added ERROR keyword - Sep 09
;-

; Restore the data

if keyword_set(hrsky) then begin
	linehr, /quiet, /hrsky 
	skydir = 'hrsky/'
endif else begin
	linehr, /quiet
	skydir = ''
endelse

if keyword_set(control) then begin
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
	fluxerr = line.fluxerr
	conind = where(strmid(tags,0,3) eq 'con')
	if conind[0] ne -1 then begin
		tags = tags[conind]
		flux = flux[conind]
		fluxerr = fluxerr[conind]
	endif else message,'No data found for '+ion+' in control files'
endif else if keyword_set(cso) then begin
	dir = '~/Astronomy/Research/Spitzer/cso/linedata/'
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
	fluxerr = line.fluxerr
endif else begin
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
	fluxerr = line.fluxerr
	ohmind = where(strmid(tags,0,4) eq 'mega' or strmid(tags,0,4) eq 'arch', ocount)

	if ocount gt 0 then begin
		tags = tags[ohmind]
		flux = flux[ohmind]
		fluxerr = fluxerr[ohmind]
	endif else message,'No data found for '+ion+' in megamaser files'

endelse

tagarr = transpose(tags)
fluxarr = transpose([string(flux * 1d21,format='(f6.2)')])
fluxerrarr = transpose([string(fluxerr * 1d21,format='(f6.2)')])
pmarr = transpose(replicate('   +/- ',n_elements(tagarr)))

if keyword_set(error) then begin
	linearr = [tagarr, fluxarr, pmarr, fluxerrarr]
	line1 = 'Object   Flux           Flux err.           '+string(n_elements(tagarr))+' objects'
endif else begin
	linearr = [tagarr, fluxarr]
	line1 = 'Object   Flux                                 '+string(n_elements(tagarr))+' objects'
endelse

if not keyword_set(quiet) then begin
	print,''
	print, line1
	print,'       [10^-21 W/cm^2]'
	print,''
	print, linearr
endif

if keyword_set(stop) then stop
end
