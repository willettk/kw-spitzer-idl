pro linelumlist, ion, linearr, $
	control=control, cso = cso, stop = stop, hrsky = hrsky, quiet = quiet, numonly=numonly
;+
; NAME:
;       
;	LINELUMLIST
;
; PURPOSE:
;
;	Print list of atomic line luminosities in units of L_sun
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
;	ARCH - 		displays data for archived OHMs (default shows Darling OHMs)
;
;	CONTROL - 	displays data for control sample galaxies
;
;	CSO - 		displays data for CSOs
;
; EXAMPLE:
;
;	IDL> linelumlist,'neII'
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Aug 09
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
	dlarr = condat('dl',/v)
	dl = dlarr[1,*]
	dltag = dlarr[0,*]
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
	dlarr = csodat('dl',/v)
	dl = dlarr[1,*]
	dltag = dlarr[0,*]
	dir = '~/Astronomy/Research/Spitzer/cso/linedata/'
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
endif else begin
	dlarr1 = ohmdat('dl',/v)
	dlarr2 = archdat('dl',/v)
	dl = [transpose(dlarr1[1,*]),transpose(dlarr2[1,*])]
	dltag = transpose([transpose(dlarr1[0,*]),transpose(dlarr2[0,*])])
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	flux = line.flux
	ohmind = where(strmid(tags,0,4) eq 'mega' or strmid(tags,0,4) eq 'arch', ocount)

	if ocount gt 0 then begin
		tags = tags[ohmind]
		flux = flux[ohmind]
	endif else message,'No data found for '+ion+' in megamaser files'

endelse

tagarr = transpose(tags)
fluxarr = transpose([string(flux * 1d21,format='(f6.2)')])

; Retrieve luminosity distance

match, tags, transpose(dltag), tagind, dlind, count = dltagcount

if dltagcount eq n_elements(fluxarr) then begin
	dl = dl[dlind]
	lum = alog10(4d * !dpi * (dl * 3.086d24)^2 * fluxarr * 1d-21 * 1d7 / 3.862d33)		; log L_sun
	lumarr = transpose([string(lum, format='(f6.2)')])
endif else message, 'Did not find luminosity distances for all targets'

linearr = [tagarr, lumarr]

if not keyword_set(quiet) then begin
	print,''
	print,'Object   Luminosity                        '+string(n_elements(tagarr))+' objects'
	print,'       [log L_sun]'
	print,''
	if keyword_set(numonly) then print, linearr[1,*] else print, linearr
endif

if keyword_set(stop) then stop
end
