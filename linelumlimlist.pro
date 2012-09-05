pro linelumlimlist, ion, linelimarr, $
	control=control, cso = cso, stop = stop, hrsky = hrsky, lh = lh, lores = lores, $
	quiet = quiet, numonly=numonly
;+
; NAME:
;       
;	LINELUMLIMLIST
;
; PURPOSE:
;
;	Print list of line luminosity limits in units of log L_sun
;
; INPUTS:
;
;	ION -		string indicating ion to measure (eg, 'neII')
;
; OUTPUTS:
;
;	LINELIMARR - 	2 x N string array of the limits on the line in question [tags, luminosities]
;
; KEYWORDS:
;
;	CONTROL - 	displays data for control sample galaxies (default shows data for archived and Darling OHMs)
;
;	CSO - 		displays data for CSOs
;
;	STOP - 		stop program at end
;
;	HRSKY - 	use data from HR-subtracted sky (where available)
;
; EXAMPLE:
;
;	IDL> linelumlimlist,'neII'
;
; NOTES:
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

if keyword_set(control) then begin							;;;;; CONTROL SAMPLE ;;;;;

	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag

	conind = where(strmid(tags,0,3) eq 'con',concount)
	allcon = transpose(condat('tag'))

	dl = transpose(condat('dl'))
	dltag = allcon

	if concount gt 0 then begin

		at = tags[conind]								; At least one detection
		match, at, allcon, a, b

		if n_elements(b) lt n_elements(allcon) then begin
			conlim = allcon[setdifference(indgen(n_elements(allcon)),b)]
			ncon = n_elements(conlim)
			conarr = fltarr(ncon)
			for i = 0, ncon - 1 do begin
				temp = linelim(conlim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
				conarr[i] = temp
			endfor
		endif else message,'No limits found on '+ion+' for control sample',/info
	
	endif else begin									; No detections
		conlim = allcon[indgen(n_elements(allcon))]
		ncon = n_elements(conlim)
		conarr = fltarr(ncon)
		for i = 0, ncon - 1 do begin
			temp = linelim(conlim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
			conarr[i] = temp
		endfor
	endelse

	fluxlim_arr = conarr
	tagarr = conlim

endif else if keyword_set(cso) then begin						;;;;; CSOs ;;;;;

	dir = '~/Astronomy/Research/Spitzer/cso/linedata/'+skydir
	file = dir+ion+'.sav'
	cso_exist = file_test(file)

	dl = transpose(csodat('dl'))
	dltag = transpose(csodat('tag'))
	
	if cso_exist eq 1 then begin
		restore,file
		tags = line.tag

		csoind = where(strmid(tags,0,3) eq 'cso')
		allcso = transpose(csodat('tag'))
		at = tags[csoind]
		match, at, allcso, a, b
	
		if n_elements(b) lt n_elements(allcso) then begin
			csolim = allcso[setdifference(indgen(n_elements(allcso)),b)]
			ncso = n_elements(csolim)
			csoarr = fltarr(ncso)
			for i = 0, ncso - 1 do begin
				temp = linelim(csolim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
				csoarr[i] = temp
			endfor
		endif else message,'No limits found on '+ion+' for CSOs'
	
	endif else begin                		; If no detections, then compute limits for all targets

		csolim = transpose(csodat('tag'))
		ncso = n_elements(csolim)
		csoarr = fltarr(ncso)
			for i = 0, ncso - 1 do begin
				temp = linelim(csolim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
				csoarr[i] = temp
			endfor

	endelse

	fluxlim_arr = csoarr
	tagarr = csolim

endif else begin									;;;;; OHMs ;;;;;

	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag

	dlarr1 = ohmdat('dl',/v)
	dlarr2 = archdat('dl',/v)
	dl = [transpose(dlarr1[1,*]),transpose(dlarr2[1,*])]
	dltag = [transpose(dlarr1[0,*]),transpose(dlarr2[0,*])]

	archind = where(strmid(tags,0,4) eq 'arch',acount)
	allarch = transpose(archdat('tag'))
	if acount gt 0 then begin
		at = tags[archind]
		match, at, allarch, a, b
	endif

	if n_elements(b) lt n_elements(allarch) then begin
		archlim = allarch[setdifference(indgen(n_elements(allarch)),b)]
		narch = n_elements(archlim)
		archarr = fltarr(narch)
		for i = 0, narch - 1 do begin
			temp = linelim(archlim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
			archarr[i] = temp
		endfor
	endif else message,'No limits found on '+ion+' for archived OHMs',/info

	megaind = where(strmid(tags,0,4) eq 'mega',mcount)
	allmega = transpose(ohmdat('tag'))
	if mcount gt 0 then begin

		at = tags[megaind]								; Assumes at least one detection
		match, at, allmega, a, b

		if n_elements(b) lt n_elements(allmega) then begin
			megalim = allmega[setdifference(indgen(n_elements(allmega)),b)]
			nmega = n_elements(megalim)
			megaarr = fltarr(nmega)
			for i = 0, nmega - 1 do begin
				temp = linelim(megalim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
				megaarr[i] = temp
			endfor
		endif else message,'No limits found on '+ion+' for Darling OHMs',/info
	
	endif else begin									; No detections
		megalim = allmega[indgen(n_elements(allmega))]
		nmega = n_elements(megalim)
		megaarr = fltarr(nmega)
		for i = 0, nmega - 1 do begin
			temp = linelim(megalim[i],ion,/noplot,/quiet,lores=lores,lh=lh)
			megaarr[i] = temp
		endfor
	endelse

	if n_elements(archarr) gt 0 and n_elements(megaarr) gt 0 then begin
		fluxlim_arr = [megaarr,archarr]
		tagarr = [megalim,archlim]
	endif else if n_elements(megaarr) gt 0 then begin
		fluxlim_arr = megaarr
		tagarr = megalim
	endif else begin
		fluxlim_arr = archarr
		tagarr = archlim
	endelse

endelse

; Find luminosities based on DL

match, tagarr, dltag, tagind, dlind, count = dltagcount

if dltagcount eq n_elements(fluxlim_arr) then begin
	dl = dl[dlind]
	lum_lim = alog10(4d * !dpi * (dl * 3.086d24)^2 * fluxlim_arr * 1d-21 * 1d7 / 3.862d33)		; log L_sun
	lumlim_arr = transpose([string(lum_lim, format='(f6.2)')])
endif else message, 'Did not find luminosity distances for all targets'

; Remove objects that do not lie in spectral range

notinmodule = where(lumlim_arr eq -1.0)
inmodule = where(lumlim_arr ne -1.0)

if not keyword_set(quiet) then begin
	
	print,''
	print,'Object   Luminosity limit                        '+string(n_elements(tagarr))+' objects'
	print,'       [log L_sun]'
	print,''

endif
	
	if inmodule[0] ne -1 then begin

		intagarr = tagarr[inmodule]
		influxlim_arr = lumlim_arr[inmodule]
		
		if n_elements(intagarr) gt 1 then $
			linelimarr = [transpose(intagarr),transpose(replicate(' <',n_elements(intagarr))),transpose(string(influxlim_arr,format='(f5.2)'))] else $
			linelimarr = [intagarr,' <',string(influxlim_arr,format='(f5.2)')]

	endif else if notinmodule[0] ne -1 then begin
	
		outtagarr = tagarr[notinmodule]
		outfluxlim_arr = lumlim_arr[notinmodule]
	
		if n_elements(outtagarr) gt 1 then linelimarr = [transpose(outtagarr)] else $
			linelimarr = [outtagarr]

		if not keyword_set(quiet) then begin
			print,'Objects in which line lies outside of spectral module'
			print,''

		endif

	endif

if not keyword_set(quiet) then begin
	if keyword_set(numonly) then print,linelimarr[2,*] else print,linelimarr
	print,''
endif

if keyword_set(stop) then stop
end
