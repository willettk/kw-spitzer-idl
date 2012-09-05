pro linehr, hrsky = hrsky, cso = cso, quiet = quiet
;+
; NAME:
;       
;	LINEHR
;
; PURPOSE:
;
;	Read HR line measurements exported from SMART into IDL and save to files
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
;	HRSKY - 	computes values for lines with HR sky subtracted; default is for no sky subtraction
;
;	QUIET - 	do not print results to screen
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	Lines measured as multiple Gaussians should be manually split into two files. Need to do this for both HR sky and no
;		sky removal for all files. 
;
; REVISION HISTORY
;       Written by K. Willett                
;	Added absorption lines - Feb 09
;-

; List of the transitions measured in HR modules

ionlist = [$
	'neII', $
	'neIII', $
	'sIII', $
	'sIII33', $
	'sIV', $
	'siII34', $
	'h2s0', $
	'h2s1', $
	'h2s2', $
	'h2s3', $
	'h2s4', $
	'h2s1_lr', $
	'h2s3_lr', $
	'h2s4_lr', $
	'arIIh2s5_lr', $
	'h2s6_lr', $
	'h2s7_lr', $
	'arIII', $
	'clII', $
	'hI76', $
	'feII', $
	'feII26', $
	'oIV', $
	'neV', $
	'neV24', $
	'c2h2', $
	'hcn', $
	'co2']

; Run modules on each and record data

nion = n_elements(ionlist)

; Loop over each transition

for i = 0, nion - 1 do begin

	if strmid(ionlist[i],strlen(ionlist[i])-2,2) eq 'lr' then begin				; LR modules
		megadir = '~/Astronomy/Research/Spitzer/ohm/lines/nosky/lores/'			
		archdir = '~/Astronomy/Research/Spitzer/archived/lines/nosky/lores/'
		controldir = '~/Astronomy/Research/Spitzer/control/lines/nosky/lores/'
		csodir = '~/Astronomy/Research/Spitzer/cso/lines/hrsky/lores/'
	endif else begin

		absorb = 0
		if ionlist[i] eq 'co2' or ionlist[i] eq 'hcn' or ionlist[i] eq 'c2h2' then absorb = 1

		if absorb eq 1 then begin

			if keyword_set(hrsky) then begin						; HR modules
				megadir = '~/Astronomy/Research/Spitzer/ohm/lines/hrsky/hires/abs/'		
				archdir = '~/Astronomy/Research/Spitzer/archived/lines/hrsky/hires/abs/'
				controldir = '~/Astronomy/Research/Spitzer/control/lines/hrsky/hires/'
				skydir = 'hrsky/'
			endif else begin
				megadir = '~/Astronomy/Research/Spitzer/ohm/lines/nosky/hires/abs/'			
				archdir = '~/Astronomy/Research/Spitzer/archived/lines/nosky/hires/abs/'
				controldir = '~/Astronomy/Research/Spitzer/control/lines/nosky/hires/abs/'
				skydir = ''
			endelse
			
			csodir = '~/Astronomy/Research/Spitzer/cso/lines/hrsky/hires/'
	
		endif else begin

			if keyword_set(hrsky) then begin						; HR modules
				megadir = '~/Astronomy/Research/Spitzer/ohm/lines/hrsky/hires/'		
				archdir = '~/Astronomy/Research/Spitzer/archived/lines/hrsky/hires/'
				controldir = '~/Astronomy/Research/Spitzer/control/lines/hrsky/hires/'
				skydir = 'hrsky/'
			endif else begin
				megadir = '~/Astronomy/Research/Spitzer/ohm/lines/nosky/hires/'			
				archdir = '~/Astronomy/Research/Spitzer/archived/lines/nosky/hires/'
				controldir = '~/Astronomy/Research/Spitzer/control/lines/nosky/hires/'
				skydir = ''
			endelse
			
			csodir = '~/Astronomy/Research/Spitzer/cso/lines/hrsky/hires/'
	
		endelse

	endelse

	megafiles = file_search(megadir+'*_'+ionlist(i)+'_lines.txt')
	archfiles = file_search(archdir+'*_'+ionlist(i)+'_lines.txt')
	confiles  = file_search(controldir+'*_'+ionlist(i)+'_lines.txt')
	csofiles  = file_search(csodir+'*_'+ionlist(i)+'_lines.txt')

	; Combine all OHM data

	allfiles = [megafiles,archfiles,confiles]
	if keyword_set(cso) then allfiles = csofiles
	nofileind = where(allfiles ne '')

	; Record data in structure for transition

	if nofileind(0) ne -1 then begin

		allfiles = allfiles(nofileind)
		nfiles = n_elements(allfiles)

		line = {ion:'', $
			restwave:0d, $
			tag:strarr(nfiles), $
			fwhm:fltarr(nfiles), $
			fwhm_vel:fltarr(nfiles), $
			fwhm_vel_err:fltarr(nfiles), $
			ew:fltarr(nfiles), $
			flux:fltarr(nfiles), $
			fluxerr:fltarr(nfiles), $
			mcenter:fltarr(nfiles)}

		for j = 0, nfiles - 1 do begin
		
			t1 = reverse(strsplit(allfiles[j],'/',/ex))
			t2 = strsplit(t1(0),'.',/ex)
			t3 = strsplit(t2(0),'_',/ex)

			tag = t3(0)
			
			if strmid(ionlist[i],strlen(ionlist[i])-2,2) eq 'lr' then begin		
				if ionlist[i] eq 'arIIh2s5_lr' then ion = 'arII' else $
					ion = strmid(ionlist[i],0,strlen(ionlist[i])-3)
			endif else ion = ionlist[i]

			line.ion = ion
			line.restwave = hrwavelength(ion)

			line.fwhm[j] = line_fwhm(allfiles[j])
				fwhm_vel_temp = line_fwhm_vel(allfiles[j])
			line.fwhm_vel[j] = fwhm_vel_temp[0]
			line.fwhm_vel_err[j] = fwhm_vel_temp[1]
			line.ew[j] = line_ew(allfiles[j])
			line.flux[j] = line_flux(allfiles[j])
			line.fluxerr[j] = line_fluxerr(allfiles[j])
			line.mcenter[j] = line_mcenter(allfiles[j])
			line.tag[j] = tag

		endfor

		; Print confirmation to screen

		if not keyword_set(quiet) then begin
			print,''
			if nfiles gt 1 then print,ionlist[i],nfiles,' detections' else print,ionlist[i],nfiles,' detection'
		endif

		if keyword_set(cso) $
			then linesavdir = '~/Astronomy/Research/Spitzer/cso/linedata/'+skydir $
			else linesavdir = '~/Astronomy/Research/Spitzer/linedata/'+skydir

		save, line, filename = linesavdir+ionlist[i]+'.sav'

	; Message if no files are found

	endif else if not keyword_set(quiet) then begin
		print,''
		print,'No line measurements found for '+ionlist[i]
		print,''
		wait,2
	endif

endfor

if not keyword_set(quiet) then begin
	if keyword_set(cso) then begin
		ncso = n_elements(csodat('tag'))
		print,''
		print,strtrim(ncso,2)+' CSOs'
	endif else begin
		nmega = n_elements(ohmdat('tag'))
		narch = n_elements(archdat('tag'))
		ncon  = n_elements(condat('tag'))
		print,''
		print,strtrim(nmega,2)+' OHMs (Darling)'
		print,strtrim(narch,2)+' OHMs (archived)'
		print,strtrim(ncon,2)+' control sample'
		print,'Total of '+strtrim(nmega+narch+ncon,2)+' objects'
		print,''
	endelse
endif

end
