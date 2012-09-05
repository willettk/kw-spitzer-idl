pro lines_resolved, ion, linearr, control=control, cso = cso, stop = stop, hrsky = hrsky, quiet = quiet

;+
; NAME:
;       
;	LINES_RESOLVED
;
; PURPOSE:
;
;	Find out if HR lines are resolved in HR spectra
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
;	
;
; EXAMPLE:
;
;	IDL> 
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
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	fwhm_vel = line.fwhm_vel
	fwhm_vel_err = line.fwhm_vel_err
	conind = where(strmid(tags,0,3) eq 'con')
	if conind[0] ne -1 then begin
		tags = tags[conind]
		fwhm_vel = fwhm_vel[conind]
	endif else message,'No data found for '+ion+' in control files'
endif else if keyword_set(cso) then begin
	dir = '~/Astronomy/Research/Spitzer/cso/linedata/'
	restore,dir+ion+'.sav'
	tags = line.tag
	fwhm_vel = line.fwhm_vel
	fwhm_vel_err = line.fwhm_vel_err
endif else begin
	dir = '~/Astronomy/Research/Spitzer/linedata/'+skydir
	restore,dir+ion+'.sav'
	tags = line.tag
	fwhm_vel = line.fwhm_vel
	fwhm_vel_err = line.fwhm_vel_err
	ohmind = where(strmid(tags,0,4) eq 'mega' or strmid(tags,0,4) eq 'arch', ocount)

	if ocount gt 0 then begin
		tags = tags[ohmind]
		fwhm_vel = fwhm_vel[ohmind]
		fwhm_vel_err = fwhm_vel_err[ohmind]
	endif else message,'No data found for '+ion+' in megamaser files'

endelse


tagarr = transpose(tags)
sigma_arr = transpose([string(fwhm_vel / 2.35,format='(i7)')])
sigma_err_arr = transpose([string(fwhm_vel_err / 2.35,format='(i7)')])

; Test if the line is resolved - assume sigma \simeq 500 \pm 60 km/s

vel_res = 500. / 2.35
vel_res_err = 60. / 2.35

resolved_ind   = where(float(sigma_arr) - float(sigma_err_arr) gt float(vel_res) + float(vel_res_err), res_count)
unresolved_ind = where(float(sigma_arr) - float(sigma_err_arr) lt float(vel_res) + float(vel_res_err), unres_count)

linearr_res = [tagarr, sigma_arr]
linearr_unres = [tagarr, sigma_arr]

print,''
print,'Velocity resolution: '+string(vel_res+vel_res_err,format='(i5)')+' km/s'

if not keyword_set(quiet) then begin
	if res_count gt 0 then begin
		print,''
		print,'Resolved lines: '
		print,''
		print,'Object   sigma_i                              '+string(n_elements(linearr_res[0,*]))+' objects'
		print,'       [km/s]'
		print,''
		print, [transpose(tagarr[resolved_ind]),transpose(sigma_arr[resolved_ind])]
	endif else print,'No lines resolved for '+ion
	if unres_count gt 0 then begin
		print,''
		print,'Unresolved lines: '
		print,''
		print,'Object   sigma_i                              '+string(n_elements(linearr_res[0,*]))+' objects'
		print,'       [km/s]'
		print,''
		print, [transpose(tagarr[unresolved_ind]),transpose(sigma_arr[unresolved_ind])]
	endif else print, 'All lines resolved for '+ion
endif

if keyword_set(stop) then stop

end
