pro zspitzer, fname, all=all, stop = stop, verbose = verbose
;+
; NAME:
;       
;	ZSPITZER
;
; PURPOSE:
;
;	Measure the ensemble redshift of the high-res lines from Spitzer IRS data
;
; INPUTS:
;
;	FNAME - 	tag for object to measure (ex., 'arch001')
;
; OUTPUTS:
;
;
;
; KEYWORDS:
;
;	ALL - 		run ZSPITZER on all OHMs and non-masing objects
;
; EXAMPLE:
;
;	IDL> zspitzer, 'arch001'
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Apr 08
;	Added archived OHM data - Oct 08
; 	Added ALL keyword, revamped to work in velocity space		Apr 09
;-

c = 299792.458d

; Read in SMART data for object

if keyword_set(all) then begin
	otags = ohmdat('tag') & otags = otags[where(otags ne 'mega016')]
	atags = archdat('tag')
	ctags = condat('tag')
	fname = [otags,transpose(atags),transpose(ctags)]
endif else if n_elements(fname) eq 0 then begin
	message,'Prompt for file tag of object (ex., mega001)'
	stop
endif

for j = 0, n_elements(fname) - 1 do begin
	
	tag, fname[j], dirtag
	dir = '~/Astronomy/Research/Spitzer/'+dirtag+'/lines/hrsky/hires/'
	
	filelist = file_search(dir+'*'+fname[j]+'*txt')
	nfiles = n_elements(filelist)
	
	if nfiles le 0 or filelist[0] eq '' then begin
		dir = '~/Astronomy/Research/Spitzer/'+dirtag+'/lines/nosky/hires/'
		filelist = file_search(dir+'*'+fname[j]+'*txt')
		nfiles = n_elements(filelist)
		if nfiles le 0 or filelist[0] eq '' then begin
			message,'No HR line measurements found for object '+fname[j]
			stop
		endif
	endif
	
	ion = strarr(nfiles)

	realcenter = dblarr(nfiles)
	meascenter = dblarr(nfiles)
	expcenter = dblarr(nfiles)

	zir_arr = dblarr(nfiles)
	zir_errarr = dblarr(nfiles)

	velir_arr = dblarr(nfiles)
	velir_errarr = dblarr(nfiles)
	
	print,''
	print,fname[j]

	for i = 0, nfiles - 1 do begin
	
		p1 = strsplit(filelist(i),/ex,'/')
		p1 = strsplit(p1(n_elements(p1)-1),/ex,'.')
		p1 = strsplit(p1(0),/ex,'_')
		p1 = strsplit(p1(1),/ex)
		ionname = strupcase(strmid(p1,0,1))+strmid(p1,1,strlen(p1)-1)
	
		nlines = file_lines(filelist(i))
		emptyarr = strarr(nlines)
		openr,/get_lun, lun1, filelist(i)
		readf, lun1, emptyarr
		close, lun1 & free_lun, lun1
	
		meascen_loc = where(strmid(emptyarr, 0, 13) eq  ' Line Center:')
		meascen_string = strmid(emptyarr(meascen_loc),33,8)
		meascen_err_string = strmid(emptyarr(meascen_loc),46,10)
		meascen = double(meascen_string)
		meascen_err = double(meascen_err_string)
		
		restcen_loc = where(strmid(emptyarr, 0, 10) eq  'Wavelength')
		restcen_string = strsplit(emptyarr(restcen_loc+2),/ex)
		restcen = double(restcen_string[0])
	
		realcenter(i) = restcen
		meascenter(i) = meascen
		ion(i) = ionname
	
		; Calculate redshift and velocity necessary to fit line
	
		z_ir = meascen / restcen - 1d
		z_ir_err = (meascen_err / restcen - 1d) + 1d

		vel_ir = c * z_ir
		vel_ir_err = c * z_ir_err

		if keyword_set(verbose) then begin
			print,ionname
			print,'   Meas. center:   ',string(meascen,format='(f7.4)')+' +- '+string(meascen_err,format='(f6.4)')
			print,'   Line redshift:   ',string(z_ir,format='(f6.4)')+' +- '+string(z_ir_err,format='(f7.5)')
		endif
	
		zir_arr(i) = z_ir
		zir_errarr(i) = z_ir_err
		velir_arr[i] = vel_ir
		velir_errarr[i] = vel_ir_err

	endfor
	
	
	; Weighted average of VELOCITIES (changed from redshifts)
	
	wtdmean, velir_arr, velir_errarr, wavg_vel

	; Error on IR velocity is taken from the scatter around the mean, NOT the sum of 1/sum(1/sig^2) (which assumes
	;	that all lines are an estimator of the same quantity)

	if n_elements(velir_arr) gt 1 then wstddev_vel = stddev(velir_arr) else wstddev_vel = velir_arr[0]
	
	; Plot data 
	
	!p.multi=[0,1,1]
	
	red = fsc_color("Red")
	yellow = fsc_color("Yellow")
	defcolor = fsc_color("White")
	
	plot, velir_arr, /nodata, $
		xrange=[-2,n_elements(velir_arr)+4], /xstyle, $
		yrange=[min(velir_arr)-200,max(velir_arr)+200], /ystyle, $
		xtitle='Line', $
		ytitle='IR velocity [km/s]', $
		charsize = 1.5, $
		title='IR velocities for '+fname[j]
	
	oploterror,velir_arr,velir_errarr,psym=symcat(16)
	xyouts,indgen(n_elements(velir_arr)),velir_arr,'   '+ion
	hor,wavg_vel,linestyle=1
	
	if dirtag eq 'OHM' then begin
		tagarr = ohmdat('tag')
		tagind = where(tagarr eq fname[j])
		obj = ohmdat('obj') & obj = obj[tagind]
		czopt_arr = long(ohmdat('cz_opt')) & cz_opt = czopt_arr[tagind]
		czopt_err_arr = long(ohmdat('cz_opt_err')) & cz_opt_err = czopt_err_arr[tagind]
		czoh_arr = long(ohmdat('czoh')) & czoh = czoh_arr[tagind]
		czoh_err_arr = long(ohmdat('czoh_err')) & czoh_err = czoh_err_arr[tagind]
	endif else if dirtag eq 'archived' then begin
		tagarr = archdat('tag')
		tagind = where(tagarr eq fname[j])
		obj = archdat('obj') & obj = obj[tagind]
		czopt_arr = long(archdat('cz_opt')) & cz_opt = czopt_arr[tagind]
		czopt_err_arr = long(archdat('cz_opt_err')) & cz_opt_err = czopt_err_arr[tagind]
		czoh_arr = long(archdat('czoh')) & czoh = czoh_arr[tagind]			; Must add velocity of archived OHMs to archdat
		czoh_err_arr = long(archdat('czoh_err')) & czoh_err = czoh_err_arr[tagind]
	endif else if dirtag eq 'control' then begin
		tagarr = condat('tag')
		tagind = where(tagarr eq fname[j])
		obj = condat('obj') & obj = obj[tagind]
		czopt_arr = long(condat('cz_opt')) & cz_opt = czopt_arr[tagind]
		czopt_err_arr = long(condat('cz_opt_err')) & cz_opt_err = czopt_err_arr[tagind]
	endif
	
	hor, cz_opt, color=yellow
	
	xyouts, n_elements(velir_arr), wavg_vel - 3d-4, 'v!Iavg!N = '+string(wavg_vel,format='(f7.1)'),charsize=1.5
	xyouts, n_elements(velir_arr), cz_opt - 100, 'v!Iopt!N = ' +string(cz_opt,format='(f7.1)'),charsize=1.5, color = yellow

	if dirtag ne 'control' then begin
		hor, czoh,   color=red
		xyouts, n_elements(velir_arr), czoh + 100, 'v!I1667!N = '+string(czoh,format='(f7.1)'),charsize=1.5, color = red
		legend,/bottom,/left,['Weighted avg. of v_IR','Optical velocity', '1667 MHz velocity'],linestyle=[1,0,0], $
			color = [defcolor,yellow,red]
		delv_ohm = abs(wavg_vel - czoh)
	endif else begin
		legend,/bottom,/left,['Weighted avg. of v_IR','Optical velocity'],linestyle=[1,0], color = [defcolor,yellow]
	endelse
	
	
	; Compute change in velocity between optical and infrared redshifts
	
	delv_opt = abs(cz_opt - wavg_vel)
	
	; Print to screen
	
	print,''
	print,'Weighted average: ',string(wavg_vel,format='(f7.1)'),' +- ',string(wstddev_vel,format='(f7.1)')
	print,'Optical velocity: ',string(cz_opt,format='(f7.1)'),' +- ',string(cz_opt_err,format='(f7.1)')
	if dirtag ne 'control' then print,'OHM velocity:     ',string(czoh,format='(f7.1)'),' +- ',string(czoh_err,format='(f7.1)')
	
	print,''
	print,'del v/v = '+string(abs(wavg_vel - cz_opt)/cz_opt,format='(e10.3)')
	print,'IR-opt del v = '+string(delv_opt,format='(f8.2)')+' km/s'
	if dirtag ne 'control' then print,'IR-OHM del v = '+string(delv_ohm,format='(f8.2)')+' km/s'
	
	; Save results to IDL file
	
	if dirtag ne 'control' then $
		save,filename=dir+'saved/'+fname[j]+'_zlines.sav',$
		realcenter,expcenter,meascenter,ion,zir_arr,zir_errarr, velir_arr, velir_errarr, $
			wavg_vel,wstddev_vel,cz_opt,cz_opt_err, czoh, czoh_err, delv_opt, delv_ohm else $
		save,filename=dir+'saved/'+fname[j]+'_zlines.sav',$
		realcenter,expcenter,meascenter,ion,zir_arr,zir_errarr, velir_arr, velir_errarr, $
			wavg_vel,wstddev_vel,cz_opt,cz_opt_err, delv_opt, delv_ohm

endfor

if keyword_set(stop) then stop

end
