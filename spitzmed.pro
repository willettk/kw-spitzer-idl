pro spitzmed, specname, shsky=shsky, lhsky = lhsky, ss = ss, modules = modules, lores = lores, stop=stop
;+
; NAME:
;       SPITZMED
;
; PURPOSE:
;
; 	Median Spitzer BCD and droop images for all modules in order to remove cosmic rays and bad pixels 
;	and sky subtract
;
; INPUTS:
;
;	SPECNAME - string giving designation of target (eg, 'mega001')
;
;	MODULES - four element integer array designating which modules to run. 
;			1 = yes, 0 = no
;			ex. [LH, SH, LL, SL]
;
; OUTPUTS:
;
;	- Medianed FITS files from all cycles for each module and nod position, written to spitzmed directory
;
; KEYWORDS:
;
;	SHSKY - set if images have separate sky backgrounds for SH module
;
;	LHSKY - set if images have separate sky backgrounds for LH module
;
;	SS - set to subtract co-added 'superskies' for low-res modules and write to different directory
;
; EXAMPLE:
;	IDL> spitzmed, 'mega023'
;
; NOTES:
;
;	Note that the sky does not necessarily have the same number of exposures
;	as the target image (and usually doesn't, in the case of CSOs).
;
;	This program combines code from earlier routines SPITZMED_HR and SPITZMED_LR.
;
; REVISION HISTORY
;       Written by K. Willett                May 2007
;	Added directory switch for OHM/CSO - KW, May 07
;	Added NO_LHSKY keyword		- KW, Jun 2007
;	Added DROOPMED call		- KW, Jul 07
;	Added MODULES keyword to allow specific modules to be run - KW, Aug 07
;	Changed numbering scheme for reading LR files; previous version gave incorrect
;		files if the number of cycles was different from SL1 to SL2, for example - KW, Nov 07
; 	Changed default to no sky subtraction; added SHSKY, LHSKY keywords	- KW, Mar 08
;	Changed SL, LL to i mod 4 (multiple exposures at one nod position) - KW, Mar 08
;	Added loop for DROOPMED; now allows single cycles for MEDIAN function - KW, Jul 08
;	Changed 'bcd' filetype to 'bcd.' - new format has 'bcdb' and 'bcdr' files in same directory - KW, Jun 08
;-

; Run twice on BCD and droop files
	
for m=0,1 do begin

	if m eq 0 then filetype='bcd.' else filetype='droop'

	; Locate directory from which to read 
	
	tag, specname, dirtag
	
	; Specific modules to median (some spectra only have lores or hires, for example)
	
	if n_elements(modules) eq 0 then modules = [1,1,1,1]
	if keyword_set(lores) then modules = [0,0,1,1]
	
	do_lh = modules(0)
	do_sh = modules(1)
	do_ll = modules(2)
	do_sl = modules(3)
	
	
	; #########
	; ##  LH ##
	; #########
	
	if do_lh eq 1 then begin
	
	; Designate directory out of which to read BCD images (downloaded from Leopard)
	
	fpath_lh = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch3/bcd/'
	files_in_lh = fpath_lh + '*'+filetype+'*fits'
	
	files_lh = findfile(files_in_lh)
	count_lh = n_elements(files_lh)
	
	; Should be even number of files in directory (two nod positions)
	
	if count_lh mod 2 ne 0 then print,'Wrong number of files in directory for LH module'
	
	; Use first file to find dimensions of FITS images
	
	junk = readfits(files_lh(0), hdr, /silent)
	sj = size(junk)
	if sj(0) ne 2 then begin
		print,'Not a 2D FITS file'
	endif
	
	; Create empty arrays
	
	lh_1p     = fltarr(sj(1),sj(2),count_lh / 2)
	lh_2p     = fltarr(sj(1),sj(2),count_lh / 2)
	
	; Read in files
	
	for i = 0, count_lh/2 - 1 do begin
	
		lh_1p(*,*,i)     = readfits(files_lh(i),                hdr_lh_1p, /silent)
		lh_2p(*,*,i)     = readfits(files_lh(count_lh/2 + i),   hdr_lh_2p, /silent)
	
	endfor
	
	; Median the files
	
	lh_size = size(lh_1p)
	if lh_size[0] eq 3 then begin
		medarr_lh_1p = median(lh_1p,dim = 3, /even)
		medarr_lh_2p = median(lh_2p,dim = 3, /even)
	endif else begin
		if filetype eq 'bcd.' then print,'LH module has only 1 cycle'
		medarr_lh_1p = lh_1p
		medarr_lh_2p = lh_2p
	endelse
	
	; Repeat all steps for the background sky, if it exists
	
	if keyword_set(lhsky) then begin
		
		skypath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_bkgrd/ch3/bcd/'
		skyfiles_in_lh = skypath+'*'+filetype+'*fits'
		
		skyfiles = findfile(skyfiles_in_lh)
		skycount = n_elements(skyfiles)
		
		if skycount mod 2 ne 0 then print,'Wrong number of files in directory for LH module (no sky)'
		skyjunk = readfits(skyfiles(0), skyhdr, /silent)
		skysize = size(skyjunk)
		if skysize(0) ne 2 then begin
			print,'Not a 2D FITS file'
		endif
		lh_1p_sky = fltarr(skysize(1),skysize(2),skycount / 2)
		lh_2p_sky = fltarr(skysize(1),skysize(2),skycount / 2)
		
		for i = 0, skycount/2 - 1 do begin
		
			lh_1p_sky(*,*,i) = readfits(skyfiles(i),                hdr_lh_1p_sky, /silent)
			lh_2p_sky(*,*,i) = readfits(skyfiles(skycount/2 + i),   hdr_lh_2p_sky, /silent)
		
		endfor
	
		lhsky_size = size(lh_1p_sky)
		if lhsky_size[0] eq 3 then begin
			medarr_lh_1p_sky = median(lh_1p_sky, dim = 3, /even)
			medarr_lh_2p_sky = median(lh_2p_sky, dim = 3, /even)
		endif else begin
			if filetype eq 'bcd.' then print,'LH sky has only 1 cycle'
			medarr_lh_1p_sky = lh_1p_sky
			medarr_lh_2p_sky = lh_2p_sky
		endelse
	
		; Subtract sky from science targets
		
		sub_medarr_lh_1p = medarr_lh_1p - medarr_lh_1p_sky
		sub_medarr_lh_2p = medarr_lh_2p - medarr_lh_2p_sky
	
	endif else begin	; NOTE - this means that the 'sub' part of the ch3 filenames is a blatant lie if /LHSKY is not set,
				; but it's not worth changing this in all downstream programs. 
		sub_medarr_lh_1p = medarr_lh_1p
		sub_medarr_lh_2p = medarr_lh_2p
	
	endelse
	
	; Write to FITS files
	
	if filetype eq 'bcd.' then file_ext = 'fits' else file_ext='.fits'
	writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch3_sub_medarr_1p_'+filetype+file_ext, sub_medarr_lh_1p, hdr_lh_1p
	writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch3_sub_medarr_2p_'+filetype+file_ext, sub_medarr_lh_2p, hdr_lh_2p
	
	endif
	
	; #########
	; ##  SH ##
	; #########
	
	if do_sh eq 1 then begin
	
	; Designate directory out of which to read images
	
	fpath_sh = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch1/bcd/'
	files_in_sh = fpath_sh + '*'+filetype+'*fits'
	
	files_sh = findfile(files_in_sh)
	count_sh = n_elements(files_sh)
	
	; Should be even number of files in directory (two nod positions)
	
	if count_sh mod 2 ne 0 then print,'Wrong number of files in directory for SH module'
	
	; Use first file to find dimensions of FITS images
	
	junk = readfits(files_sh(0), hdr, /silent)
	sj = size(junk)
	if sj(0) ne 2 then begin
		print,'Not a 2D FITS file'
	endif
	
	; Create empty arrays
	
	sh_1p     = fltarr(sj(1),sj(2),count_sh / 2)
	sh_2p     = fltarr(sj(1),sj(2),count_sh / 2)
	
	; Read in files
	
	for i = 0, count_sh/2 - 1 do begin
	
		sh_1p(*,*,i)     = readfits(files_sh(i),                hdr_sh_1p, /silent)
		sh_2p(*,*,i)     = readfits(files_sh(count_sh/2 + i),   hdr_sh_2p, /silent)
	
	endfor
	
	; Median the files
	
	sh_size = size(sh_1p)
	if sh_size[0] eq 3 then begin
		medarr_sh_1p = median(sh_1p,dim = 3, /even)
		medarr_sh_2p = median(sh_2p,dim = 3, /even)
	endif else begin
		if filetype eq 'bcd.' then print,'SH module has only 1 cycle'
		medarr_sh_1p = sh_1p
		medarr_sh_2p = sh_1p
	endelse
	
	if keyword_set(shsky) then begin
		
		skypath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_bkgrd/ch1/bcd/'
		skyfiles_in_sh = skypath+'*'+filetype+'*fits'
		
		skyfiles = findfile(skyfiles_in_sh)
		skycount = n_elements(skyfiles)
		
		if skycount mod 2 ne 0 then print,'Wrong number of files in directory for SH module (no sky)'
		skyjunk = readfits(skyfiles(0), skyhdr, /silent)
		skysize = size(skyjunk)
		if skysize(0) ne 2 then begin
			print,'Not a 2D FITS file'
		endif
		sh_1p_sky = fltarr(skysize(1),skysize(2),skycount / 2)
		sh_2p_sky = fltarr(skysize(1),skysize(2),skycount / 2)
		
		for i = 0, skycount/2 - 1 do begin
		
			sh_1p_sky(*,*,i) = readfits(skyfiles(i),                hdr_sh_1p_sky, /silent)
			sh_2p_sky(*,*,i) = readfits(skyfiles(skycount/2 + i),   hdr_sh_2p_sky, /silent)
		
		endfor
	
		shsky_size = size(sh_1p_sky)
		if shsky_size[0] eq 3 then begin
			medarr_sh_1p_sky = median(sh_1p_sky, dim = 3, /even)
			medarr_sh_2p_sky = median(sh_2p_sky, dim = 3, /even)
		endif else begin
			if filetype eq 'bcd.' then print,'SH sky has only 1 cycle'
			medarr_sh_1p_sky = sh_1p_sky
			medarr_sh_2p_sky = sh_2p_sky
		endelse
	
		; Subtract sky from science targets
		
		sub_medarr_sh_1p = medarr_sh_1p - medarr_sh_1p_sky
		sub_medarr_sh_2p = medarr_sh_2p - medarr_sh_2p_sky
	
	endif else begin	; NOTE - this means that the 'sub' part of the ch3 filenames is not true if /SHSKY is not set,
				; but it's not worth changing this in all downstream programs. 
		sub_medarr_sh_1p = medarr_sh_1p
		sub_medarr_sh_2p = medarr_sh_2p
	
	endelse
	
	; Write to FITS files - again, note that having '_sub_' in the filename for HR modules
	;			does NOT indicate whether or not the target has had backgrounds subtracted, unfortunately. 
	
	if filetype eq 'bcd.' then file_ext = 'fits' else file_ext='.fits'
	writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch1_medarr_1p_'+filetype+file_ext, sub_medarr_sh_1p, hdr_sh_1p
	writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch1_medarr_2p_'+filetype+file_ext, sub_medarr_sh_2p, hdr_sh_2p
	
	endif
	
	; #########
	; ##  LL ##
	; #########
	
	if do_ll eq 1 then begin
	
	; Designate directory out of which to read images
	
	fpath_ll = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch2/bcd/'
	files_in_ll = fpath_ll + '*'+filetype+'*fits'
	
	files_ll=findfile(files_in_ll)
	count_ll = n_elements(files_ll)
	digits_ll=strarr(count_ll)
	
	; Use first file to find dimensions of FITS images
	
	junk = readfits(files_ll(0), hdr, /silent)
	sj = size(junk)
	if sj(0) ne 2 then begin
		print,'Not a 2D FITS file'
	endif
	
	; General procedure for arbitrary number of cycles
	
	for i = 0, count_ll-1 do begin
		splitup = strsplit(files_ll(i),'/',/extract)
		lastname = splitup(n_elements(splitup)-1)
		split_underscore = strsplit(lastname,'_',/extract)
		fourdigit = split_underscore(3)
		digits_ll(i) = fourdigit
	endfor
	
	diff_dig_ll = digits_ll(rem_dup(digits_ll))
	
	if n_elements(diff_dig_ll)  mod 4 ne 0 then print,'Wrong number of files in directory for LL module'
	
	for i = 0, n_elements(diff_dig_ll) - 1 do begin
		whichfiles_ll = where(digits_ll eq diff_dig_ll(i))
		case i mod 4 of
			0: ll2_1p = fltarr(sj(1),sj(2),n_elements(whichfiles_ll))
			1: ll2_2p = fltarr(sj(1),sj(2),n_elements(whichfiles_ll))
			2: ll1_1p = fltarr(sj(1),sj(2),n_elements(whichfiles_ll))
			3: ll1_2p = fltarr(sj(1),sj(2),n_elements(whichfiles_ll))
		endcase
		for j = 0, n_elements(whichfiles_ll) - 1 do begin
			case i mod 4 of
				0: ll2_1p(*,*,j) = readfits(files_ll(whichfiles_ll(j)), hdr_ll2_1p, /silent)
				1: ll2_2p(*,*,j) = readfits(files_ll(whichfiles_ll(j)), hdr_ll2_2p, /silent)
				2: ll1_1p(*,*,j) = readfits(files_ll(whichfiles_ll(j)), hdr_ll1_1p, /silent)
				3: ll1_2p(*,*,j) = readfits(files_ll(whichfiles_ll(j)), hdr_ll1_2p, /silent)
			endcase
		endfor
	endfor
	
	; Median the files
	
	ll_size = size(ll2_1p)
	if ll_size[0] eq 3 then begin
		medarr_ll2_1p = median(ll2_1p,dim = 3, /even)
		medarr_ll2_2p = median(ll2_2p,dim = 3, /even)
		medarr_ll1_1p = median(ll1_1p,dim = 3, /even)
		medarr_ll1_2p = median(ll1_2p,dim = 3, /even)
	endif else begin
		if filetype eq 'bcd.' then print,'LL module has only 1 cycle'
		medarr_ll2_1p = ll2_1p
		medarr_ll2_2p = ll2_2p
		medarr_ll1_1p = ll1_1p
		medarr_ll1_2p = ll1_2p
	endelse
	
	; Create quasi-superskies by coadding the empty modules and off-position nod
	
	ll1_1p_supersky = (medarr_ll2_1p + medarr_ll2_2p + medarr_ll1_2p) / 3d
	ll1_2p_supersky = (medarr_ll2_1p + medarr_ll2_2p + medarr_ll1_1p) / 3d
	ll2_1p_supersky = (medarr_ll1_1p + medarr_ll1_2p + medarr_ll2_2p) / 3d
	ll2_2p_supersky = (medarr_ll1_1p + medarr_ll1_2p + medarr_ll2_2p) / 3d
	
	; Subtract empty modules and write FITS files
	
	if keyword_set(ss) then begin
	
		sub_medarr_ll2_1p = medarr_ll2_1p - ll2_1p_supersky
		sub_medarr_ll2_2p = medarr_ll2_2p - ll2_2p_supersky
		sub_medarr_ll1_1p = medarr_ll1_1p - ll1_1p_supersky
		sub_medarr_ll1_2p = medarr_ll1_2p - ll1_2p_supersky
	
		if filetype eq 'bcd.' then file_ext = 'fits' else file_ext='.fits'
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch2_sub_medarr_2o_1p_'+filetype+file_ext, sub_medarr_ll2_1p, hdr_ll2_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch2_sub_medarr_2o_2p_'+filetype+file_ext, sub_medarr_ll2_2p, hdr_ll2_2p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch2_sub_medarr_1o_1p_'+filetype+file_ext, sub_medarr_ll1_1p, hdr_ll1_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch2_sub_medarr_1o_2p_'+filetype+file_ext, sub_medarr_ll1_2p, hdr_ll1_2p
	
	endif else begin
	
		sub_medarr_ll2_1p = medarr_ll2_1p - medarr_ll2_2p
		sub_medarr_ll2_2p = medarr_ll2_2p - medarr_ll2_1p
		sub_medarr_ll1_1p = medarr_ll1_1p - medarr_ll1_2p
		sub_medarr_ll1_2p = medarr_ll1_2p - medarr_ll1_1p
	
		if filetype eq 'bcd.' then file_ext = 'fits' else file_ext='.fits'
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch2_sub_medarr_2o_1p_'+filetype+file_ext, sub_medarr_ll2_1p, hdr_ll2_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch2_sub_medarr_2o_2p_'+filetype+file_ext, sub_medarr_ll2_2p, hdr_ll2_2p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch2_sub_medarr_1o_1p_'+filetype+file_ext, sub_medarr_ll1_1p, hdr_ll1_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch2_sub_medarr_1o_2p_'+filetype+file_ext, sub_medarr_ll1_2p, hdr_ll1_2p
	
	endelse
	
	endif
	
	; #########
	; ##  SL ##
	; #########
	
	if do_sl eq 1 then begin
	
	; Designate directory out of which to read images
	
	fpath_sl = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch0/bcd/'
	files_in_sl = fpath_sl + '*'+filetype+'*fits'
	
	files_sl=findfile(files_in_sl)
	if filetype eq 'bcd.' and n_elements(files_sl) gt 8 then $
		files_sl = files_sl(8:n_elements(files_sl)-1)		; Skip peakup files for BCD?
	count_sl = n_elements(files_sl)
	digits_sl=strarr(count_sl)
	
	; Use first file to find dimensions of FITS images
	
	junk = readfits(files_sl(0), hdr, /silent)
	sj = size(junk)
	if sj(0) ne 2 then begin
		print,'Not a 2D FITS file'
	endif
	
	; General procedure for arbitrary number of cycles
	
	for i = 0, count_sl-1 do begin
		splitup = strsplit(files_sl(i),'/',/extract)
		lastname = splitup(n_elements(splitup)-1)
		split_underscore = strsplit(lastname,'_',/extract)
		fourdigit = split_underscore(3)
		digits_sl(i) = fourdigit
	endfor
	
	diff_dig_sl = digits_sl(rem_dup(digits_sl))
	
	if n_elements(diff_dig_sl) mod 4 ne 0 then begin
		print,m
		print,'Wrong number of files in directory for SL module'
		print,diff_dig_sl
	endif
	
	for i = 0, n_elements(diff_dig_sl) - 1 do begin
		whichfiles_sl = where(digits_sl eq diff_dig_sl(i))
		case i mod 4 of
			0: sl2_1p = fltarr(sj(1),sj(2),n_elements(whichfiles_sl))
			1: sl2_2p = fltarr(sj(1),sj(2),n_elements(whichfiles_sl))
			2: sl1_1p = fltarr(sj(1),sj(2),n_elements(whichfiles_sl))
			3: sl1_2p = fltarr(sj(1),sj(2),n_elements(whichfiles_sl))
		endcase
		for j = 0, n_elements(whichfiles_sl) - 1 do begin
			case i mod 4 of
				0: sl2_1p[*,*,j] = readfits(files_sl(whichfiles_sl(j)), hdr_sl2_1p, /silent)
				1: sl2_2p[*,*,j] = readfits(files_sl(whichfiles_sl(j)), hdr_sl2_2p, /silent)
				2: sl1_1p[*,*,j] = readfits(files_sl(whichfiles_sl(j)), hdr_sl1_1p, /silent)
				3: sl1_2p[*,*,j] = readfits(files_sl(whichfiles_sl(j)), hdr_sl1_2p, /silent)
			endcase
		endfor
	endfor
	
	; Median the files (if multiple cycles exist)
	
	if count_sl / 4 gt 1 then begin
	
		medarr_sl2_1p = median(sl2_1p,dim = 3, /even)
		medarr_sl2_2p = median(sl2_2p,dim = 3, /even)
		medarr_sl1_1p = median(sl1_1p,dim = 3, /even)
		medarr_sl1_2p = median(sl1_2p,dim = 3, /even)
	
	endif else begin
	
		medarr_sl2_1p = sl2_1p
		medarr_sl2_2p = sl2_2p
		medarr_sl1_1p = sl1_1p
		medarr_sl1_2p = sl1_2p
	
		if filetype eq 'bcd.' then print,'SL module has only 1 cycle'
	endelse
	
	; Create quasi-superskies by coadding the empty modules and off-position nod
	
	sl1_1p_supersky = (medarr_sl2_1p + medarr_sl2_2p + medarr_sl1_2p) / 3d
	sl1_2p_supersky = (medarr_sl2_1p + medarr_sl2_2p + medarr_sl1_1p) / 3d
	sl2_1p_supersky = (medarr_sl1_1p + medarr_sl1_2p + medarr_sl2_2p) / 3d
	sl2_2p_supersky = (medarr_sl1_1p + medarr_sl1_2p + medarr_sl2_2p) / 3d
	
	; Subtract empty modules and write FITS files
	
	if keyword_set(ss) then begin
	
		sub_medarr_sl2_1p = medarr_sl2_1p - sl2_1p_supersky
		sub_medarr_sl2_2p = medarr_sl2_2p - sl2_2p_supersky
		sub_medarr_sl1_1p = medarr_sl1_1p - sl1_1p_supersky
		sub_medarr_sl1_2p = medarr_sl1_2p - sl1_2p_supersky
	
		if filetype eq 'bcd.' then file_ext = 'fits' else file_ext='.fits'
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch0_sub_medarr_2o_1p_'+filetype+file_ext, sub_medarr_sl2_1p, hdr_sl2_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch0_sub_medarr_2o_2p_'+filetype+file_ext, sub_medarr_sl2_2p, hdr_sl2_2p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch0_sub_medarr_1o_1p_'+filetype+file_ext, sub_medarr_sl1_1p, hdr_sl1_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/ss/'+specname+'_ch0_sub_medarr_1o_2p_'+filetype+file_ext, sub_medarr_sl1_2p, hdr_sl1_2p
	
	endif else begin
	
		sub_medarr_sl2_1p = medarr_sl2_1p - medarr_sl2_2p
		sub_medarr_sl2_2p = medarr_sl2_2p - medarr_sl2_1p
		sub_medarr_sl1_1p = medarr_sl1_1p - medarr_sl1_2p
		sub_medarr_sl1_2p = medarr_sl1_2p - medarr_sl1_1p
	
		if filetype eq 'bcd.' then file_ext = 'fits' else file_ext='.fits'
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch0_sub_medarr_2o_1p_'+filetype+file_ext, sub_medarr_sl2_1p, hdr_sl2_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch0_sub_medarr_2o_2p_'+filetype+file_ext, sub_medarr_sl2_2p, hdr_sl2_2p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch0_sub_medarr_1o_1p_'+filetype+file_ext, sub_medarr_sl1_1p, hdr_sl1_1p
		writefits,'~/Astronomy/Research/Spitzer/spitzmed/'+specname+'_ch0_sub_medarr_1o_2p_'+filetype+file_ext, sub_medarr_sl1_2p, hdr_sl1_2p
	
	endelse
	
	endif
	
endfor
; Median droop files

;droopmed, specname, shsky = shsky, lhsky = lhsky, ss = ss, modules = modules

print,'SPITZMED completed for ',specname

; End program

if keyword_set(stop) then stop
end

