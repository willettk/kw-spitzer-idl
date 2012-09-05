pro spiceprep, specname, ss = ss, modules = modules, lores = lores
;+
; NAME:
;       SPICEPREP
;
; PURPOSE:
; 	Copy pixel mask and uncertainty files from Spitzer data and rename them for
;	batch processing in SPICE
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
;	Copies mask and uncertainty files to spitzmed/optimal directory with medianed BCD files
;	from SPITZMED, and writes text batch files to same directory for use with SPICE. 
;
; KEYWORDS:
;
;	SS - if set, reads the supersky-subtracted files from the directory
;		~/Astronomy/Research/Spitzer/spitzmed/ss/ for the low-res modules
;		and writes to a separate directory
;
; EXAMPLE:
;	IDL> spiceprep, 'mega023'
;
; NOTES:
;
;	Note that the sky does not necessarily have the same number of exposures
;	as the target image (and usually doesn't, in the case of CSOs).
;
;	Should add errors from sky in quadrature for best statistics - otherwise, just
;	using the noise from the raw data files will likely be too low. 
;
; REVISION HISTORY
;       Written by K. Willett                May 2007
;	Added directory switch for OHM/CSO, modified creation of batch files - KW, May 07
;	Changed batch file to make four (need to run each module separately) - KW, Jun 07
;	Added copying for droop files (used for IRSCLEAN_MASK), added _CLEAN
;		extensions to new batch textfiles - KW, Jul 07
;	Added MODULES keyword		- KW, Aug 07
;	Added i mod 4; use the first mask and uncertainty file if multiple ones exist - KW, Mar 08
;	Replaced all of the spawned copy commands with FILE_COPY; much faster response. - KW, Jul 08
;-

; Locate directory from which to read (either megamasers or CSOs)

tag, specname, dirtag

; Specific modules to median (some spectra only have lores or hires, for example)

if n_elements(modules) eq 0 then modules = [1,1,1,1]
if keyword_set(lores) then modules = [0,0,1,1]

do_lh = modules(0)
do_sh = modules(1)
do_ll = modules(2)
do_sl = modules(3)

; Set directory paths for reading and writing

if keyword_set(ss) then begin
	wpath = '/Users/willettk/Astronomy/Research/Spitzer/spitzmed/optimal/ss/' 
	spitzmed_dir = '/Users/willettk/Astronomy/Research/Spitzer/spitzmed/ss/' 
endif else begin 
	wpath = '/Users/willettk/Astronomy/Research/Spitzer/spitzmed/optimal/'
	spitzmed_dir = '/Users/willettk/Astronomy/Research/Spitzer/spitzmed/'
endelse

cd, spitzmed_dir


; #########
; ##  LH ##
; #########

if do_lh eq 1 then begin

; Designate directory out of which to read images (downloaded from Leopard)

fpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch3/bcd/'

files_in_lh_bmask = fpath + '*bmask*fits'
files_in_lh_func  = fpath + '*func*fits'
files_in_lh_dmask = fpath + '*dmask*fits'
files_in_lh_drunc = fpath + '*drunc*fits'

files_lh_bmask=findfile(files_in_lh_bmask)
files_lh_func=findfile(files_in_lh_func)
files_lh_dmask=findfile(files_in_lh_dmask)
files_lh_drunc=findfile(files_in_lh_drunc)

count_lh_bmask = n_elements(files_lh_bmask)
count_lh_func = n_elements(files_lh_func)
count_lh_dmask = n_elements(files_lh_dmask)
count_lh_drunc = n_elements(files_lh_drunc)

; Should be even number of files in directory (two nod positions)

if count_lh_bmask mod 2 ne 0 then print,'Wrong number of files in directory for LH bmask'
if count_lh_func mod 2 ne 0 then print,'Wrong number of files in directory for LH func'
if count_lh_dmask mod 2 ne 0 then print,'Wrong number of files in directory for LH dmask'
if count_lh_drunc mod 2 ne 0 then print,'Wrong number of files in directory for LH drunc'

; Copy mask, uncertainty files

file_copy,/overwrite,files_lh_bmask(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_1p_bmask.fits'
file_copy,/overwrite,files_lh_bmask(count_lh_bmask/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_2p_bmask.fits'

file_copy,/overwrite,files_lh_func(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_1p_func.fits'
file_copy,/overwrite,files_lh_func(count_lh_func/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_2p_func.fits'

file_copy,/overwrite,files_lh_dmask(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_1p_dmask.fits'
file_copy,/overwrite,files_lh_dmask(count_lh_dmask/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_2p_dmask.fits'

file_copy,/overwrite,files_lh_drunc(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_1p_drunc.fits'
file_copy,/overwrite,files_lh_drunc(count_lh_drunc/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch3_sub_medarr_2p_drunc.fits'

; Copy the bcd medianed files from the spitzmed directory to the optimal extraction directory

cd,'~/Astronomy/Research/Spitzer/spitzmed/'
file_copy,/overwrite, '*'+specname+'*ch3*bcd*','optimal/'
file_copy,/overwrite, '*'+specname+'*ch3*droop*','optimal/'

endif

; #########
; ##  SH ##
; #########

if do_sh eq 1 then begin

; Designate directory out of which to read images

fpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch1/bcd/'
files_in_sh_bmask = fpath + '*bmask*fits'
files_in_sh_func = fpath + '*func*fits'
files_in_sh_dmask = fpath + '*dmask*fits'
files_in_sh_drunc = fpath + '*drunc*fits'

files_sh_bmask = findfile(files_in_sh_bmask)
files_sh_func = findfile(files_in_sh_func)
files_sh_dmask = findfile(files_in_sh_dmask)
files_sh_drunc = findfile(files_in_sh_drunc)

count_sh_bmask = n_elements(files_sh_bmask)
count_sh_func = n_elements(files_sh_func)
count_sh_dmask = n_elements(files_sh_dmask)
count_sh_drunc = n_elements(files_sh_drunc)

; Should be even number of files in directory (two nod positions)

if count_sh_bmask mod 2 ne 0 then print,'Wrong number of files in directory for SH bmask'
if count_sh_func mod 2 ne 0 then print,'Wrong number of files in directory for SH func'
if count_sh_dmask mod 2 ne 0 then print,'Wrong number of files in directory for SH dmask'
if count_sh_drunc mod 2 ne 0 then print,'Wrong number of files in directory for SH drunc'

file_copy,/overwrite,files_sh_bmask(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_1p_bmask.fits'
file_copy,/overwrite,files_sh_bmask(count_sh_bmask/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_2p_bmask.fits'

file_copy,/overwrite,files_sh_func(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_1p_func.fits'
file_copy,/overwrite,files_sh_func(count_sh_func/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_2p_func.fits'

file_copy,/overwrite,files_sh_dmask(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_1p_dmask.fits'
file_copy,/overwrite,files_sh_dmask(count_sh_dmask/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_2p_dmask.fits'

file_copy,/overwrite,files_sh_drunc(0),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_1p_drunc.fits'
file_copy,/overwrite,files_sh_drunc(count_sh_drunc/2 - 1),'~/Astronomy/Research/Spitzer/spitzmed/optimal/'+specname+'_ch1_medarr_2p_drunc.fits'

; Copy the bcd medianed files from the spitzmed directory to theoptimal extraction directory

cd,'~/Astronomy/Research/Spitzer/spitzmed/'
file_copy,/overwrite, '*'+specname+'*ch1*bcd*','optimal/'
file_copy,/overwrite, '*'+specname+'*ch1*droop*','optimal/'

endif

; #########
; ##  LL ##
; #########

if do_ll eq 1 then begin

; Designate directory out of which to read images

fpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch2/bcd/'

files_in_ll_bmask = fpath + '*bmask*fits'
files_in_ll_func = fpath + '*func*fits'
files_in_ll_dmask = fpath + '*dmask*fits'
files_in_ll_drunc = fpath + '*drunc*fits'

files_ll_bmask=findfile(files_in_ll_bmask)
files_ll_func=findfile(files_in_ll_func)
files_ll_dmask=findfile(files_in_ll_dmask)
files_ll_drunc=findfile(files_in_ll_drunc)

count_ll_bmask = n_elements(files_ll_bmask)
count_ll_func = n_elements(files_ll_func)
count_ll_dmask = n_elements(files_ll_dmask)
count_ll_drunc = n_elements(files_ll_drunc)

digits_ll_bmask=strarr(count_ll_bmask)
digits_ll_func=strarr(count_ll_func)
digits_ll_dmask=strarr(count_ll_dmask)
digits_ll_drunc=strarr(count_ll_drunc)

; General procedure for arbitrary number of cycles

for i = 0, count_ll_bmask-1 do begin
	splitup_bmask = strsplit(files_ll_bmask(i),'/',/extract)
	lastname_bmask = splitup_bmask(n_elements(splitup_bmask)-1)
	split_underscore_bmask = strsplit(lastname_bmask,'_',/extract)
	fourdigit_bmask = split_underscore_bmask(3)
	digits_ll_bmask(i) = fourdigit_bmask
endfor
for i = 0, count_ll_func-1 do begin
	splitup_func = strsplit(files_ll_func(i),'/',/extract)
	lastname_func = splitup_func(n_elements(splitup_func)-1)
	split_underscore_func = strsplit(lastname_func,'_',/extract)
	fourdigit_func = split_underscore_func(3)
	digits_ll_func(i) = fourdigit_func
endfor
for i = 0, count_ll_dmask-1 do begin
	splitup_dmask = strsplit(files_ll_dmask(i),'/',/extract)
	lastname_dmask = splitup_dmask(n_elements(splitup_dmask)-1)
	split_underscore_dmask = strsplit(lastname_dmask,'_',/extract)
	fourdigit_dmask = split_underscore_dmask(3)
	digits_ll_dmask(i) = fourdigit_dmask
endfor
for i = 0, count_ll_drunc-1 do begin
	splitup_drunc = strsplit(files_ll_drunc(i),'/',/extract)
	lastname_drunc = splitup_drunc(n_elements(splitup_drunc)-1)
	split_underscore_drunc = strsplit(lastname_drunc,'_',/extract)
	fourdigit_drunc = split_underscore_drunc(3)
	digits_ll_drunc(i) = fourdigit_drunc
endfor

diff_dig_ll_bmask = digits_ll_bmask(rem_dup(digits_ll_bmask))
file_index_bmask = intarr(n_elements(diff_dig_ll_bmask))
diff_dig_ll_func = digits_ll_func(rem_dup(digits_ll_func))
file_index_func = intarr(n_elements(diff_dig_ll_func))
diff_dig_ll_dmask = digits_ll_dmask(rem_dup(digits_ll_dmask))
file_index_dmask = intarr(n_elements(diff_dig_ll_dmask))
diff_dig_ll_drunc = digits_ll_drunc(rem_dup(digits_ll_drunc))
file_index_drunc = intarr(n_elements(diff_dig_ll_drunc))

for i = 0, n_elements(diff_dig_ll_bmask) - 1 do begin
	whichfiles_ll_bmask = where(digits_ll_bmask eq diff_dig_ll_bmask(i))
	file_index_bmask(i) = whichfiles_ll_bmask(0)
endfor
for i = 0, n_elements(diff_dig_ll_func) - 1 do begin
	whichfiles_ll_func = where(digits_ll_func eq diff_dig_ll_func(i))
	file_index_func(i) = whichfiles_ll_func(0)
endfor
for i = 0, n_elements(diff_dig_ll_dmask) - 1 do begin
	whichfiles_ll_dmask = where(digits_ll_dmask eq diff_dig_ll_dmask(i))
	file_index_dmask(i) = whichfiles_ll_dmask(0)
endfor
for i = 0, n_elements(diff_dig_ll_drunc) - 1 do begin
	whichfiles_ll_drunc = where(digits_ll_drunc eq diff_dig_ll_drunc(i))
	file_index_drunc(i) = whichfiles_ll_drunc(0)
endfor

; General procedure for arbitrary number of cycles

if n_elements(file_index_bmask) mod 4 ne 0 then print,'Wrong number of files in directory for LL bmask'
if n_elements(file_index_func)  mod 4 ne 0 then print,'Wrong number of files in directory for LL func'
if n_elements(file_index_dmask) mod 4 ne 0 then print,'Wrong number of files in directory for LL dmask'
if n_elements(file_index_drunc) mod 4 ne 0 then print,'Wrong number of files in directory for LL drunc'

file_copy,/overwrite,files_ll_bmask(file_index_bmask(0)),wpath+specname+'_ch2_sub_medarr_2o_1p_bmask.fits'
file_copy,/overwrite,files_ll_bmask(file_index_bmask(1)),wpath+specname+'_ch2_sub_medarr_2o_2p_bmask.fits'
file_copy,/overwrite,files_ll_bmask(file_index_bmask(2)),wpath+specname+'_ch2_sub_medarr_1o_1p_bmask.fits'
file_copy,/overwrite,files_ll_bmask(file_index_bmask(3)),wpath+specname+'_ch2_sub_medarr_1o_2p_bmask.fits'

file_copy,/overwrite,files_ll_func(file_index_func(0)),wpath+specname+'_ch2_sub_medarr_2o_1p_func.fits'
file_copy,/overwrite,files_ll_func(file_index_func(1)),wpath+specname+'_ch2_sub_medarr_2o_2p_func.fits'
file_copy,/overwrite,files_ll_func(file_index_func(2)),wpath+specname+'_ch2_sub_medarr_1o_1p_func.fits'
file_copy,/overwrite,files_ll_func(file_index_func(3)),wpath+specname+'_ch2_sub_medarr_1o_2p_func.fits'

file_copy,/overwrite,files_ll_dmask(file_index_dmask(0)),wpath+specname+'_ch2_sub_medarr_2o_1p_dmask.fits'
file_copy,/overwrite,files_ll_dmask(file_index_dmask(1)),wpath+specname+'_ch2_sub_medarr_2o_2p_dmask.fits'
file_copy,/overwrite,files_ll_dmask(file_index_dmask(2)),wpath+specname+'_ch2_sub_medarr_1o_1p_dmask.fits'
file_copy,/overwrite,files_ll_dmask(file_index_dmask(3)),wpath+specname+'_ch2_sub_medarr_1o_2p_dmask.fits'

file_copy,/overwrite,files_ll_drunc(file_index_drunc(0)),wpath+specname+'_ch2_sub_medarr_2o_1p_drunc.fits'
file_copy,/overwrite,files_ll_drunc(file_index_drunc(1)),wpath+specname+'_ch2_sub_medarr_2o_2p_drunc.fits'
file_copy,/overwrite,files_ll_drunc(file_index_drunc(2)),wpath+specname+'_ch2_sub_medarr_1o_1p_drunc.fits'
file_copy,/overwrite,files_ll_drunc(file_index_drunc(3)),wpath+specname+'_ch2_sub_medarr_1o_2p_drunc.fits'

; Copy the bcd medianed files from the spitzmed directory to theoptimal extraction directory

file_copy,/overwrite, '*'+specname+'*ch2*bcd*','optimal/'
file_copy,/overwrite, '*'+specname+'*ch2*droop*','optimal/'

endif

; #########
; ##  SL ##
; #########

if do_sl eq 1 then begin

; Designate directory out of which to read images

fpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/'+specname+'_spectra/ch0/bcd/'

files_in_sl_bmask = fpath + '*bmask*fits'
files_in_sl_func = fpath + '*func*fits'
files_in_sl_dmask = fpath + '*dmask*fits'
files_in_sl_drunc = fpath + '*drunc*fits'

files_sl_bmask=findfile(files_in_sl_bmask)
files_sl_func=findfile(files_in_sl_func)
files_sl_dmask=findfile(files_in_sl_dmask)		; DMASK files are associated with peakups and need to be skipped
files_sl_drunc=findfile(files_in_sl_drunc)

count_sl_bmask = n_elements(files_sl_bmask)
count_sl_func = n_elements(files_sl_func)
count_sl_dmask = n_elements(files_sl_dmask)
count_sl_drunc = n_elements(files_sl_drunc)

digits_sl_bmask=strarr(count_sl_bmask)
digits_sl_func=strarr(count_sl_func)
digits_sl_dmask=strarr(count_sl_dmask)
digits_sl_drunc=strarr(count_sl_drunc)

; General procedure for arbitrary number of cycles

for i = 0, count_sl_bmask-1 do begin
	splitup_bmask = strsplit(files_sl_bmask(i),'/',/extract)
	lastname_bmask = splitup_bmask(n_elements(splitup_bmask)-1)
	split_underscore_bmask = strsplit(lastname_bmask,'_',/extract)
	fourdigit_bmask = split_underscore_bmask(3)
	digits_sl_bmask(i) = fourdigit_bmask
endfor
for i = 0, count_sl_func-1 do begin
	splitup_func = strsplit(files_sl_func(i),'/',/extract)
	lastname_func = splitup_func(n_elements(splitup_func)-1)
	split_underscore_func = strsplit(lastname_func,'_',/extract)
	fourdigit_func = split_underscore_func(3)
	digits_sl_func(i) = fourdigit_func
endfor
for i = 0, count_sl_dmask-1 do begin
	splitup_dmask = strsplit(files_sl_dmask(i),'/',/extract)
	lastname_dmask = splitup_dmask(n_elements(splitup_dmask)-1)
	split_underscore_dmask = strsplit(lastname_dmask,'_',/extract)
	fourdigit_dmask = split_underscore_dmask(3)
	digits_sl_dmask(i) = fourdigit_dmask
endfor
for i = 0, count_sl_drunc-1 do begin
	splitup_drunc = strsplit(files_sl_drunc(i),'/',/extract)
	lastname_drunc = splitup_drunc(n_elements(splitup_drunc)-1)
	split_underscore_drunc = strsplit(lastname_drunc,'_',/extract)
	fourdigit_drunc = split_underscore_drunc(3)
	digits_sl_drunc(i) = fourdigit_drunc
endfor

diff_dig_sl_bmask = digits_sl_bmask(rem_dup(digits_sl_bmask))
file_index_bmask = intarr(n_elements(diff_dig_sl_bmask))
diff_dig_sl_func = digits_sl_func(rem_dup(digits_sl_func))
file_index_func = intarr(n_elements(diff_dig_sl_func))
diff_dig_sl_dmask = digits_sl_dmask(rem_dup(digits_sl_dmask))
diff_dig_sl_dmask = diff_dig_sl_dmask(2:n_elements(diff_dig_sl_dmask)-1)
file_index_dmask = intarr(n_elements(diff_dig_sl_dmask))
diff_dig_sl_drunc = digits_sl_drunc(rem_dup(digits_sl_drunc))
file_index_drunc = intarr(n_elements(diff_dig_sl_drunc))

for i = 0, n_elements(diff_dig_sl_bmask) - 1 do begin
	whichfiles_sl_bmask = where(digits_sl_bmask eq diff_dig_sl_bmask(i))
	file_index_bmask(i) = whichfiles_sl_bmask(0)
endfor
for i = 0, n_elements(diff_dig_sl_func) - 1 do begin
	whichfiles_sl_func = where(digits_sl_func eq diff_dig_sl_func(i))
	file_index_func(i) = whichfiles_sl_func(0)
endfor
for i = 0, n_elements(diff_dig_sl_dmask) - 1 do begin
	whichfiles_sl_dmask = where(digits_sl_dmask eq diff_dig_sl_dmask(i))
	file_index_dmask(i) = whichfiles_sl_dmask(0)
endfor
for i = 0, n_elements(diff_dig_sl_drunc) - 1 do begin
	whichfiles_sl_drunc = where(digits_sl_drunc eq diff_dig_sl_drunc(i))
	file_index_drunc(i) = whichfiles_sl_drunc(0)
endfor

; General procedure for arbitrary number of cycles

if n_elements(file_index_bmask) mod 4 ne 0 then print,'Wrong number of files in directory for SL bmask'
if n_elements(file_index_func)  mod 4 ne 0 then print,'Wrong number of files in directory for SL func'
if n_elements(file_index_dmask) mod 4 ne 0 then print,'Wrong number of files in directory for SL dmask' 
if n_elements(file_index_drunc) mod 4 ne 0 then print,'Wrong number of files in directory for SL drunc'

file_copy,/overwrite,files_sl_bmask(file_index_bmask(0)),wpath+specname+'_ch0_sub_medarr_2o_1p_bmask.fits'
file_copy,/overwrite,files_sl_bmask(file_index_bmask(1)),wpath+specname+'_ch0_sub_medarr_2o_2p_bmask.fits'
file_copy,/overwrite,files_sl_bmask(file_index_bmask(2)),wpath+specname+'_ch0_sub_medarr_1o_1p_bmask.fits'
file_copy,/overwrite,files_sl_bmask(file_index_bmask(3)),wpath+specname+'_ch0_sub_medarr_1o_2p_bmask.fits'

file_copy,/overwrite,files_sl_func(file_index_func(0)),wpath+specname+'_ch0_sub_medarr_2o_1p_func.fits'
file_copy,/overwrite,files_sl_func(file_index_func(1)),wpath+specname+'_ch0_sub_medarr_2o_2p_func.fits'
file_copy,/overwrite,files_sl_func(file_index_func(2)),wpath+specname+'_ch0_sub_medarr_1o_1p_func.fits'
file_copy,/overwrite,files_sl_func(file_index_func(3)),wpath+specname+'_ch0_sub_medarr_1o_2p_func.fits'

file_copy,/overwrite,files_sl_dmask(file_index_dmask(0)),wpath+specname+'_ch0_sub_medarr_2o_1p_dmask.fits'
file_copy,/overwrite,files_sl_dmask(file_index_dmask(1)),wpath+specname+'_ch0_sub_medarr_2o_2p_dmask.fits'
file_copy,/overwrite,files_sl_dmask(file_index_dmask(2)),wpath+specname+'_ch0_sub_medarr_1o_1p_dmask.fits'
file_copy,/overwrite,files_sl_dmask(file_index_dmask(3)),wpath+specname+'_ch0_sub_medarr_1o_2p_dmask.fits'

file_copy,/overwrite,files_sl_drunc(file_index_drunc(0)),wpath+specname+'_ch0_sub_medarr_2o_1p_drunc.fits'
file_copy,/overwrite,files_sl_drunc(file_index_drunc(1)),wpath+specname+'_ch0_sub_medarr_2o_2p_drunc.fits'
file_copy,/overwrite,files_sl_drunc(file_index_drunc(2)),wpath+specname+'_ch0_sub_medarr_1o_1p_drunc.fits'
file_copy,/overwrite,files_sl_drunc(file_index_drunc(3)),wpath+specname+'_ch0_sub_medarr_1o_2p_drunc.fits'

; Copy the bcd medianed files from the spitzmed directory to theoptimal extraction directory

file_copy,/overwrite, '*'+specname+'*ch0*bcd*','optimal/'
file_copy,/overwrite, '*'+specname+'*ch0*droop*','optimal/'

endif

;
; ############################################
;

; Create list of files to use in SPICE batch processing - this list still has odd characters in it
; from the ls UNIX command, so it needs to be manually edited before it goes into SPICE

lopath = wpath
hipath = '/Users/willettk/Astronomy/Research/Spitzer/spitzmed/optimal/'

batch_sl = ['_ch0_sub_medarr_2o_1p_bcd_clean.fits', $
		'_ch0_sub_medarr_2o_2p_bcd_clean.fits', $
		'_ch0_sub_medarr_1o_1p_bcd_clean.fits', $
		'_ch0_sub_medarr_1o_2p_bcd_clean.fits', $
		'_ch0_sub_medarr_2o_1p_bmask_clean.fits', $
		'_ch0_sub_medarr_2o_2p_bmask_clean.fits', $
		'_ch0_sub_medarr_1o_1p_bmask_clean.fits', $
		'_ch0_sub_medarr_1o_2p_bmask_clean.fits', $
		'_ch0_sub_medarr_2o_1p_func_clean.fits', $
		'_ch0_sub_medarr_2o_2p_func_clean.fits', $
		'_ch0_sub_medarr_1o_1p_func_clean.fits', $
		'_ch0_sub_medarr_1o_2p_func_clean.fits']

batch_ll = [ '_ch2_sub_medarr_2o_1p_bcd_clean.fits', $
		'_ch2_sub_medarr_2o_2p_bcd_clean.fits', $
		'_ch2_sub_medarr_1o_1p_bcd_clean.fits', $
		'_ch2_sub_medarr_1o_2p_bcd_clean.fits', $
		'_ch2_sub_medarr_2o_1p_bmask_clean.fits', $
		'_ch2_sub_medarr_2o_2p_bmask_clean.fits', $
		'_ch2_sub_medarr_1o_1p_bmask_clean.fits', $
		'_ch2_sub_medarr_1o_2p_bmask_clean.fits', $
		'_ch2_sub_medarr_2o_1p_func_clean.fits', $
		'_ch2_sub_medarr_2o_2p_func_clean.fits', $
		'_ch2_sub_medarr_1o_1p_func_clean.fits', $
		'_ch2_sub_medarr_1o_2p_func_clean.fits']

batch_sh = ['_ch1_medarr_1p_bcd_clean.fits', $
		'_ch1_medarr_2p_bcd_clean.fits', $
		'_ch1_medarr_1p_bmask_clean.fits', $
		'_ch1_medarr_2p_bmask_clean.fits', $
		'_ch1_medarr_1p_func_clean.fits', $
		'_ch1_medarr_2p_func_clean.fits']

batch_lh = ['_ch3_sub_medarr_1p_bcd_clean.fits', $
		'_ch3_sub_medarr_2p_bcd_clean.fits', $
		'_ch3_sub_medarr_1p_bmask_clean.fits', $
		'_ch3_sub_medarr_2p_bmask_clean.fits', $
		'_ch3_sub_medarr_1p_func_clean.fits', $
		'_ch3_sub_medarr_2p_func_clean.fits']

batch_sl = lopath + specname + batch_sl
batch_ll = lopath + specname + batch_ll
batch_sh = hipath + specname + batch_sh
batch_lh = hipath + specname + batch_lh
		
if do_sl eq 1 then forprint, textout = lopath+'batch_'+specname+'_sl.txt', batch_sl, /nocomment, /silent
if do_ll eq 1 then forprint, textout = lopath+'batch_'+specname+'_ll.txt', batch_ll, /nocomment, /silent
if do_sh eq 1 then forprint, textout = hipath+'batch_'+specname+'_sh.txt', batch_sh, /nocomment, /silent
if do_lh eq 1 then forprint, textout = hipath+'batch_'+specname+'_lh.txt', batch_lh, /nocomment, /silent

print,'SPICEPREP completed for ',specname

; End program

;stop
end

