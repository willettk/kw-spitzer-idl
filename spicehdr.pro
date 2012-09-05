pro spicehdr, fname, readdir = readdir, lores = lores, no_clean = no_clean
;+
; NAME:
;       SPICEHDR
;
; PURPOSE:
; 	Extract SPICE headers and add them to the coadded and trimmed spectra as output by SPEC2PAHFIT and FULLQL_HR. 
;
; INPUTS:
;
; 	FNAME - string describing object in Spitzer nomenclature (eg, 'mega005')
;
;	READDIR - string giving directory in which to place the files; 
;		typically set to either 'stitched' (when called by SMARTMAKE) or 'calibrated' (called by SPECCAL; default)
;
; OUTPUTS:
;
;	Writes SPICE-like file to data directory with header; suitable for importing into SMART for
;	measurements of line fluxes.
;
; KEYWORDS:
;
;	SPECCAL - use the PU-calibrated and stitched files from the second run through.
;
;	NO_CLEAN - operates on data not cleaned by IRSCLEAN_MASK
;
; EXAMPLE:
;	IDL> spicehdr, 'mega005', readdir = 'stitched'
;
; REQUIRES:
;
;	TAG.pro
;
; NOTES:
;
; REVISION HISTORY
;       Written by K. Willett                Jun 2007
;	Add individual nods as well as coadded spectra - KW, Jun 07
;	Changed files to peakup-calibrated versions written by SPEC2PUFLUX.pro - KW, Aug 07
;	Added bonus order - KW, Aug 07
;	Modified program to be called by either SMARTMAKE or SPECCAL - KW, Feb 08
;	Updated keywords and documentation - KW, May 08
;	Changed location of writepath - KW, Jul 08
;	Was missing the LH_1p write command, for some reason - Aug 08
;-

; Close all units to free up memory

close, /all

; Define directory paths used

tag, fname, dirtag		; CSO vs. OHM

	hdrpath_lo = '~/spice/output/lores/optimal/all/'
	hdrpath_hi = '~/spice/output/hires/optimal/all/'
	trimpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/'+readdir+'/coadd/'
	nodpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/'+readdir+'/nods/'
	writepath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/smartfiles/coadd/'
	nwritepath = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/smartfiles/nods/'

if n_elements(readdir) eq 0 then begin
	message,'Must specify readdir (eg, stitched or calibrated)'
	stop
endif

if not keyword_set(lores) then begin

; ######
; # SH #
; ######

; Open spectrum from SPICE

if keyword_set(no_clean) then file = hdrpath_hi+fname+'_ch1_medarr_1p_spect.tbl' $
	else file = hdrpath_hi+fname+'_ch1_medarr_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

; Read in entire file

openr, sh_lun, file, /get_lun
readf, sh_lun, emptyarr
close, sh_lun & free_lun, sh_lun

; Isolate only the SPICE header

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
sh_arr = temporary(emptyarr(0:eoh))

; Read in the data from the trimmed and coadded spectrum

shco_file = '_sh_cal.tbl'
sh1p_file = '_sh_1p_cal.tbl'
sh2p_file = '_sh_2p_cal.tbl'


readcol, trimpath+fname+shco_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sh1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sh2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent

; Print the header and spectrum to separate files (easiest to keep formatting this way)

forprint, sh_arr, format = 'a', textout = writepath+fname+'_ch1.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch1.dat', /nocomment, /silent

forprint, sh_arr, format = 'a', textout = nwritepath+fname+'_ch1.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch1_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch1_2p.dat', /nocomment, /silent

; Combine the header and data, remove temp. files

cd, writepath

spawn, 'cat '+fname+'_ch1.hdr'+' '+fname+'_ch1.dat'+' > '+fname+'_sh_hdr.tbl'
file_delete,fname+'_ch1.hdr'
file_delete,fname+'_ch1.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch1.hdr'+' '+fname+'_ch1_1p.dat'+' > '+fname+'_sh_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch1.hdr'+' '+fname+'_ch1_2p.dat'+' > '+fname+'_sh_2p_hdr.tbl'
file_delete,fname+'_ch1.hdr'
file_delete,fname+'_ch1_1p.dat'
file_delete,fname+'_ch1_2p.dat'

;
; All other modules follow identical procedure, with small changes due to filenaming scheme
;

; ######
; # LH #
; ######

if keyword_set(no_clean) then file = hdrpath_hi+fname+'_ch3_sub_medarr_1p_spect.tbl' $
	else file = hdrpath_hi+fname+'_ch3_sub_medarr_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lh_lun, file, /get_lun
readf, lh_lun, emptyarr
close, lh_lun & free_lun, lh_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
lh_arr = temporary(emptyarr(0:eoh))

; Read in the data from the trimmed and coadded spectrum

lhco_file = '_lh_cal.tbl'
lh1p_file = '_lh_1p_cal.tbl'
lh2p_file = '_lh_2p_cal.tbl'


readcol, trimpath+fname+lhco_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+lh1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+lh2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent


forprint, lh_arr, format = 'a', textout = writepath+fname+'_ch3.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch3.dat', /nocomment, /silent

forprint, lh_arr, format = 'a', textout = nwritepath+fname+'_ch3.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch3_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch3_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch3.hdr'+' '+fname+'_ch3.dat'+' > '+fname+'_lh_hdr.tbl'
file_delete,fname+'_ch3.hdr'
file_delete,fname+'_ch3.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch3.hdr'+' '+fname+'_ch3_1p.dat'+' > '+fname+'_lh_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch3.hdr'+' '+fname+'_ch3_2p.dat'+' > '+fname+'_lh_2p_hdr.tbl'
file_delete,fname+'_ch3.hdr'
file_delete,fname+'_ch3_1p.dat'
file_delete,fname+'_ch3_2p.dat'

; End LORES keyword

endif

; #######
; # SL1 #
; #######

if keyword_set(no_clean) then file = hdrpath_lo+fname+'_ch0_sub_medarr_1o_1p_spect.tbl' $
	else file = hdrpath_lo+fname+'_ch0_sub_medarr_1o_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, sl1_lun, file, /get_lun
readf, sl1_lun, emptyarr
close, sl1_lun & free_lun, sl1_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
sl1_arr = temporary(emptyarr(0:eoh))

sl1_co_file = '_sl1_cal.tbl'
sl1_1p_file = '_sl1_1p_cal.tbl'
sl1_2p_file = '_sl1_2p_cal.tbl'

readcol, trimpath+fname+sl1_co_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sl1_1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sl1_2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent

forprint, sl1_arr, format = 'a', textout = writepath+fname+'_ch0.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch0.dat', /nocomment, /silent

forprint, sl1_arr, format = 'a', textout = nwritepath+fname+'_ch0.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch0_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch0_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0.dat'+' > '+fname+'_sl1_hdr.tbl'
file_delete,fname+'_ch0.hdr'
file_delete,fname+'_ch0.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0_1p.dat'+' > '+fname+'_sl1_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0_2p.dat'+' > '+fname+'_sl1_2p_hdr.tbl'
file_delete,fname+'_ch0.hdr'
file_delete,fname+'_ch0_1p.dat'
file_delete,fname+'_ch0_2p.dat'


; #######
; # SL3 #
; #######

if keyword_set(no_clean) then file = hdrpath_lo+fname+'_ch0_sub_medarr_2o_1p_spect.tbl' $
	else file = hdrpath_lo+fname+'_ch0_sub_medarr_2o_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, sl3_lun, file, /get_lun
readf, sl3_lun, emptyarr
close, sl3_lun & free_lun, sl3_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
sl3_arr = temporary(emptyarr(0:eoh))

sl3_co_file = '_sl3_cal.tbl'
sl3_1p_file = '_sl3_1p_cal.tbl'
sl3_2p_file = '_sl3_2p_cal.tbl'

readcol, trimpath+fname+sl3_co_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sl3_1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sl3_2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent


forprint, sl3_arr, format = 'a', textout = writepath+fname+'_ch0.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch0.dat', /nocomment, /silent

forprint, sl3_arr, format = 'a', textout = nwritepath+fname+'_ch0.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch0_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch0_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0.dat'+' > '+fname+'_sl3_hdr.tbl'
file_delete,fname+'_ch0.hdr'
file_delete,fname+'_ch0.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0_1p.dat'+' > '+fname+'_sl3_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0_2p.dat'+' > '+fname+'_sl3_2p_hdr.tbl'
file_delete,fname+'_ch0.hdr'
file_delete,fname+'_ch0_1p.dat'
file_delete,fname+'_ch0_2p.dat'

; #######
; # SL2 #
; #######

if keyword_set(no_clean) then file = hdrpath_lo+fname+'_ch0_sub_medarr_2o_1p_spect.tbl' $
	else file = hdrpath_lo+fname+'_ch0_sub_medarr_2o_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, sl2_lun, file, /get_lun
readf, sl2_lun, emptyarr
close, sl2_lun & free_lun, sl2_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
sl2_arr = temporary(emptyarr(0:eoh))

sl2_co_file = '_sl2_cal.tbl'
sl2_1p_file = '_sl2_1p_cal.tbl'
sl2_2p_file = '_sl2_2p_cal.tbl'

readcol, trimpath+fname+sl2_co_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sl2_1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+sl2_2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent


forprint, sl2_arr, format = 'a', textout = writepath+fname+'_ch0.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch0.dat', /nocomment, /silent

forprint, sl2_arr, format = 'a', textout = nwritepath+fname+'_ch0.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch0_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch0_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0.dat'+' > '+fname+'_sl2_hdr.tbl'
file_delete,fname+'_ch0.hdr'
file_delete,fname+'_ch0.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0_1p.dat'+' > '+fname+'_sl2_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch0.hdr'+' '+fname+'_ch0_2p.dat'+' > '+fname+'_sl2_2p_hdr.tbl'
file_delete,fname+'_ch0.hdr'
file_delete,fname+'_ch0_1p.dat'
file_delete,fname+'_ch0_2p.dat'

; #######
; # LL1 #
; #######

if keyword_set(no_clean) then file = hdrpath_lo+fname+'_ch2_sub_medarr_1o_1p_spect.tbl' $
	else file = hdrpath_lo+fname+'_ch2_sub_medarr_1o_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, ll1_lun, file, /get_lun
readf, ll1_lun, emptyarr
close, ll1_lun & free_lun, ll1_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
ll1_arr = temporary(emptyarr(0:eoh))

ll1_co_file = '_ll1_cal.tbl'
ll1_1p_file = '_ll1_1p_cal.tbl'
ll1_2p_file = '_ll1_2p_cal.tbl'

readcol, trimpath+fname+ll1_co_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+ll1_1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+ll1_2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent

forprint, ll1_arr, format = 'a', textout = writepath+fname+'_ch2.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch2.dat', /nocomment, /silent

forprint, ll1_arr, format = 'a', textout = nwritepath+fname+'_ch2.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch2_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch2_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2.dat'+' > '+fname+'_ll1_hdr.tbl'
file_delete,fname+'_ch2.hdr'
file_delete,fname+'_ch2.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2_1p.dat'+' > '+fname+'_ll1_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2_2p.dat'+' > '+fname+'_ll1_2p_hdr.tbl'
file_delete,fname+'_ch2.hdr'
file_delete,fname+'_ch2_1p.dat'
file_delete,fname+'_ch2_2p.dat'


; #######
; # LL3 #
; #######

if keyword_set(no_clean) then file = hdrpath_lo+fname+'_ch2_sub_medarr_2o_1p_spect.tbl' $
	else file = hdrpath_lo+fname+'_ch2_sub_medarr_2o_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, ll3_lun, file, /get_lun
readf, ll3_lun, emptyarr
close, ll3_lun & free_lun, ll3_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
ll3_arr = temporary(emptyarr(0:eoh))

ll3_co_file = '_ll3_cal.tbl'
ll3_1p_file = '_ll3_1p_cal.tbl'
ll3_2p_file = '_ll3_2p_cal.tbl'

readcol, trimpath+fname+ll3_co_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+ll3_1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+ll3_2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent


forprint, ll3_arr, format = 'a', textout = writepath+fname+'_ch2.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch2.dat', /nocomment, /silent

forprint, ll3_arr, format = 'a', textout = nwritepath+fname+'_ch2.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch2_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch2_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2.dat'+' > '+fname+'_ll3_hdr.tbl'
file_delete,fname+'_ch2.hdr'
file_delete,fname+'_ch2.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2_1p.dat'+' > '+fname+'_ll3_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2_2p.dat'+' > '+fname+'_ll3_2p_hdr.tbl'
file_delete,fname+'_ch2.hdr'
file_delete,fname+'_ch2_1p.dat'
file_delete,fname+'_ch2_2p.dat'

; #######
; # LL2 #
; #######

if keyword_set(no_clean) then file = hdrpath_lo+fname+'_ch2_sub_medarr_2o_1p_spect.tbl' $
	else file = hdrpath_lo+fname+'_ch2_sub_medarr_2o_1p_bcd_clean.spect.tbl' 
nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, ll2_lun, file, /get_lun
readf, ll2_lun, emptyarr
close, ll2_lun & free_lun, ll2_lun

eoh = where(strmid(emptyarr, 0, 2) eq '| ')
eoh = eoh(n_elements(eoh)-1)
ll2_arr = temporary(emptyarr(0:eoh))

ll2_co_file = '_ll2_cal.tbl'
ll2_1p_file = '_ll2_1p_cal.tbl'
ll2_2p_file = '_ll2_2p_cal.tbl'

readcol, trimpath+fname+ll2_co_file, $
	det, wave, flux, err, bit, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+ll2_1p_file, $
	det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', skipline = 1, /silent

readcol, nodpath+fname+ll2_2p_file, $
	det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', skipline = 1, /silent


forprint, ll2_arr, format = 'a', textout = writepath+fname+'_ch2.hdr', /nocomment, /silent
forprint, det, wave, flux, err, bit, format = 'i,f,f,f,i', textout = writepath+fname+'_ch2.dat', /nocomment, /silent

forprint, ll2_arr, format = 'a', textout = nwritepath+fname+'_ch2.hdr', /nocomment, /silent
forprint, det_1p, wave_1p, flux_1p, err_1p, bit_1p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch2_1p.dat', /nocomment, /silent
forprint, det_2p, wave_2p, flux_2p, err_2p, bit_2p, format = 'i,f,f,f,i', textout = nwritepath+fname+'_ch2_2p.dat', /nocomment, /silent

cd, writepath

spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2.dat'+' > '+fname+'_ll2_hdr.tbl'
file_delete,fname+'_ch2.hdr'
file_delete,fname+'_ch2.dat'

cd, nwritepath

spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2_1p.dat'+' > '+fname+'_ll2_1p_hdr.tbl'
spawn, 'cat '+fname+'_ch2.hdr'+' '+fname+'_ch2_2p.dat'+' > '+fname+'_ll2_2p_hdr.tbl'
file_delete,fname+'_ch2.hdr'
file_delete,fname+'_ch2_1p.dat'
file_delete,fname+'_ch2_2p.dat'

print,'Completed SPICEHDR for PU-calibrated '+fname

;####################
;####################

; End program

end
