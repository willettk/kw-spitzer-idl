pro spec2pahfit, fname, bw = bw, nods=nods, ps = ps, nolines = nolines, $
	xr = xr, yr = yr, pahfit = pahfit, write = write, notrim = notrim, ss = ss, noclean = noclean, $
	regular = regular, nobonus = nobonus, panel=panel, nolabel = nolabel, nolegend = nolegend, $
	no_offset = no_offset, quiet = quiet
;+
; NAME:
;       SPEC2PAHFIT
;
; PURPOSE:
; 	Display trimmed lo-res Spitzer IRS spectra with option of fitting ULIRG template via PAHFIT
;
; INPUTS:
;	XR - plot range in x-direction (log scale)
;
;	YR - plot range in y-direction (log scale)
;
; OUTPUTS:
; 	- Plots low-res spectra sorted by order in rest frame
;
; KEYWORDS:
;
;	BW - displays spectra in black and white. Default is to display each order in a different color.
;
;	NOLINES - removes common mid-IR lines with labels. Default is to display them. 
;
;	NODS - displays individual nods vertically offset from each other (nod 2 on top). Default is
;		to display coadded spectra weighted by flux uncertainty in each point.
;
;	PS - plots the spectrum to a hardcopy postscript file
;
;	WRITE - write full, trimmed spectrum to ASCII file (used for SPEC2PUFLUX)
;
;	NOCLEAN - view data (if available) that has not gone through IRSCLEAN_MASK. Default
;			is to view cleaned data. 
;
;	NOTRIM - display all orders without trimming on any ends. Default is to use custom templates
;		for each object (if available) or a standard template based on IRAS 16255+2801 (mega023). 
;
;	REGULAR - view and write data extracted using nominal settings in SPICE (if present). Default
;			is to use the optimally extracted data. 
;
;	PAHFIT - run JD Smith's PAHFIT routine to fit PAH features, emission lines, BB, and continua as
;		method of decomposing SED.
;
;	SS - runs the program on spectra that have been extracted with a supersky background (obsolete)
;
;	NOBONUS - does not plot the bonus orders SL3 and LL3
;
; EXAMPLE:
;	IDL> spec2pahfit, 'mega023'
;
; REQUIRES:
;
;	PAHFIT.pro (and associated routines)
;	READCOL.pro
;	TAG.pro
;	TARGETS.pro
;
; NOTES:
;
;	This program displays the spectra as extracted from SPICE following trimming - the fluxes have
;	not yet been corrected using the peakups.
;
;	Note that PAHFIT.pro must be compiled before running this routine, regardless of whether the keyword PAHFIT
;	is called (somewhat annoying). 
;
; REVISION HISTORY
;       Written by K. Willett                Apr 2007
;	Added comments, COADD, COLOR, PS, LINES keywords - KW, May 07
;	Added WRITE, PAHFIT keywords - KW, May 07
;	Reference lookup table in TARGETS.pro instead of hard-coding it - KW, May 07
;	Added more lines at shortest wavelengths - KW, May 07
; 	Added NOTRIM keyword - KW, Jun 07
;	Added CLEAN/NOCLEAN keyword - KW, Jul 07
;	Added lookup tables for trim templates at order edges - KW, Jul 07
;	Added bonus order - KW, Aug 07
;	Changed defaults to run coadded, colored, line-IDed, no stop spectra - KW, Sep 07
;	Added NOBONUS keyword - KW, Dec 07
;	Negative pixels in only one nod are not used in the co-added spectrum - KW, Jul 09
;-

; Set device to read in colors

;if not keyword_set(ps) then device, decomposed = 1, window_state = state 
;device, window_state = state

; Locate directory from which to read (either megamasers or CSOs)

tag, fname, dirtag 

; Read in the optimally extracted data from SPICE directory

if keyword_set(regular) then begin
	writedir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/regular/'
	plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/regular/'
	opath = '~/spice/output/lores/regular/all/' 
	pahpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/regular/'
endif else begin
	writedir = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_spectra/'
	plotdir = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/'
	opath = '~/spice/output/lores/optimal/all/'
	pahpath = '~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/'
endelse

if keyword_set(noclean) then files_opt = opath+'*'+fname+'*p_spect.tbl'$
	else files_opt = opath+'*'+fname+'*clean*spect.tbl'

ofiles = findfile(files_opt)
onum = n_elements(ofiles)

; Stored object designations, redshifts for the Spitzer sample

targets, fname, redshift, obj

; Read information in from tables for both nod positions

readcol, ofiles(0), sl1_1p_order, sl1_1p_wave, sl1_1p_flux, sl1_1p_error, sl1_1p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(1), sl1_2p_order, sl1_2p_wave, sl1_2p_flux, sl1_2p_error, sl1_2p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(2), sl2_1p_order, sl2_1p_wave, sl2_1p_flux, sl2_1p_error, sl2_1p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(3), sl2_2p_order, sl2_2p_wave, sl2_2p_flux, sl2_2p_error, sl2_2p_bit, format = 'i,f,f,f,i', /silent

readcol, ofiles(4), ll1_1p_order, ll1_1p_wave, ll1_1p_flux, ll1_1p_error, ll1_1p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(5), ll1_2p_order, ll1_2p_wave, ll1_2p_flux, ll1_2p_error, ll1_2p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(6), ll2_1p_order, ll2_1p_wave, ll2_1p_flux, ll2_1p_error, ll2_1p_bit, format = 'i,f,f,f,i', /silent
readcol, ofiles(7), ll2_2p_order, ll2_2p_wave, ll2_2p_flux, ll2_2p_error, ll2_2p_bit, format = 'i,f,f,f,i', /silent

; Remove bonus order from SL2, LL2 nods while keeping the SL3 data

nosl3_1p = where(sl2_1p_order eq 2)
nosl3_2p = where(sl2_2p_order eq 2)

sl3_1p_ind = where(sl2_1p_order eq 3)
sl3_2p_ind = where(sl2_2p_order eq 3)

sl3_1p_wave = sl2_1p_wave(sl3_1p_ind) & sl3_1p_flux = sl2_1p_flux(sl3_1p_ind) & sl3_1p_error = sl2_1p_error(sl3_1p_ind) 
sl3_1p_bit = sl2_1p_bit(sl3_1p_ind) & sl3_1p_order = sl2_1p_order(sl3_1p_ind)
sl3_2p_wave = sl2_2p_wave(sl3_2p_ind) & sl3_2p_flux = sl2_2p_flux(sl3_2p_ind) & sl3_2p_error = sl2_2p_error(sl3_2p_ind) 
sl3_2p_bit = sl2_2p_bit(sl3_2p_ind) & sl3_2p_order = sl2_2p_order(sl3_2p_ind)

sl2_1p_wave = sl2_1p_wave(nosl3_1p) & sl2_1p_flux = sl2_1p_flux(nosl3_1p) & sl2_1p_error = sl2_1p_error(nosl3_1p) 
sl2_1p_bit = sl2_1p_bit(nosl3_1p) & sl2_1p_order = sl2_1p_order(nosl3_1p)
sl2_2p_wave = sl2_2p_wave(nosl3_2p) & sl2_2p_flux = sl2_2p_flux(nosl3_2p) & sl2_2p_error = sl2_2p_error(nosl3_2p) 
sl2_2p_bit = sl2_2p_bit(nosl3_2p) & sl2_2p_order = sl2_2p_order(nosl3_2p)

noll3_1p = where(ll2_1p_order eq 2)
noll3_2p = where(ll2_2p_order eq 2)

ll3_1p_ind = where(ll2_1p_order eq 3)
ll3_2p_ind = where(ll2_2p_order eq 3)

ll3_1p_wave = ll2_1p_wave(ll3_1p_ind) & ll3_1p_flux = ll2_1p_flux(ll3_1p_ind) & ll3_1p_error = ll2_1p_error(ll3_1p_ind) 
ll3_1p_bit = ll2_1p_bit(ll3_1p_ind) & ll3_1p_order = ll2_1p_order(ll3_1p_ind)
ll3_2p_wave = ll2_2p_wave(ll3_2p_ind) & ll3_2p_flux = ll2_2p_flux(ll3_2p_ind) & ll3_2p_error = ll2_2p_error(ll3_2p_ind) 
ll3_2p_bit = ll2_2p_bit(ll3_2p_ind) & ll3_2p_order = ll2_2p_order(ll3_2p_ind)

ll2_1p_wave = ll2_1p_wave(noll3_1p) & ll2_1p_flux = ll2_1p_flux(noll3_1p) & ll2_1p_error = ll2_1p_error(noll3_1p) 
ll2_1p_bit = ll2_1p_bit(noll3_1p) & ll2_1p_order = ll2_1p_order(noll3_1p)
ll2_2p_wave = ll2_2p_wave(noll3_2p) & ll2_2p_flux = ll2_2p_flux(noll3_2p) & ll2_2p_error = ll2_2p_error(noll3_2p) 
ll2_2p_bit = ll2_2p_bit(noll3_2p) & ll2_2p_order = ll2_2p_order(noll3_2p)


; Sort the pixels so that they are in order by wavelength

sl1_1p_index = sort(sl1_1p_wave)
sl2_1p_index = sort(sl2_1p_wave)
sl3_1p_index = sort(sl3_1p_wave)
ll1_1p_index = sort(ll1_1p_wave)
ll2_1p_index = sort(ll2_1p_wave)
ll3_1p_index = sort(ll3_1p_wave)

sl1_2p_index = sort(sl1_2p_wave)
sl2_2p_index = sort(sl2_2p_wave)
sl3_2p_index = sort(sl3_2p_wave)
ll1_2p_index = sort(ll1_2p_wave)
ll2_2p_index = sort(ll2_2p_wave)
ll3_2p_index = sort(ll3_2p_wave)

; Trim the ends of the modules - template based on SL spectrum for mega023 (IRAS 16255+2801)
; No trimming of the bonus orders for now

trim_file = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/trim_templates/'+fname+'_lores.trim')
	if trim_file(0) eq '' then begin

		sl2_edges = [2,2]
		sl3_edges = [0,2]
		sl1_edges = [2,15]
		ll2_edges = [3,15]
		ll3_edges = [0,1]
		ll1_edges = [3,25]

		if not keyword_set(quiet) then print,'Using standard template for lo-res modules on '+fname

	endif else begin

		if not keyword_set(quiet) then print,'Using saved template for lo-res modules on '+fname
		readcol, trim_file(0), trimstart, trimend, format = 'x,i,i', skipline = 2, /silent

		sl2_edges = [trimstart(0),trimend(0)]	
		sl3_edges = [trimstart(1),trimend(1)]	
		sl1_edges = [trimstart(2),trimend(2)]	
		ll2_edges = [trimstart(3),trimend(3)]	
		ll3_edges = [trimstart(4),trimend(4)]	
		ll1_edges = [trimstart(5),trimend(5)]	

	endelse


if keyword_set(notrim) then begin
	sl1_edges = [0,0]
	sl2_edges = [0,0]
	sl3_edges = [0,0]
	ll1_edges = [0,0]
	ll2_edges = [0,0]
	ll3_edges = [0,0]
endif

sl1_1p_index = sl1_1p_index(sl1_edges(0):n_elements(sl1_1p_index)-sl1_edges(1)-1)	
sl2_1p_index = sl2_1p_index(sl2_edges(0):n_elements(sl2_1p_index)-sl2_edges(1)-1)
sl3_1p_index = sl3_1p_index(sl3_edges(0):n_elements(sl3_1p_index)-sl3_edges(1)-1)
ll1_1p_index = ll1_1p_index(ll1_edges(0):n_elements(ll1_1p_index)-ll1_edges(1)-1)
ll2_1p_index = ll2_1p_index(ll2_edges(0):n_elements(ll2_1p_index)-ll2_edges(1)-1)
ll3_1p_index = ll3_1p_index(ll3_edges(0):n_elements(ll3_1p_index)-ll3_edges(1)-1)

sl1_2p_index = sl1_2p_index(sl1_edges(0):n_elements(sl1_2p_index)-sl1_edges(1)-1)	
sl2_2p_index = sl2_2p_index(sl2_edges(0):n_elements(sl2_2p_index)-sl2_edges(1)-1)
sl3_2p_index = sl3_2p_index(sl3_edges(0):n_elements(sl3_2p_index)-sl3_edges(1)-1)
ll1_2p_index = ll1_2p_index(ll1_edges(0):n_elements(ll1_2p_index)-ll1_edges(1)-1)
ll2_2p_index = ll2_2p_index(ll2_edges(0):n_elements(ll2_2p_index)-ll2_edges(1)-1)
ll3_2p_index = ll3_2p_index(ll3_edges(0):n_elements(ll3_2p_index)-ll3_edges(1)-1)

; New orders with trimmed edges

sl1_1p_wave = sl1_1p_wave(sl1_1p_index) & sl1_1p_flux = sl1_1p_flux(sl1_1p_index) & sl1_1p_error = sl1_1p_error(sl1_1p_index) & $
	sl1_1p_order = sl1_1p_order(sl1_1p_index) & sl1_1p_bit = sl1_1p_bit(sl1_1p_index)
sl2_1p_wave = sl2_1p_wave(sl2_1p_index) & sl2_1p_flux = sl2_1p_flux(sl2_1p_index) & sl2_1p_error = sl2_1p_error(sl2_1p_index) & $
	sl2_1p_order = sl2_1p_order(sl2_1p_index) & sl2_1p_bit = sl2_1p_bit(sl2_1p_index)
sl3_1p_wave = sl3_1p_wave(sl3_1p_index) & sl3_1p_flux = sl3_1p_flux(sl3_1p_index) & sl3_1p_error = sl3_1p_error(sl3_1p_index) & $
	sl3_1p_order = sl3_1p_order(sl3_1p_index) & sl3_1p_bit = sl3_1p_bit(sl3_1p_index)
ll1_1p_wave = ll1_1p_wave(ll1_1p_index) & ll1_1p_flux = ll1_1p_flux(ll1_1p_index) & ll1_1p_error = ll1_1p_error(ll1_1p_index) & $
	ll1_1p_order = ll1_1p_order(ll1_1p_index) & ll1_1p_bit = ll1_1p_bit(ll1_1p_index)
ll2_1p_wave = ll2_1p_wave(ll2_1p_index) & ll2_1p_flux = ll2_1p_flux(ll2_1p_index) & ll2_1p_error = ll2_1p_error(ll2_1p_index) & $
	ll2_1p_order = ll2_1p_order(ll2_1p_index) & ll2_1p_bit = ll2_1p_bit(ll2_1p_index)
ll3_1p_wave = ll3_1p_wave(ll3_1p_index) & ll3_1p_flux = ll3_1p_flux(ll3_1p_index) & ll3_1p_error = ll3_1p_error(ll3_1p_index) & $
	ll3_1p_order = ll3_1p_order(ll3_1p_index) & ll3_1p_bit = ll3_1p_bit(ll3_1p_index)

sl1_2p_wave = sl1_2p_wave(sl1_2p_index) & sl1_2p_flux = sl1_2p_flux(sl1_2p_index) & sl1_2p_error = sl1_2p_error(sl1_2p_index) & $
	sl1_2p_order = sl1_2p_order(sl1_2p_index) & sl1_2p_bit = sl1_2p_bit(sl1_2p_index)
sl2_2p_wave = sl2_2p_wave(sl2_2p_index) & sl2_2p_flux = sl2_2p_flux(sl2_2p_index) & sl2_2p_error = sl2_2p_error(sl2_2p_index) & $
	sl2_2p_order = sl2_2p_order(sl2_2p_index) & sl2_2p_bit = sl2_2p_bit(sl2_2p_index)
sl3_2p_wave = sl3_2p_wave(sl3_2p_index) & sl3_2p_flux = sl3_2p_flux(sl3_2p_index) & sl3_2p_error = sl3_2p_error(sl3_2p_index) & $
	sl3_2p_order = sl3_2p_order(sl3_2p_index) & sl3_2p_bit = sl3_2p_bit(sl3_2p_index)
ll1_2p_wave = ll1_2p_wave(ll1_2p_index) & ll1_2p_flux = ll1_2p_flux(ll1_2p_index) & ll1_2p_error = ll1_2p_error(ll1_2p_index) & $
	ll1_2p_order = ll1_2p_order(ll1_2p_index) & ll1_2p_bit = ll1_2p_bit(ll1_2p_index)
ll2_2p_wave = ll2_2p_wave(ll2_2p_index) & ll2_2p_flux = ll2_2p_flux(ll2_2p_index) & ll2_2p_error = ll2_2p_error(ll2_2p_index) & $
	ll2_2p_order = ll2_2p_order(ll2_2p_index) & ll2_2p_bit = ll2_2p_bit(ll2_2p_index)
ll3_2p_wave = ll3_2p_wave(ll3_2p_index) & ll3_2p_flux = ll3_2p_flux(ll3_2p_index) & ll3_2p_error = ll3_2p_error(ll3_2p_index) & $
	ll3_2p_order = ll3_2p_order(ll3_2p_index) & ll3_2p_bit = ll3_2p_bit(ll3_2p_index)

; Collapse orders into one large spectrum sorted by wavelength

allflux_1p   = [sl2_1p_flux, sl3_1p_flux, sl1_1p_flux, ll2_1p_flux, ll3_1p_flux, ll1_1p_flux]
allwave_1p   = [sl2_1p_wave, sl3_1p_wave, sl1_1p_wave, ll2_1p_wave, ll3_1p_wave, ll1_1p_wave]
allorder_1p  = [sl2_1p_order, sl3_1p_order, sl1_1p_order, ll2_1p_order, ll3_1p_order, ll1_1p_order]
allerror_1p  = [sl2_1p_error, sl3_1p_error, sl1_1p_error, ll2_1p_error, ll3_1p_error, ll1_1p_error]
allbit_1p    = [sl2_1p_bit, sl3_1p_bit, sl1_1p_bit, ll2_1p_bit, ll3_1p_bit, ll1_1p_bit]

allflux_2p   = [sl2_2p_flux, sl3_2p_flux, sl1_2p_flux, ll2_2p_flux, ll3_2p_flux, ll1_2p_flux]
allwave_2p   = [sl2_2p_wave, sl3_2p_wave, sl1_2p_wave, ll2_2p_wave, ll3_2p_wave, ll1_2p_wave]
allorder_2p  = [sl2_2p_order, sl3_2p_order, sl1_2p_order, ll2_2p_order, ll3_2p_order, ll1_2p_order]
allerror_2p  = [sl2_2p_error, sl3_2p_error, sl1_2p_error, ll2_2p_error, ll3_2p_error, ll1_2p_error]
allbit_2p    = [sl2_2p_bit, sl3_2p_bit, sl1_2p_bit, ll2_2p_bit, ll3_2p_bit, ll1_2p_bit]

allwave_1p  = allwave_1p(sort(allwave_1p))
allflux_1p  = allflux_1p(sort(allwave_1p))
allorder_1p = allorder_1p(sort(allwave_1p))
allerror_1p = allerror_1p(sort(allwave_1p))
allbit_1p   = allbit_1p(sort(allwave_1p))

allwave_2p  = allwave_2p(sort(allwave_2p))
allflux_2p  = allflux_2p(sort(allwave_2p))
allorder_2p = allorder_2p(sort(allwave_2p))
allerror_2p = allerror_2p(sort(allwave_2p))
allbit_2p   = allbit_2p(sort(allwave_2p))

; Coadd the nod positions using a weighted mean 

weights_sl1_1p = 1d-6/sl1_1p_error^2
weights_sl2_1p = 1d-6/sl2_1p_error^2
weights_sl3_1p = 1d-6/sl3_1p_error^2
weights_ll1_1p = 1d-6/ll1_1p_error^2
weights_ll2_1p = 1d-6/ll2_1p_error^2
weights_ll3_1p = 1d-6/ll3_1p_error^2

weights_sl1_2p = 1d-6/sl1_2p_error^2
weights_sl2_2p = 1d-6/sl2_2p_error^2
weights_sl3_2p = 1d-6/sl3_2p_error^2
weights_ll1_2p = 1d-6/ll1_2p_error^2
weights_ll2_2p = 1d-6/ll2_2p_error^2
weights_ll3_2p = 1d-6/ll3_2p_error^2

sl1_coadd = (weights_sl1_1p * sl1_1p_flux + weights_sl1_2p * sl1_2p_flux) / (weights_sl1_1p + weights_sl1_2p)
sl2_coadd = (weights_sl2_1p * sl2_1p_flux + weights_sl2_2p * sl2_2p_flux) / (weights_sl2_1p + weights_sl2_2p)
sl3_coadd = (weights_sl3_1p * sl3_1p_flux + weights_sl3_2p * sl3_2p_flux) / (weights_sl3_1p + weights_sl3_2p)
ll1_coadd = (weights_ll1_1p * ll1_1p_flux + weights_ll1_2p * ll1_2p_flux) / (weights_ll1_1p + weights_ll1_2p)
ll2_coadd = (weights_ll2_1p * ll2_1p_flux + weights_ll2_2p * ll2_2p_flux) / (weights_ll2_1p + weights_ll2_2p)
ll3_coadd = (weights_ll3_1p * ll3_1p_flux + weights_ll3_2p * ll3_2p_flux) / (weights_ll3_1p + weights_ll3_2p)

sl1_coadd_err = sqrt(sl1_1p_error^2 + sl1_2p_error^2)
sl2_coadd_err = sqrt(sl2_1p_error^2 + sl2_2p_error^2)
sl3_coadd_err = sqrt(sl3_1p_error^2 + sl3_2p_error^2)
ll1_coadd_err = sqrt(ll1_1p_error^2 + ll1_2p_error^2)
ll2_coadd_err = sqrt(ll2_1p_error^2 + ll2_2p_error^2)
ll3_coadd_err = sqrt(ll3_1p_error^2 + ll3_2p_error^2)

; If pixel flux for one nod is negative, use only positive flux from the other nod

posflux_sl1_1p = where(sl1_1p_flux gt 0. and sl1_2p_flux lt 0.,sl1_1p_poscount)
posflux_sl1_2p = where(sl1_1p_flux lt 0. and sl1_2p_flux gt 0.,sl1_2p_poscount)
if sl1_1p_poscount gt 0 then begin
	sl1_coadd[posflux_sl1_1p] = sl1_1p_flux[posflux_sl1_1p]
	sl1_coadd_err[posflux_sl1_1p] = sl1_1p_error[posflux_sl1_1p]
endif
if sl1_2p_poscount gt 0 then begin
	sl1_coadd[posflux_sl1_2p] = sl1_2p_flux[posflux_sl1_2p]
	sl1_coadd_err[posflux_sl1_2p] = sl1_2p_error[posflux_sl1_2p]
endif

posflux_sl2_1p = where(sl2_1p_flux gt 0. and sl2_2p_flux lt 0.,sl2_1p_poscount)
posflux_sl2_2p = where(sl2_1p_flux lt 0. and sl2_2p_flux gt 0.,sl2_2p_poscount)
if sl2_1p_poscount gt 0 then begin
	sl2_coadd[posflux_sl2_1p] = sl2_1p_flux[posflux_sl2_1p]
	sl2_coadd_err[posflux_sl2_1p] = sl2_1p_error[posflux_sl2_1p]
endif
if sl2_2p_poscount gt 0 then begin
	sl2_coadd[posflux_sl2_2p] = sl2_2p_flux[posflux_sl2_2p]
	sl2_coadd_err[posflux_sl2_2p] = sl2_2p_error[posflux_sl2_2p]
endif

posflux_sl3_1p = where(sl3_1p_flux gt 0. and sl3_2p_flux lt 0.,sl3_1p_poscount)
posflux_sl3_2p = where(sl3_1p_flux lt 0. and sl3_2p_flux gt 0.,sl3_2p_poscount)
if sl3_1p_poscount gt 0 then begin
	sl3_coadd[posflux_sl3_1p] = sl3_1p_flux[posflux_sl3_1p]
	sl3_coadd_err[posflux_sl3_1p] = sl3_1p_error[posflux_sl3_1p]
endif
if sl3_2p_poscount gt 0 then begin
	sl3_coadd[posflux_sl3_2p] = sl3_2p_flux[posflux_sl3_2p]
	sl3_coadd_err[posflux_sl3_2p] = sl3_2p_error[posflux_sl3_2p]
endif

posflux_ll1_1p = where(ll1_1p_flux gt 0. and ll1_2p_flux lt 0.,ll1_1p_poscount)
posflux_ll1_2p = where(ll1_1p_flux lt 0. and ll1_2p_flux gt 0.,ll1_2p_poscount)
if ll1_1p_poscount gt 0 then begin
	ll1_coadd[posflux_ll1_1p] = ll1_1p_flux[posflux_ll1_1p]
	ll1_coadd_err[posflux_ll1_1p] = ll1_1p_error[posflux_ll1_1p]
endif
if ll1_2p_poscount gt 0 then begin
	ll1_coadd[posflux_ll1_2p] = ll1_2p_flux[posflux_ll1_2p]
	ll1_coadd_err[posflux_ll1_2p] = ll1_2p_error[posflux_ll1_2p]
endif

posflux_ll2_1p = where(ll2_1p_flux gt 0. and ll2_2p_flux lt 0.,ll2_1p_poscount)
posflux_ll2_2p = where(ll2_1p_flux lt 0. and ll2_2p_flux gt 0.,ll2_2p_poscount)
if ll2_1p_poscount gt 0 then begin
	ll2_coadd[posflux_ll2_1p] = ll2_1p_flux[posflux_ll2_1p]
	ll2_coadd_err[posflux_ll2_1p] = ll2_1p_error[posflux_ll2_1p]
endif
if ll2_2p_poscount gt 0 then begin
	ll2_coadd[posflux_ll2_2p] = ll2_2p_flux[posflux_ll2_2p]
	ll2_coadd_err[posflux_ll2_2p] = ll2_2p_error[posflux_ll2_2p]
endif

posflux_ll3_1p = where(ll3_1p_flux gt 0. and ll3_2p_flux lt 0.,ll3_1p_poscount)
posflux_ll3_2p = where(ll3_1p_flux lt 0. and ll3_2p_flux gt 0.,ll3_2p_poscount)
if ll3_1p_poscount gt 0 then begin
	ll3_coadd[posflux_ll3_1p] = ll3_1p_flux[posflux_ll3_1p]
	ll3_coadd_err[posflux_ll3_1p] = ll3_1p_error[posflux_ll3_1p]
endif
if ll3_2p_poscount gt 0 then begin
	ll3_coadd[posflux_ll3_2p] = ll3_2p_flux[posflux_ll3_2p]
	ll3_coadd_err[posflux_ll3_2p] = ll3_2p_error[posflux_ll3_2p]
endif


allflux_coadd = [sl2_coadd, sl3_coadd, sl1_coadd, ll2_coadd, ll3_coadd, ll1_coadd]
allerror_coadd = sqrt(allerror_1p^2 + allerror_2p^2)

; Create grid of likely suspects for emission/absorption lines 
; (using references from Armus et al 2004, Armus et al 2006)

linelist = [5.511, 6.1088, 6.2, 6.909, 6.98, 7.7, 8.025, 8.6, 8.99138, 9.665, $
	 10.511, 11.3, 12.279, 12.6, 12.814, $
	  14.2,  14.322, 15.555, 16.4, 17.035, 17.4, $
	18.713,  24.318, 25.890, 28.218]

line_id = ['H!I2!N S(7)', 'H!I2!N S(6)', 'PAH', 'H!I2!N S(5)', 'ArII', 'PAH', 'H!I2!N S(4)', 'PAH', 'ArIII', 'H!I2!N S(3)', $
	'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
	  'PAH',  'NeV', 'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', $
	'SIII', 'NeV', 'OIV', 'H!I2!N S(0)']

; Old lineset - includes many weaker medium-ionization lines usually not resolvable in lo-res spectra

;linelist = [5.34, 5.511, 6.2, 6.909, 6.98, 7.65, 7.7, 8.025, 8.6, 9.665, 10.511, 11.3, 12.279, 12.6, 12.814, 13.7, $
;	14.0, 14.2, 14.322, 14.368, 15.0, 15.555, 16.4, 17.035, 17.4, 17.934, 18.713, 24.318, 25.890, $
;	25.988, 28.218, 33.481, 34.815]

;line_id = ['FeII', 'H!I2!N S(7)', 'PAH', 'H!I2!N S(5)', 'ArII', 'NeVI', 'PAH', 'H!I2!N S(4)', 'PAH', 'H!I2!N S(3)', $
;	'SIV', 'PAH', 'H!I2!N S(2)', 'PAH', 'NeII', $
;	'C!I2!NH!I2!N', 'HCN', 'PAH', 'NeV', 'ClII', 'CO!I2!N', 'NeIII', 'PAH', 'H!I2!N S(1)', 'PAH', 'FeII', 'SIII', $
;	'NeV', 'OIV', 'FeII', 'H!I2!N S(0)', 'SIII', 'SiII']

; Offset factor to plot 2nd nod position

if not keyword_set(no_offset) then offset = 4. else offset = 1.

; Plot data

;if state(0) eq 0 then window,0 else wset,0
;!p.multi = [0,1,1]

; Hard copy option

if keyword_set(ps) then begin
	set_plot,'ps'
	device, /landscape, filename = plotdir+strjoin(strsplit(obj,' ',/extract))+'_spect_lr.ps', /color
	defcolor = fsc_color('Black')
	lthick = 2
	cs = 1
endif else begin
	defcolor = fsc_color('White')
	lthick = 1
	cs = 2
endelse

if not keyword_set(panel) then !p.multi=[0,1,1]

red = fsc_color("Red")
blue = fsc_color("Blue")
green = fsc_color("Green")
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")

if not keyword_set(xr) then xr = [4,42]
if not keyword_set(yr) then yr = [1d-4,max(allflux_2p*offset)]

plot, sl1_1p_wave, sl1_1p_flux, $
	/xlog, $
	/ylog, $
;	xticks = 13, $
;	xtickv = [5,6,7,8,9,10,12,14,16,18,20,25,30,35], $
	xrange = xr, /xstyle, $
	yrange = yr, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [Jy]', $
	title = obj, $
	charsize = cs, $
	color = defcolor, $
	thick = lthick, $
	charthick = lthick, $
	/nodata

if not keyword_set(bw) then begin
	if not keyword_set(nods) then begin
		oplot, sl1_1p_wave / (redshift + 1.),sl1_coadd, psym = 10, color = red, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.),sl2_coadd, psym = 10, color = blue, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.),ll1_coadd, psym = 10, color = yellow, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.),ll2_coadd, psym = 10, color = green, thick = lthick
		if not keyword_set(nobonus) then begin
			oplot, sl3_1p_wave / (redshift + 1.),sl3_coadd, psym = 10, color = orange, thick = lthick
			oplot, ll3_1p_wave / (redshift + 1.),ll3_coadd, psym = 10, color = orange, thick = lthick
		endif
		if not keyword_set(nolabel) then xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
		if not keyword_set(nolegend) then legend, ['SL1', 'SL2', 'LL1', 'LL2', 'Bonus orders'], linestyle = intarr(5), $
			color = [red, blue, yellow, green, orange], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
	endif else begin
		if not keyword_set(nolegend) then legend, ['SL1, pos 1', 'SL2, pos 1', 'LL1, pos 1', 'LL2, pos 1', 'Bonus orders, pos 1', $
			'SL1, pos 2', 'SL2, pos 2', 'LL1, pos 2', 'LL2, pos 2', 'Bonus orders, pos 2'], $
			linestyle = intarr(10), $
			color = [red, blue, yellow, green, orange, red, blue, yellow, green, orange], $
			/top, /left, charsize = cs, textcolors = defcolor, outline_color = defcolor
		oplot, sl1_1p_wave / (redshift + 1.), sl1_1p_flux, psym = 10, color = red, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.), sl2_1p_flux, psym = 10, color = blue, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.), ll1_1p_flux, psym = 10, color = yellow, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_1p_flux, psym = 10, color = green, thick = lthick
		oplot, sl1_2p_wave / (redshift + 1.), sl1_2p_flux*offset, psym = 10, color = red, linestyle = 0, thick = lthick
		oplot, sl2_2p_wave / (redshift + 1.), sl2_2p_flux*offset, psym = 10, color = blue, linestyle = 0, thick = lthick
		oplot, ll1_2p_wave / (redshift + 1.), ll1_2p_flux*offset, psym = 10, color = yellow, linestyle = 0, thick = lthick
		oplot, ll2_2p_wave / (redshift + 1.), ll2_2p_flux*offset, psym = 10, color = green, linestyle = 0, thick = lthick
		if not keyword_set(nobonus) then begin
			oplot, sl3_1p_wave / (redshift + 1.), sl3_1p_flux, psym = 10, color = orange, thick = lthick
			oplot, ll3_1p_wave / (redshift + 1.), ll3_1p_flux, psym = 10, color = orange, thick = lthick
			oplot, sl3_2p_wave / (redshift + 1.), sl3_2p_flux*offset, psym = 10, color = orange, linestyle = 0, thick = lthick
			oplot, ll3_2p_wave / (redshift + 1.), ll3_2p_flux*offset, psym = 10, color = orange, linestyle = 0, thick = lthick
		endif
	endelse
endif else begin
	if not keyword_set(nods) then begin
		oplot, sl1_1p_wave / (redshift + 1.), sl1_coadd, psym = 10, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.), sl2_coadd, psym = 10, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.), ll1_coadd, psym = 10, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_coadd, psym = 10, thick = lthick
		if not keyword_set(nobonus) then begin
			oplot, sl3_1p_wave / (redshift + 1.), sl3_coadd, psym = 10, thick = lthick
			oplot, ll3_1p_wave / (redshift + 1.), ll3_coadd, psym = 10, thick = lthick
		endif
		if not keyword_set(nolabel) then xyouts, 0.3, 0.9, /normal, 'Coadded', charsize = cs
	endif else begin
		oplot, sl1_1p_wave / (redshift + 1.), sl1_1p_flux, psym = 10, thick = lthick
		oplot, sl2_1p_wave / (redshift + 1.), sl2_1p_flux, psym = 10, thick = lthick
		oplot, ll1_1p_wave / (redshift + 1.), ll1_1p_flux, psym = 10, thick = lthick
		oplot, ll2_1p_wave / (redshift + 1.), ll2_1p_flux, psym = 10, thick = lthick
		oplot, sl1_2p_wave / (redshift + 1.), sl1_2p_flux*offset, psym = 10, linestyle = 0, thick = lthick
		oplot, sl2_2p_wave / (redshift + 1.), sl2_2p_flux*offset, psym = 10, linestyle = 0, thick = lthick
		oplot, ll1_2p_wave / (redshift + 1.), ll1_2p_flux*offset, psym = 10, linestyle = 0, thick = lthick
		oplot, ll2_2p_wave / (redshift + 1.), ll2_2p_flux*offset, psym = 10, linestyle = 0, thick = lthick
		if not keyword_set(nobonus) then begin
			oplot, sl3_1p_wave / (redshift + 1.), sl3_1p_flux, psym = 10, thick = lthick
			oplot, ll3_1p_wave / (redshift + 1.), ll3_1p_flux, psym = 10, thick = lthick
			oplot, sl3_2p_wave / (redshift + 1.), sl3_2p_flux*offset, psym = 10, linestyle = 1, thick = lthick
			oplot, ll3_2p_wave / (redshift + 1.), ll3_2p_flux*offset, psym = 10, linestyle = 0, thick = lthick
		endif
	endelse
endelse

if not keyword_set(nolines) then begin
	for i = 0, n_elements(linelist) - 1 do begin
		off = i mod 2
		ver, linelist(i), linestyle = 1, color = defcolor
		xyouts, linelist(i), 0.35*yr(1) + 1.05*off, line_id(i), orientation = 90, charsize = cs, /data
	endfor
endif

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif

; Write the new, trimmed spectra to file

if keyword_set(write) then begin

	; All modules

	forprint, allorder_1p, allwave_1p, allflux_1p, allerror_1p, allbit_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]      Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_lr_1p.tbl', /silent

	forprint, allorder_2p, allwave_2p, allflux_2p, allerror_2p, allbit_2p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]      Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_lr_2p.tbl', /silent

	forprint, allorder_1p, allwave_1p, allflux_coadd, allerror_coadd, allbit_1p, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]      Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_lr_coadd.tbl', /silent

	; Individual modules

	forprint, sl1_1p_order, sl1_1p_wave, sl1_1p_flux, sl1_1p_error, sl1_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sl1_1p.tbl', /silent

	forprint, sl2_1p_order, sl2_1p_wave, sl2_1p_flux, sl2_1p_error, sl2_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sl2_1p.tbl', /silent

	forprint, sl3_1p_order, sl3_1p_wave, sl3_1p_flux, sl3_1p_error, sl3_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sl3_1p.tbl', /silent

	forprint, ll1_1p_order, ll1_1p_wave, ll1_1p_flux, ll1_1p_error, ll1_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll1_1p.tbl', /silent

	forprint, ll2_1p_order, ll2_1p_wave, ll2_1p_flux, ll2_1p_error, ll2_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll2_1p.tbl', /silent

	forprint, ll3_1p_order, ll3_1p_wave, ll3_1p_flux, ll3_1p_error, ll3_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll3_1p.tbl', /silent

	forprint, sl1_2p_order, sl1_2p_wave, sl1_2p_flux, sl1_2p_error, sl1_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sl1_2p.tbl', /silent

	forprint, sl2_2p_order, sl2_2p_wave, sl2_2p_flux, sl2_2p_error, sl2_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sl2_2p.tbl', /silent

	forprint, sl3_2p_order, sl3_2p_wave, sl3_2p_flux, sl3_2p_error, sl3_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_sl3_2p.tbl', /silent

	forprint, ll1_2p_order, ll1_2p_wave, ll1_2p_flux, ll1_2p_error, ll1_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll1_2p.tbl', /silent

	forprint, ll2_2p_order, ll2_2p_wave, ll2_2p_flux, ll2_2p_error, ll2_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll2_2p.tbl', /silent

	forprint, ll3_2p_order, ll3_2p_wave, ll3_2p_flux, ll3_2p_error, ll3_2p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'nods/'+fname+'_ll3_2p.tbl', /silent

	; Coadded

	forprint, sl1_1p_order, sl1_1p_wave, sl1_coadd, sl1_1p_error, sl1_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_sl1_coadd.tbl', /silent

	forprint, sl2_1p_order, sl2_1p_wave, sl2_coadd, sl2_1p_error, sl2_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_sl2_coadd.tbl', /silent

	forprint, sl3_1p_order, sl3_1p_wave, sl3_coadd, sl3_1p_error, sl3_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_sl3_coadd.tbl', /silent

	forprint, ll1_1p_order, ll1_1p_wave, ll1_coadd, ll1_1p_error, ll1_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_ll1_coadd.tbl', /silent

	forprint, ll2_1p_order, ll2_1p_wave, ll2_coadd, ll2_1p_error, ll2_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_ll2_coadd.tbl', /silent

	forprint, ll3_1p_order, ll3_1p_wave, ll3_coadd, ll3_1p_error, ll3_1p_bit, format = '(i2,f15.4,e20.4,e20.4,i10)', $
		comment = 'Order   Wavelength [um]      Flux [e/s]         Uncert.        Bit type', $
		textout = writedir+'coadd/'+fname+'_ll3_coadd.tbl', /silent

endif


;stop
end
