pro pu_arch, stop = stop, quiet=quiet
;+
; NAME: 
;       PU_ARCH
;
; PURPOSE:
;
;	Extract the offset peakup fluxes from a subset of targets in the archived OHM samples
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
;
; EXAMPLE:
;
;	IDL> pu_arch
;         
; MODIFICATION HISTORY:
;
;	Written by KW (adapted from PU_READ.pro by V. Charmandaris) - Nov 07
;	Updated the conversion factors with most recent nos. from IRS Data Handbook - KW, Dec 07
; 	Fixed bug - misordering of files due to non-sorted nature of archfiles array - KW, Aug 08
;-

; List of all targets to be extracted (those with position peakups are zeroed at the end) 

archfiles = [$
'arch004', 'arch005', $
'arch008', 'arch010', 'arch013', 'arch016', $
'arch018', 'arch020', $
'arch023', 'arch024', 'arch025', $
'arch026', 'arch029', 'arch031', $
'arch032', 'arch033', $
'arch036', 'arch040', 'arch045', 'arch048', $
'arch003','arch007','arch009','arch014','arch017','arch030','arch034','arch035','arch039']
;'arch012']		No peakup on arch012

filelist=['~/Astronomy/Research/Spitzer/archived/data/arch0??_spectra/ch0/bcd/*_000[0,1]_000[0,1,2,3]_*bcd.fits']

files=file_search(filelist)


; Cull list

newfiles = strarr(n_elements(files))
m = 0
for k = 0, n_elements(files) - 1 do begin
	temptag = strmid(files(k), 57, 7)
	tempindex = where(temptag eq archfiles)
	if tempindex ge 0 then begin 
		newfiles(m) = files(k)
		m = m+1
	endif
endfor
zerofiles = where(newfiles eq '')
newfiles = newfiles(0:zerofiles(0)-1)
n_targets = n_elements(newfiles)
pu_flux = fltarr(n_targets/8)
pu_filter = strarr(n_targets/8)
pu_type = bytarr(n_targets/8)
pufilelist = strarr(n_targets/8)

if not keyword_set(quiet) then begin
	print, 'Targets w/peakups found: ',n_targets / 8
	print,''
endif

dir_out='~/Astronomy/Research/Spitzer/archived/pu_arch/'
fitsdir_out='~/Astronomy/Research/Spitzer/archived/pu_arch/fits/'
phot_file='~/Astronomy/Research/Spitzer/archived/pu_arch/pu_arch_phot.txt'
centroid_file='~/Astronomy/Research/Spitzer/archived/pu_arch/pu_arch_cntrds.txt'

; Conversion factors

bf = 7.5d5			; New conversion factors from the IRS Data Handbook, v.3.1
rf = 6.8d5

;bf=1.18e6                       ; BLUE PU: [DN in 8.3888sec] /Jy : 3pix aperture
;rf=1.15e6                       ; RED PU:  [DN in 8.3888sec] /Jy : 4pix aperture

    openw, 5, phot_file
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'File created on: ',systime()
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG ', 'TARGET  ', 'FILTR', 'GRPTIME', 'X-CNTR', 'Y-CNTR', '5 px [DN]', 'Sky [DN]', 'f(mJy)',$
      format='(A7, 4X, A8, 3X, A12, 3X, A8, 3X, A14, 1X, A7, 2X, A6, 2X, A7, 2X, A6, 5X, A9, 2X, A9, 2X, A6)'
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'

   
    openw, 6, centroid_file
    printf, 6, '------------------------------------------------------------------------------------------------------------------------'
    printf, 6, 'File created on: ',systime()
    printf, 6, '------------------------------------------------------------------------------------------------------------------------'
    printf, 6, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG ', 'TARGET  ', 'FILTER',  'X-CNTR', 'Y-CNTR', 'RA (J2000)', 'Dec(J2000)', $
      format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A12, 1X, A6, 2X, A6, 2X, A6, 5X, A12,3X, A12)'
    printf, 6, '--------------------------------------------------------------------------------------------------------------'

  

;    print, '-------------------------------------------------------------------------------------------------------------------------------'
;    print, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG ', 'TARGET  ', 'FILTR', 'GRPTIME', 'X-CNTR', 'Y-CNTR', '5 px [DN]', 'Sky [DN]', 'f(mJy)',$
;      format='(A7, 4X, A8, 3X, A12, 3X, A8, 3X, A14, 1X,A6, 2X, A7, 2X, A6, 2X, A9, 2X, A9, 2X, A6)'
;    print, '-------------------------------------------------------------------------------------------------------------------------------'
 
   
for i=0, n_targets-7 do begin

    temp0=readfits(newfiles(i),hdr_0, /silent)
    temp1=readfits(newfiles(i+4),hdr_1, /silent)

; Get the AORKEY
    aorkey=sxpar(hdr_0, 'AORKEY')
; Get the Program ID
    progid='P_'+strn(sxpar(hdr_0, 'PROGID'))
; Get the target name
    object_long=sxpar(hdr_0, 'OBJECT')
; Find units used
    bunit = sxpar(hdr_0, 'BUNIT')
; Get the filter used
    fovname=sxpar(hdr_0, 'FOVNAME')
    filter = STRMID(fovname,4,3)
; Get the exposure times 
    grptime=sxpar(hdr_0, 'GRPTIME')
; Log the day of observations
    julianday=sxpar(hdr_0, 'MJD_OBS')+2400000.5
    daycnv, julianday, year,month,day
    date=strn(YMD2DATE(year,month,day))
; Find whether the PU target is absolute or relative
	ispupos = sxpar(hdr_0,'ISPUPOS')
; Get tag from filename
	tag = strmid(newfiles(i),57,7)

;        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),''), '?',/extract),'')

	targets, tag, r, o
	object = strjoin(strsplit(o,' ',/extract))
	nameout=progid+'_'+STRMID(object,0,18)+'_'+STRTRIM(filter,2)+'_'+STRMID(newfiles(i), 57, 7)



; Check if there's WCS in the header
    is_wcs_in=sxpar(hdr_0, 'CD1_1')

    
; read images from each dither position

        dither0_0=readfits(newfiles(i), /silent)
        dither0_1=readfits(newfiles(i+1), /silent)
        dither0_2=readfits(newfiles(i+2), /silent)

        dither1_0=readfits(newfiles(i+4), /silent)
        dither1_1=readfits(newfiles(i+5), /silent)
        dither1_2=readfits(newfiles(i+6), /silent)

	; Medianing for CR rejection

        pu_d0_clean=median([[[dither0_0]], [[dither0_1]], [[dither0_2]]], dimension = 3)
        pu_d1_clean=median([[[dither1_0]], [[dither1_1]], [[dither1_2]]], dimension = 3)
        

; Eyeball glance seems to indicate that the image has been dejailbar-ed already
; Look at image processing in SOM, or think of a quantitative test I could do 
; to check on the bias levels (mean value of different sky pixel channels?)

        dejailbar, pu_d0_clean, pu_d0_dj
        dejailbar, pu_d1_clean, pu_d1_dj


; Photometry

		jyconv = 0.01617d

            if strn(filter) eq 'Red' then begin				; Red filter
                pos_x_0=106
                pos_y_0=93
		
                pos_x_1=107
                pos_y_1=95

	        aperture=[4]
		
                dn_to_mJy=4.6 / (rf * grptime) * 1d3			; See p. 43 of the IRS Data Handbook, v. 3.1
            endif else begin						; Blue filter
                pos_x_0=108
                pos_y_0=31

                pos_x_1=109
                pos_y_1=29

        	aperture=[3]

                dn_to_mJy=4.6 / (bf * grptime) * 1d3
            endelse

	    if tag eq 'arch033' then begin				; Misses the centroiding for one frame
		    pos_x_0 = 107 
		    pos_y_0 = 29
	    endif

                
; Sky annulus [8 - 14]
            skyrad=[8,14]
            badpix=[-32765,65000]
            phpadu=1
;  Centroid is computed using a box of half width equal to 
; 1.5 sigma = 0.637* FWHM
            fwhm=2.

; Photometry of ACQ image

            cntrd, pu_d0_dj, pos_x_0, pos_y_0, xcen_0, ycen_0, fwhm
            aper, pu_d0_dj, xcen_0, ycen_0, mags_0, errap_0, sky_0, skyerr_0, phpadu, aperture, skyrad, badpix, $
              /flux, /silent

            cntrd, pu_d1_dj, pos_x_1, pos_y_1, xcen_1, ycen_1, fwhm
            aper, pu_d1_dj, xcen_1, ycen_1, mags_1, errap_1, sky_1, skyerr_1, phpadu, aperture, skyrad, badpix, $
              /flux, /silent



		; Convert DN to flux

		flux_0 = mags_0 * dn_to_mJy
		flux_1 = mags_1 * dn_to_mJy

		;flux_0 = mags_0(0) * 4.6 / grptime * jyconv * 1d9 * omega 	; DN * e/s * (MJy/sr)/(e/s) * 1d9 mJy/MJy * sr
		;flux_1 = mags_1(0) * 4.6 / grptime * jyconv * 1d9 * omega 

;	      	print,  progid, aorkey, date, tag, object, filter, grptime, xcen_0, ycen_0, mags_0, sky_0, flux_0, $
;                 format='(A7, 4X, I8, 3X, A12, 3X,A8, 3X, A14, 3X,A4, 2X, F6.3, 4X, F6.2, 2X, F6.2, 2X, F10.1, 2X, F6.0, F7.2)'
;
;                print,  progid, aorkey, date, tag, object, filter, grptime, xcen_1, ycen_1, mags_1, sky_1, flux_1, $
;                  format='(A7, 4X, I8, 3X, A12, 3X,A8, 3X, A14, 3X,A4, 2X, F6.3, 4X, F6.2, 2X, F6.2, 2X, F10.1, 2X, F6.0, F7.2)'



;	  	print,''
;		print, 'Average flux = ',string(mean([flux_0, flux_1]),format='(f5.2)'),' mJy'

                printf, 5,  progid, aorkey, date, tag, object, filter, grptime, xcen_0, ycen_0, mags_0, sky_0, flux_0, $
                  format='(A7, 4X, I8, 3X, A12, 3X,A8, 3X, A14, 3X,A4, 2X, F6.3, 4X, F6.2, 2X, F6.2, 2X, F10.1, 2X, F6.0, F7.2)'

                printf, 5,  progid, aorkey, date, tag, object, filter, grptime, xcen_1, ycen_1, mags_1, sky_1, flux_1, $
                  format='(A7, 4X, I8, 3X, A12, 3X,A8, 3X, A14, 3X,A4, 2X, F6.3, 4X, F6.2, 2X, F6.2, 2X, F10.1, 2X, F6.0, F7.2)'

                extast, hdr_0, astr_0
                xy2ad, xcen_0, ycen_0, astr_0, a_0, d_0
                radec, a_0, d_0, ihr_0, imin_0, xsec_0, ideg_0, imn_0, xsc_0
                extast, hdr_1, astr_1
                xy2ad, xcen_1, ycen_1, astr_0, a_1, d_1
                radec, a_1, d_1, ihr_1, imin_1, xsec_1, ideg_1, imn_1, xsc_1
                             
		if ideg_0 ge 0 then $
                printf, 6, progid, 'ACQ_'+strn(aorkey), date, tag, object, filter, xcen_0, ycen_0, ihr_0, imin_0, xsec_0, ideg_0, imn_0, xsc_0, $
                  format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X, F6.2, 4X, F6.2, 3X, i02,":", i02, ":", f05.2, 4X, i02,":", i02, ":", f05.2)' $
		else $
                printf, 6, progid, 'ACQ_'+strn(aorkey), date, tag, object, filter, xcen_0, ycen_0, ihr_0, imin_0, xsec_0, ideg_0, imn_0, xsc_0, $
                  format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X, F6.2, 4X, F6.2, 3X, i02,":", i02, ":", f05.2, 4X, i03,":", i02, ":", f05.2)'

		if ideg_1 ge 0 then $
                printf, 6, progid, 'ACQ_'+strn(aorkey), date, tag, object, filter, xcen_1, ycen_1, ihr_1, imin_1, xsec_1, ideg_1, imn_1, xsc_1, $
                  format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X, F6.2, 4X, F6.2, 3X, i02,":", i02, ":", f05.2, 4X, i02,":", i02, ":", f05.2)' $
		else $
                printf, 6, progid, 'ACQ_'+strn(aorkey), date, tag, object, filter, xcen_1, ycen_1, ihr_1, imin_1, xsec_1, ideg_1, imn_1, xsc_1, $
                  format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X, F6.2, 4X, F6.2, 3X, i02,":", i02, ":", f05.2, 4X, i03,":", i02, ":", f05.2)'

; Write out the files after modifying the header to fix the WCS bug of S10.5


        fxaddpar, hdr_0, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_0, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_0, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_0, 'HISTORY', 'Coadded PU acquisition DCEs 0:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_0, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_0, 'HISTORY', 'Created by pu_arch program of V. Charmandaris/K. Willett '

        fxaddpar, hdr_1, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_1, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_1, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_1, 'HISTORY', 'Coadded PU acquisition DCEs 0:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_1, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_1, 'HISTORY', 'Created by pu_arch program of V. Charmandaris/K. Willett '


        writefits, fitsdir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d0.fits', pu_d0_dj, hdr_0
        writefits, fitsdir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d1.fits', pu_d1_dj, hdr_1

;	print,''

; Save the average of the two fluxes to an IDL file

	if finite(flux_0,/nan) eq 1 then pu_flux(i/8) = flux_1 else $
		if finite(flux_1,/nan) eq 1 then pu_flux(i/8) = flux_0 else $
		pu_flux(i/8) = mean([flux_0,flux_1])

	pu_filter(i/8) = filter
	pu_type(i/8) = ispupos
	pufilelist[i/8] = tag
    i=i+7

    if not keyword_set(quiet) then print,tag,filter,ispupos,pu_flux[i/8]
endfor

    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    close, 5


    printf, 6, '-------------------------------------------------------------------------------------------------------'
    close,6

	pu16 = fltarr(n_elements(pu_flux))
	pu22 = fltarr(n_elements(pu_flux))
	pu_flux(where(pu_type eq 1)) = 0				; Eliminate fluxes of peakup stars
	pu_red = where(pu_filter eq 'Red' and pu_flux ne 0)
	pu_blue = where(pu_filter eq 'Blu' and pu_flux ne 0)
	pu16(pu_blue) = pu_flux(pu_blue)
	pu22(pu_red) = pu_flux(pu_red)

	pu16_err = 0.20 * pu16
	pu22_err = 0.20 * pu22

	pu16 = pu16 * 1d-3 & pu22 = pu22 * 1d-3				; Convert from mJy to Jy
	pu16_err = pu16_err * 1d-3 & pu22_err = pu22_err * 1d-3

	; Append zero values for arch012, which does not have a peakup

	pufilelist_temp = strarr(n_elements(pufilelist)+1) & pufilelist_temp[0:n_elements(pufilelist)-1] = pufilelist
	pu_flux_temp = fltarr(n_elements(pu_flux)+1) & pu_flux_temp[0:n_elements(pu_flux)-1] = pu_flux
	pu_filter_temp = pufilelist_temp & pu_filter_temp[0:n_elements(pu_filter)-1] = pu_filter
	pu_type_temp = intarr(n_elements(pu_type)+1) & pu_type_temp[0:n_elements(pu_type)-1] = pu_type
	pu16_temp = pu_flux_temp & pu16_temp[0:n_elements(pu16)-1] = pu16
	pu22_temp = pu_flux_temp & pu22_temp[0:n_elements(pu22)-1] = pu22
	pu16_err_temp = pu_flux_temp & pu16_err_temp[0:n_elements(pu16_err)-1] = pu16_err
	pu22_err_temp = pu_flux_temp & pu22_err_temp[0:n_elements(pu22_err)-1] = pu22_err

	pufilelist = pufilelist_temp
	pu_flux = pu_flux_temp
	pu_filter = pu_filter_temp
	pu_type = pu_type_temp
	pu16 = pu16_temp
	pu22 = pu22_temp
	pu16_err = pu16_err_temp
	pu22_err = pu22_err_temp

	lastind = n_elements(pufilelist)-1
	pufilelist[lastind] = 'arch012'
	pu_filter[lastind] = 'Blu'

    save, file='~/Astronomy/Research/Spitzer/archived/data/idl_sav/pu_fluxes.sav', $
	    pufilelist, pu_flux, pu_filter, pu_type, pu16, pu22, pu16_err, pu22_err

    if keyword_set(stop) then stop
end
