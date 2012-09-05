pro pu_vas_sur, write=write, phot=phot, wphot=wphot, wcs=wcs
;+
; NAME: 
;       PU_DCS 
;
; PURPOSE:
; 	Read Spitzer DCS peakup images, combine them doing cosmic ray rejection,
;	jailbar correct, and produce final image for both acquisition and
;	sweet spot location.
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
;
; CALLING SEQUENCE:
; 	pu_dcs [, /write, /wphot, /phot, /wcs]
;
; INPUTS:
;	Spitzer DCS peakup FITS images with header information
;	
; OUTPUTS:
;	- two mosaicked FITS files per target with cleaned image
;	- pu_phot.txt: file containing photometry, centroid, flux values 
;	- pu_centroid.txt: file containing centroid, RA and DEC values
;
; KEYWORDS:
;
;    	/WRITE  : 	Create the combined files
;    	/PHOT   : 	Do the photometry and display it on the screen
;    	/WPHOT  : 	Do the photometry and write the output into a file
;    	/WCS    : 	Compute the RA/DEC of the PU source from the WCS
;              			of the FITS header and write it into a file
;
; REQUIRES:
;	DEJAILBAR.PRO
;
; EXAMPLE:
;
;    - See the output on screen
;             IDL> pu_read, /write, /phot
;    - Write the output to a file
;             IDL> pu_read, /write, /phot, /wphot
;         
; MODIFICATION HISTORY:
;	written by V. Charmandaris (originally PU_READ.pro):
;
; 	Version: 1.2  (15 July. 2006) - update conversion factors based on
;                                 the numbers calculated in July 2005
; 	Version: 1.1  (24 Nov. 2004) - update conversion factors & take into
;                                account the various exptimes of PUs
; 	Version: 1.0  (20 Sep. 2004)
;
;	Written to test Vassilis' method of data reduction on the SUR images
;-

write = 1
wcs = 1
phot = 1
wphot = 1

; Directories for OHMs in which to search for DCS peakup files

files_in='~/Astronomy/Research/Spitzer/OHM/data/mega0??_peakup/ch0/bcd/*_000[0,1,2]_000?_*bcdb.fits'
dir_out='~/Astronomy/Research/Spitzer/OHM/data/pu_vas_sur/'
phot_file='~/Astronomy/Research/Spitzer/OHM/data/pu_vas_sur/pu_phot.txt'
centroid_file='~/Astronomy/Research/Spitzer/OHM/data/pu_vas_sur/pu_cntrds.txt'

; Directories for CSO DCS peakups

;files_in='~/Astronomy/Research/Spitzer/cso/data/cso0??_spectra/ch0/bcd/*_000[0,1]_000?_*bcd.fits'
;dir_out='~/Astronomy/Research/Spitzer/cso/data/pu_read/'
;phot_file='~/Astronomy/Research/Spitzer/cso/data/pu_read/pu_phot.txt'
;centroid_file='~/Astronomy/Research/Spitzer/cso/data/pu_read/pu_cntrds.txt'

print, 'Template files: ',files_in

files=file_search(files_in)
num= n_elements(files)

; Conversion factors
bf=1.18e6                       ; BLUE PU: [DN in 8.3888sec] /Jy : 3pix aperture
rf=1.15e6                       ; RED PU:  [DN in 8.3888sec] /Jy : 4pix aperture

; Up to Jun. 20, 2005
; bf=1.33e6  ; BLUE PU: DN/Jy 3pix aperture
; rf=1.13e6  ; RED PU: DN/Jy 4pix aperture
; Up to Sep. 20, 2004
; bf=1.34e6  ; BLUE PU: DN/Jy 3pix aperture
; rf=1.16e6  ; RED PU: DN/Jy 4pix aperture

;n_targets = num/8
n_targets = num

print, 'Images found: ',n_targets

if keyword_set(wphot) then begin
    openw, 5, phot_file, /append
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'File created on: ',systime()
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'Flux Conversion factors:'
    printf, 5, '16um:', bf, 'DN/Jy, GRPTIME=8.3888sec, 3pix', format='(A5, X, E8.2, X, A30)'
    printf, 5, '22um:', rf, 'DN/Jy, GRPTIME=8.3888sec, 4pix', format='(A5, X, E8.2, X, A30)'
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'PROGID', 'AORKEY', 'OBS DATE', 'TARGET  ', 'FILTR', 'GRPTIME','X-CNTR', 'Y-CNTR', '[3/4pix]DN','[5/7pix]DN', 'SKY_DN', 'f(mJy)',$
      format='(A7, 4X, A8, 3X, A12, 3X,A12, 1X,A6, 2X, A7, 2X, A6, 2X, A6, 2X, A10,2X, A10, 2X, A6, 2X, A6)'
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
endif

if keyword_set(wcs) then begin
    openw, 6, centroid_file, /append
    printf, 6, '------------------------------------------------------------------------------------------------------------------------'
    printf, 6, 'File created on: ',systime()
    printf, 6, '------------------------------------------------------------------------------------------------------------------------'
    printf, 6, 'PROGID', 'AORKEY', 'OBS DATE', 'TARGET  ', 'FILTER', 'X-CNTR', 'Y-CNTR', 'RA (J2000)', 'Dec(J2000)', $
      format='(A7, 4X, A12, 3X, A12, 3X,A12, 1X,A6, 2X, A6, 2X, A6, 2X, A12,2X, A12)'
    printf, 6, '--------------------------------------------------------------------------------------------------------------'
endif

for i=0, n_targets-1 do begin
    temp=readfits(files(i),hdr, /silent)
    readmode=sxpar(hdr, 'READMODE')
;    if strn(readmode) ne 'PEAKUP' then begin
;        print, i, 'is not a DCS PeakUp', strn(readmode), 'skip...', format='(I4,X,A20,X,A4,X,A7)'
;        goto, JUMP1
;    endif
;    print, i, 'is a DCS PeakUp', format='(I4,X,A20)'
    a=readfits(files(i),hdr_0, /silent)
    b=readfits(files(i+2),hdr_1, /silent)
    c=readfits(files(i+4),hdr_2, /silent)
;    print, files(i), files(i+4)
;    wait, 1
; Get the AORKEY
    aorkey=sxpar(hdr_0, 'AORKEY')
; Get the Program ID
    progid='P_'+strn(sxpar(hdr_0, 'PROGID'))
; Get the target name
    object_long=sxpar(hdr_0, 'OBJECT')
; Get the filter used
;   filter=sxpar(hdr_0, 'PKUPFILT')
    fovname=sxpar(hdr_0, 'FOVNAME')
    filter = STRMID(fovname,4,3)
; Get the exposure times 
    time_acq=sxpar(hdr_0, 'GRPTIME')
; Log the day of observations
    julianday=sxpar(hdr_0, 'MJD_OBS')+2400000.5
    daycnv, julianday, year,month,day
    date=strn(YMD2DATE(year,month,day))

; Check if on source PU  
;    pu_status=sxpar(hdr_acq, 'ISPUPOS')
;    if pu_status then begin
;        object='OFF_'+STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),''),'?',/extract),'')
;        nameout=progid+'_'+STRMID(object,0,18)+'_'+STRTRIM(filter,2)+'_'+STRMID(files(i), 52, 7)
;    endif else begin
;        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),'')
        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),''), '?',/extract),'')
        nameout=progid+'_'+STRMID(object,0,18)+'_'+STRTRIM(filter,2)+'_'+STRMID(files(i), 52, 7)
;    endelse

; Check if there's WCS in the header
    is_wcs_in=sxpar(hdr_0, 'CD1_1')

    if keyword_set(write) then begin
    
; read images from dither 0

        d00=readfits(files(i), /silent)
        d01=readfits(files(i+1), /silent)
	print, 'Dither pos. 0 Images:'
	print, files(i)
	print, files(i+1)
;	wait, 1

; read images from dither 1

        d10=readfits(files(i+2), /silent)
        d11=readfits(files(i+3), /silent)
	print, 'Dither pos. 1 Images:'
	print, 'ACQ Images:'
	print, files(i+2)
	print, files(i+3)
;	wait, 1

; read images from dither 2

        d20=readfits(files(i+4), /silent)
        d21=readfits(files(i+5), /silent)
	print, 'Dither pos. 2 Images:'
	print, files(i+4)
	print, files(i+5)
;	wait, 1

; cosmic ray rejection
        pu_d0_clean= d00 < d01
        pu_d1_clean= d10 < d11
        pu_d2_clean= d20 < d21
        
        if keyword_set(phot) then begin

; photometry
            if strn(filter) eq 'Red' then begin
                pos_x_0=22
                pos_y_0=28

                pos_x_1=20
                pos_y_1=30
                
		pos_x_2=19
                pos_y_2=29

                aperture=[4,7]
                dn_to_mJy_acq=(rf/1000.) * (time_acq/8.3888) /4.	; Adjust for new grptime

            endif else begin
                pos_x_0=25
                pos_y_0=30

                pos_x_1=20
                pos_y_1=29
                
		pos_x_2=26
                pos_y_2=28

                aperture=[3,5]
                dn_to_mJy_acq=(bf/1000.)* (time_acq/8.3888 ) /4.	; Adjust for new grptime

            endelse

; Sky annulus [8 - 14]
            skyrad=[8,14]
            badpix=[-32765,65000]
            phpadu=1
;  Centroid is computed using a box of half width equal to 
; 1.5 sigma = 0.637* FWHM
            fwhm=2.

; Photometry 
            cntrd, pu_d0_clean, pos_x_0, pos_y_0, xcen_0, ycen_0, fwhm
            aper, pu_d0_clean, xcen_0, ycen_0, mags_0, errap_0, sky_0, skyerr_0, phpadu, aperture, skyrad, badpix, $
              /flux, /silent

            cntrd, pu_d1_clean, pos_x_1, pos_y_1, xcen_1, ycen_1, fwhm
            aper, pu_d1_clean, xcen_1, ycen_1, mags_1, errap_1, sky_1, skyerr_1, phpadu, aperture, skyrad, badpix, $
              /flux, /silent

            cntrd, pu_d2_clean, pos_x_2, pos_y_2, xcen_2, ycen_2, fwhm
            aper, pu_d2_clean, xcen_2, ycen_2, mags_2, errap_2, sky_2, skyerr_2, phpadu, aperture, skyrad, badpix, $
              /flux, /silent

            print, progid, aorkey, date, object, filter, time_acq, xcen_0, ycen_0, mags_0(0), mags_0(1), sky_0, mags_0(0)/dn_to_mjy_acq, $
              format='(A7, 4X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            print, progid, aorkey, date, object, filter, time_acq, xcen_1, ycen_1, mags_1(0), mags_1(1), sky_1, mags_1(0)/dn_to_mjy_acq, $
              format='(A7, 4X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            print, progid, aorkey, date, object, filter, time_acq, xcen_2, ycen_2, mags_2(0), mags_2(1), sky_2, mags_2(0)/dn_to_mjy_acq, $
              format='(A7, 4X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            if keyword_set(wphot) then begin

            printf, 5, progid, aorkey, date, object, filter, time_acq, xcen_0, ycen_0, mags_0(0), mags_0(1), sky_0, mags_0(0)/dn_to_mjy_acq, $
              format='(A7, 4X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            printf, 5, progid, aorkey, date, object, filter, time_acq, xcen_1, ycen_1, mags_1(0), mags_1(1), sky_1, mags_1(0)/dn_to_mjy_acq, $
              format='(A7, 4X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            printf, 5, progid, aorkey, date, object, filter, time_acq, xcen_2, ycen_2, mags_2(0), mags_2(1), sky_2, mags_2(0)/dn_to_mjy_acq, $
              format='(A7, 4X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            endif

            if keyword_set(wcs) and is_wcs_in ne 0 then begin
                extast, hdr_0, astr_0
                xy2ad, xcen_0, ycen_0, astr_0, a_0, d_0
                radec, a_0, d_0, ihr_0, imin_0, xsec_0, ideg_0, imn_0, xsc_0
                
                extast, hdr_1, astr_1
                xy2ad, xcen_1, ycen_1, astr_1, a_1, d_1
                radec, a_1, d_1, ihr_1, imin_1, xsec_1, ideg_1, imn_1, xsc_1
                
                extast, hdr_2, astr_2
                xy2ad, xcen_2, ycen_2, astr_2, a_2, d_2
                radec, a_2, d_2, ihr_2, imin_2, xsec_2, ideg_2, imn_2, xsc_2
                
                printf, 6, progid, 'D0_'+strn(aorkey), date, object, filter, xcen_0, ycen_0, ihr_0, imin_0, xsec_0, ideg_0, imn_0, xsc_0, $
                  format='(A7, 4X, A11, 3X, A12, 3X,A12, 3X,A4, 2X, F6.2, 2X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'

                printf, 6, progid, 'D1_'+strn(aorkey), date, object, filter, xcen_1, ycen_1, ihr_1, imin_1, xsec_1, ideg_1, imn_1, xsc_1, $
                  format='(A7, 4X, A11, 3X, A12, 3X,A12, 3X,A4, 2X, F6.2, 2X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'

                printf, 6, progid, 'D0_'+strn(aorkey), date, object, filter, xcen_2, ycen_2, ihr_2, imin_2, xsec_2, ideg_2, imn_2, xsc_2, $
                  format='(A7, 4X, A11, 3X, A12, 3X,A12, 3X,A4, 2X, F6.2, 2X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'

            endif 

    endif

; Write out the files after modifying the header to fix the WCS bug of S10.5

        fxaddpar, hdr_0, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_0, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_0, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_0, 'HISTORY', 'Coadded PU acquisition DCEs 0:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_0, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_0, 'HISTORY', 'Created by pu_read program of V. Charmandaris/K. Willett '

        fxaddpar, hdr_1, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_1, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_1, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_1, 'HISTORY', 'Coadded PU acquisition DCEs 0:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_1, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_1, 'HISTORY', 'Created by pu_read program of V. Charmandaris/K. Willett '

        fxaddpar, hdr_2, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_2, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_2, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_2, 'HISTORY', 'Coadded PU acquisition DCEs 0:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_2, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_2, 'HISTORY', 'Created by pu_read program of V. Charmandaris/K. Willett '

        print, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'.fits'
        writefits, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d0.fits', pu_d0_clean, hdr_0
        writefits, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d1.fits', pu_d1_clean, hdr_1
        writefits, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d2.fits', pu_d2_clean, hdr_2

	d0name = dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d0.fits'
	d1name = dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d1.fits'
	d2name = dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d2.fits'

        print, 'Wrote ', nameout
    endif

    i=i+5
    JUMP1: 
endfor

if keyword_set(wphot) then begin
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    close, 5
endif


if keyword_set(wcs) then begin
    printf, 6, '-------------------------------------------------------------------------------------------------------'
    close,6
endif

stop
end
