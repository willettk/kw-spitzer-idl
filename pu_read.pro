pro pu_read, campaign, write=write, phot=phot, wphot=wphot, wcs=wcs

;---------------------------------------------------------------------
; Program which reads DCS PU images, combines them doing CR rejection
; and jailbar correction and produces one final image for the
; acquisition and sweet spot location. 
; 
; The final image does have WCS information.
;
; An option (phot) to perform a centroid at the expected position of
; the source and calculate the photometry over two apertures is
; included.
; 
; Options:
;    /write  : Create the combined files
;    /phot   : Do the photometry and display it on the screen
;    /wphot  : Do the photometry and write the output into a file
;    /wcs    : Compute the RA/DEC of the PU source from the WCS
;              of the FITS header and write it into a file
; 
; Example:
;    - See the output on screen
;             pu_read, /write, /phot
;    - Write the output to a file
;             pu_read, /write, /phot, /wphot
; 
; Dependencies: 
;    - dejailbar.pro
; 
; Version: 1.2  (15 July. 2006) - update conversion factors based on
;                                 the numbers calculated in July 2005
; Version: 1.1  (24 Nov. 2004) - update conversion factors & take into
;                                account the various exptimes of PUs
; Version: 1.0  (20 Sep. 2004)
; Author: Vassilis Charmandaris (Cornell Univ.)
;---------------------------------------------------------------------

; You may want to modify the following four. 
;
;files_in='/home/pipeprods/online/'+strn(campaign)+'-12.0/*/bcd/*_000[0,1]_000?_*bcd.fits'
;files_in='/home/pipeprods/online/'+strn(campaign)+'-10.0/*/bcd/*_000[0,1]_000?_*bcd.fits'
;files_in='/home/pipeprods/online/'+strn(campaign)+'-13.2/*/bcd/*_000[0,1]_000?_*bcd.fits'
files_in = '~/Astronomy/Research/Spitzer/archived/data/'+strn(campaign)+'_spectra/*/bcd/*_000[0,1]_000?_*bcd.fits'
;dir_out='/home/data/Peakup/science_new/IMAGES/'
;phot_file='/home/data/Peakup/science_new/pu_phot'+strn(campaign)+'.txt'
;centroid_file='/home/data/Peakup/science_new/pu_cntrds_'+strn(campaign)+'.txt'
dir_out='~/Desktop/'
phot_file='~/Desktop/'+strn(campaign)+'_phot.txt'
centroid_file='~/Desktop/'+strn(campaign)+'_cen.txt'


;files_in='r9*/ch0/bcd/*bcd.fits'
;dir_out='/home/vassilis/IRS_DATA/COMET/S10.5/PU/IMAGES/'
;phot_file='/home/vassilis/IRS_DATA/COMET/S10.5/PU/phot_pu'+strn(campaign)+'.txt'
;centroid_file='/home/vassilis/IRS_DATA/COMET/S10.5/PU/pu_centroids_'+strn(campaign)+'.txt'

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
    openw, 5, phot_file
    printf, 5, '------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'File created on: ',systime()
    printf, 5, '------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'Flux Conversion factors:'
    printf, 5, '16um:', bf, 'DN/Jy, GRPTIME=8.3888sec, 3pix', format='(A5, X, E8.2, X, A30)'
    printf, 5, '22um:', rf, 'DN/Jy, GRPTIME=8.3888sec, 4pix', format='(A5, X, E8.2, X, A30)'
    printf, 5, '------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'PROGID', 'AORKEY', 'OBS DATE', 'TARGET  ', 'FILTR', 'GRPTIME','X-CNTR', 'Y-CNTR', '[3/4pix]DN','[5/7pix]DN', 'SKY_DN', 'f(mJy)',$
      format='(A6, 3X,A8, 3X, A12, 3X,A12, 1X,A6, 2X, A7, 2X, A6, 2X, A6, 2X, A10,2X, A10, 2X, A6, 2X, A6)'
    printf, 5, '------------------------------------------------------------------------------------------------------------------------'
endif

if keyword_set(wcs) then begin
    openw, 6, centroid_file
    printf, 6, 'PROGID', 'AORKEY', 'OBS DATE', 'TARGET  ', 'FILTER', 'X-CNTR', 'Y-CNTR', 'RA (J2000)', 'Dec(J2000)', $
      format='(A6, 3X, A12, 3X, A12, 3X,A12, 1X,A6, 2X, A6, 2X, A6, 2X, A12,2X, A12)'
    printf, 6, '-------------------------------------------------------------------------------------------------------'
endif

for i=0, n_targets-1 do begin
    temp=readfits(files(i),hdr, /silent)
    readmode=sxpar(hdr, 'READMODE')
    if strn(readmode) ne 'PEAKUP' then begin
        print, i, 'is not a DCS PeakUp', strn(readmode), 'skip...', format='(I4,X,A20,X,A4,X,A7)'
        goto, JUMP1
    endif
    print, i, 'is a DCS PeakUp', format='(I4,X,A20)'
    a=readfits(files(i),hdr_acq, /silent)
    b=readfits(files(i+4),hdr_ss, /silent)
;    print, files(i), files(i+4)
;    wait, 1
; Get the AORKEY
    aorkey=sxpar(hdr_acq, 'AORKEY')
; Get the Program ID
    progid='P_'+strn(sxpar(hdr_acq, 'PROGID'))
; Get the target name
    object_long=sxpar(hdr_acq, 'OBJECT')
; Get the filter used
    filter=sxpar(hdr_acq, 'PKUPFILT')
; Get the exposure times 
    time_acq=sxpar(hdr_acq, 'GRPTIME')
    time_ss=sxpar(hdr_ss, 'GRPTIME')
; Log the day of observations
    julianday=sxpar(hdr_acq, 'MJD_OBS')+2400000.5
    daycnv, julianday, year,month,day
    date=strn(YMD2DATE(year,month,day))

; Check if on source PU  
    pu_status=sxpar(hdr_acq, 'ISPUPOS')
    if pu_status then begin
        object='OFF_'+STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),''),'?',/extract),'')
        nameout=progid+'_'+STRMID(object,0,15)+'_'+STRTRIM(filter,2)
    endif else begin
;        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),'')
        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),''), '?',/extract),'')
        nameout=progid+'_'+STRMID(object,0,15)+'_'+STRTRIM(filter,2)
    endelse

; Check if there's WCS in the header
    is_wcs_in=sxpar(hdr_acq, 'CD1_1')

    if keyword_set(write) then begin
    
; read acquisition images

        acq0=readfits(files(i), /silent)
        acq1=readfits(files(i+1), /silent)
        acq2=readfits(files(i+2), /silent)
	print, 'ACQ Images:'
	print, files(i)
	print, files(i+1)
	print, files(i+2)
	wait, 1
; read sweet spot images
        ss0=readfits(files(i+4), /silent)
        ss1=readfits(files(i+5), /silent)
        ss2=readfits(files(i+6), /silent)
   	print, 'SS Images:'
	print, files(i+4)
	print, files(i+5)
	print, files(i+6)
	wait,1 
; cosmic ray rejection
        pu_acq_clean=(acq0 < acq1) < acq2
        pu_ss_clean=(ss0 < ss1) < ss2
        
        dejailbar, pu_acq_clean, pu_acq_dj
        dejailbar, pu_ss_clean, pu_ss_dj

        if keyword_set(phot) then begin

; photometry
            if strn(filter) eq 'RED' then begin
                pos_x=105
                pos_y=92
                aperture=[4,7]
                dn_to_mJy_acq=(rf/1000.) * (time_acq/8.3888)
                dn_to_mJy_ss=(rf/1000.) * (time_ss/8.3888)
;        aperture=[5]
            endif else begin
                pos_x=107
                pos_y=30
                aperture=[3,5]
                dn_to_mJy_acq=(bf/1000.)* (time_acq/8.3888)
                dn_to_mJy_ss =(bf/1000.) * (time_ss/8.3888)
;        aperture=[3]
            endelse

; Sky annulus [8 - 14]
            skyrad=[8,14]
            badpix=[-32765,65000]
            phpadu=1
;  Centroid is computed using a box of half width equal to 
; 1.5 sigma = 0.637* FWHM
            fwhm=2.
        
; Photometry of AQC image
            cntrd, pu_acq_clean, pos_x, pos_y, xcen, ycen, fwhm
            aper, pu_acq_clean, xcen, ycen, mags, errap, sky, skyerr, phpadu, aperture, skyrad, badpix, $
              /flux, /silent
            print, progid, aorkey, date, object, filter, time_acq, xcen, ycen, mags(0), mags(1), sky, mags(0)/dn_to_mjy_acq, $
              format='(A6, 3X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            if keyword_set(wphot) then begin
                printf, 5,  progid, aorkey, date, object, filter, time_acq, xcen, ycen, mags(0), mags(1), sky, mags(0)/dn_to_mjy_acq, $
                  format='(A6, 3X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'
            endif

            if keyword_set(wcs) and is_wcs_in ne 0 then begin
                extast, hdr_acq, astr
                xy2ad, xcen, ycen, astr, a, d
                radec, a, d, ihr, imin, xsec, ideg, imn, xsc
                
                printf, 6, progid, 'ACQ_'+strn(aorkey), date, object, filter, xcen, ycen, ihr, imin, xsec, ideg, imn, xsc, $
                  format='(A6, 3X, A12, 3X, A12, 3X,A12, 3X,A4, 2X, F6.2, 2X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'
            endif 

; Photometry of SS image
            cntrd, pu_ss_clean, pos_x, pos_y, xcen, ycen, fwhm
            aper, pu_ss_clean, xcen, ycen, mags, errap, sky, skyerr, phpadu, aperture, skyrad, badpix, $
              /flux, /silent
            print, progid, aorkey, date, object, filter, time_ss, xcen, ycen, mags(0), mags(1), sky, mags(0)/dn_to_mjy_ss,$
              format='(A6, 3X, I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'

            if keyword_set(wphot) then begin
                printf, 5, progid, aorkey, date, object, filter, time_ss, xcen, ycen, mags(0), mags(1), sky, mags(0)/dn_to_mjy_ss, $
                  format='(A6, 3X,I8, 3X, A12, 3X,A12, 3X,A4, 2X, F6.3, 2X, F6.2, 2X, F6.2, 2X, F10.1, 2X,F10.1, 2X, F6.0, F7.1)'
            endif

            if keyword_set(wcs)  and is_wcs_in ne 0 then begin
                extast, hdr_ss, astr
                xy2ad, xcen, ycen, astr, a, d
                radec, a, d, ihr, imin, xsec, ideg, imn, xsc
                printf, 6, progid, 'SS_'+strn(aorkey), date, object, filter, xcen, ycen, ihr, imin, xsec, ideg, imn, xsc, $
                  format='(A6, 3X, A12, 3X, A12, 3X,A12, 3X,A4, 2X, F6.2, 2X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'
            endif 
        endif


; Write out the files after modifying the header to fix the WCS bug of S10.5

        fxaddpar, hdr_acq, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_acq, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_acq, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_acq, 'HISTORY', 'Coadded PU acquisition DCEs 0:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_acq, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_acq, 'HISTORY', 'Created by pu_read program of V. Charmandaris '

        fxaddpar, hdr_ss, 'CTYPE1','RA---TAN-SIP', 'RA---TAN with distortion in pixel space'
        fxaddpar, hdr_ss, 'CTYPE2','DEC--TAN-SIP', 'DEC--TAN with distortion in pixel space'
        fxaddpar, hdr_ss, 'HISTORY', 'Added the CTYPE keywords in header for WCS to work'
        fxaddpar, hdr_ss, 'HISTORY', 'Coadded PU sweet spot DCEs 1:0,1,2 on '+SYSTIME()
        fxaddpar, hdr_ss, 'HISTORY', 'Cosmic ray / jailbar removal on '+SYSTIME()
        fxaddpar, hdr_ss, 'HISTORY', 'Created by pu_read program of V. Charmandaris'

        print, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_acq.fits'
        writefits, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_acq.fits', pu_acq_dj, hdr_acq
        writefits, dir_out+nameout+'_'+strn(aorkey)+'_'+date+'_ss.fits', pu_ss_dj, hdr_ss

        print, 'Wrote ', nameout
    endif

    i=i+7
    JUMP1: 
endfor

if keyword_set(wphot) then begin
    printf, 5, '------------------------------------------------------------------------------------------------------------------------'
    close, 5
endif


if keyword_set(wcs) then begin
    printf, 6, '-------------------------------------------------------------------------------------------------------'
    close,6
endif

end
