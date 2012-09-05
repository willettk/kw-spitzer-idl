pro pu_blue_cso
;+
; NAME: 
;       PU_BLUE_CSO
;
; PURPOSE:
;
;	Extract the blue offset peakup fluxes from a subset of targets in the OHM and CSO samples
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
;
; CALLING SEQUENCE:
; 	pu_blue_cso
;
; INPUTS:
;	
; OUTPUTS:
;
; KEYWORDS:
;
; REQUIRES:
;
; EXAMPLE:
;
;	IDL> pu_blue_cso
;         
; MODIFICATION HISTORY:
;
;	Written by KW (adapted from PU_READ.pro by V. Charmandaris) - Aug 07
;-

; Program is run twice - once for blue and once for red peakups

; List of CSOs with spectra using offset peakups

cso_offset = ['cso003','cso006']

; Note - none of the second batch (7-10) have offset peakups.

filelist=['~/Astronomy/Research/Spitzer/CSO/data/cso0??_spectra/ch0/bcd/*_000[0,1]_000[0,1,2,3]_*bcd.fits']

files=file_search(filelist)


; Cull list

newfiles = strarr(n_elements(files))
m = 0
for k = 0, n_elements(files) - 1 do begin
	temptag = strmid(files(k), 52, 6)
	tempindex = where(temptag eq cso_offset)
	if tempindex ge 0 then begin 
		newfiles(m) = files(k)
		m = m+1
	endif
endfor
zerofiles = where(newfiles eq '')
newfiles = newfiles(0:zerofiles(0)-1)
n_targets = n_elements(newfiles)
print, 'Targets found: ',n_targets / 8

dir_out='~/Astronomy/Research/Spitzer/CSO/pu_blue/'
fitsdir_out='~/Astronomy/Research/Spitzer/CSO/pu_blue/fits/'
phot_file='~/Astronomy/Research/Spitzer/CSO/pu_blue/pu_blue_phot.txt'
centroid_file='~/Astronomy/Research/Spitzer/CSO/pu_blue/pu_blue_cntrds.txt'

; Conversion factors
bf=1.18e6                       ; BLUE PU: [DN in 8.3888sec] /Jy : 3pix aperture
rf=1.15e6                       ; RED PU:  [DN in 8.3888sec] /Jy : 4pix aperture

; Up to Jun. 20, 2005
; bf=1.33e6  ; BLUE PU: DN/Jy 3pix aperture
; rf=1.13e6  ; RED PU: DN/Jy 4pix aperture
; Up to Sep. 20, 2004
; bf=1.34e6  ; BLUE PU: DN/Jy 3pix aperture
; rf=1.16e6  ; RED PU: DN/Jy 4pix aperture

    openw, 5, phot_file
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'File created on: ',systime()
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG ', 'TARGET  ', 'FILTR', 'GRPTIME', 'X-CNTR', 'Y-CNTR', '3 px [DN]', 'Sky [DN]', 'f(mJy)',$
      format='(A7, 4X, A8, 3X, A12, 3X, A8, 3X, A14, 1X,A6, 2X, A7, 2X, A7, 2X, A6, 2X, A9, 2X, A9, 2X, A6)'
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'

   
    openw, 6, centroid_file
    printf, 6, '------------------------------------------------------------------------------------------------------------------------'
    printf, 6, 'File created on: ',systime()
    printf, 6, '------------------------------------------------------------------------------------------------------------------------'
    printf, 6, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG ', 'TARGET  ', 'FILTER', 'GRPTIME', 'X-CNTR', 'Y-CNTR', 'RA (J2000)', 'Dec(J2000)', $
      format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A12, 1X,A6, 2X, A7, 2X, A6, 2X, A6, 2X, A12,2X, A12)'
    printf, 6, '--------------------------------------------------------------------------------------------------------------'

  

    print, '-------------------------------------------------------------------------------------------------------------------------------'
    print, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG ', 'TARGET  ', 'FILTR', 'GRPTIME', 'X-CNTR', 'Y-CNTR', '3 px [DN]', 'Sky [DN]', 'f(mJy)',$
      format='(A7, 4X, A8, 3X, A12, 3X, A8, 3X, A14, 1X,A6, 2X, A7, 2X, A7, 2X, A6, 2X, A9, 2X, A9, 2X, A6)'
    print, '-------------------------------------------------------------------------------------------------------------------------------'
 
   
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

    tag = strmid(newfiles(i),52,6)

        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,'[', /extract),''),']', /extract),''),'/', /extract),'_'),/extract),''), '?',/extract),'')

	nameout=progid+'_'+STRMID(object,0,18)+'_'+STRTRIM(filter,2)+'_'+STRMID(newfiles(i), 52, 6)



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

                pos_x_0=108
                pos_y_0=31

                pos_x_1=109
                pos_y_1=29
                
                aperture=[3]
		jyconv = 0.01617d

                dn_to_mJy_acq=(bf/1000.)* (grptime/8.3888)
                dn_to_mJy_ss =(bf/1000.) * (grptime/8.3888)
                ;dn_to_mJy_acq=(bf/1000.)* (time_acq/8.3888)
                ;dn_to_mJy_ss =(bf/1000.) * (time_ss/8.3888)

; Sky annulus [8 - 14]
            skyrad=[12,15]
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

		flux_0 = mags_0 / dn_to_mJy_acq
		flux_1 = mags_1 / dn_to_mJy_acq

		;flux_0 = mags_0(0) * 4.6 / grptime * jyconv * 1d9 * omega 	; DN * e/s * (MJy/sr)/(e/s) * 1d9 mJy/MJy * sr
		;flux_1 = mags_1(0) * 4.6 / grptime * jyconv * 1d9 * omega 

	      	print,  progid, aorkey, date, tag, object, filter, grptime, xcen_0, ycen_0, mags_0, sky_0, flux_0, $
                  format='(A7, 4X, I8, 3X, A12, 3X,A8, 3X, A14, 3X,A4, 2X, F6.3, 4X, F6.2, 2X, F6.2, 2X, F10.1, 2X, F6.0, F7.2)'

                print,  progid, aorkey, date, tag, object, filter, grptime, xcen_1, ycen_1, mags_1, sky_1, flux_1, $
                  format='(A7, 4X, I8, 3X, A12, 3X,A8, 3X, A14, 3X,A4, 2X, F6.3, 4X, F6.2, 2X, F6.2, 2X, F10.1, 2X, F6.0, F7.2)'



	  	print,''
		print, 'Average flux = ',string(mean([flux_0, flux_1]),format='(f7.2)'),' mJy'

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
                             
                printf, 6, progid, 'ACQ_'+strn(aorkey), date, tag, object, filter, xcen_0, ycen_0, ihr_0, imin_0, xsec_0, ideg_0, imn_0, xsc_0, $
                  format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X, F6.2, 4X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'

                printf, 6, progid, 'ACQ_'+strn(aorkey), date, tag, object, filter, xcen_1, ycen_1, ihr_1, imin_1, xsec_1, ideg_1, imn_1, xsc_1, $
                  format='(A7, 4X, A12, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X, F6.2, 4X, F6.2, 2X, i3,":", i2, ":", f5.2, 3X, i3,":", i2, ":", f5.2)'

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


        writefits, fitsdir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d0.fits', pu_d0_dj, hdr_0
        writefits, fitsdir_out+nameout+'_'+strn(aorkey)+'_'+date+'_d1.fits', pu_d1_dj, hdr_1

	print,''

    i=i+7
endfor


    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    close, 5


    printf, 6, '-------------------------------------------------------------------------------------------------------'
    close,6

;    stop
end


