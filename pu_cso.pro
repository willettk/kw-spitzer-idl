pro pu_cso, stop=stop
;+
; NAME: 
;       PU_CSO 
;
; PURPOSE:
; 	Read Spitzer SUR peakup images and perform photometry on PBCD mosaicked images in CSO sample. 
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
;
; CALLING SEQUENCE:
; 	pu_cso 
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
;	IDL> .r pu_cso
;
; MODIFICATION HISTORY:
;
;	Adapted from PU_OHM.pro - KW, May 08
;-

; Program is run twice - once for blue and once for red peakups

filelist=['~/Astronomy/Research/Spitzer/CSO/data/cso0??_peakup/ch0/pbcd/*b_mos.fits', $
	'~/Astronomy/Research/Spitzer/CSO/data/cso0??_peakup/ch0/pbcd/*r_mos.fits']
filterlist = ['Blue','Red']

dir_out='~/Astronomy/Research/Spitzer/CSO/pu_cso/'
phot_file='~/Astronomy/Research/Spitzer/CSO/pu_cso/pu_phot.txt'

; Create photometry file

    openw, 5, phot_file
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'File created on: ',systime()
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
    printf, 5, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG    ', 'TARGET  ', 'FILTR', 'X-CNTR', 'Y-CNTR', '5 px [MJy/sr]', 'Sky [MJy/sr]', 'f(mJy)',$
      format='(A7, 4X, A8, 3X, A12, 3X, A8, 3X, A14, 1X,A6, 2X, A7, 2X, A6, 2X, A13, 2X, A12, 2X, A6)'
    printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'

    print, '-------------------------------------------------------------------------------------------------------------------------------'
    print, 'PROGID', 'AORKEY', 'OBS DATE', 'TAG    ', 'TARGET  ', 'FILTR', 'X-CNTR', 'Y-CNTR', '5 px [MJy/sr]', 'Sky [MJy/sr]', 'f(mJy)',$
      format='(A7, 4X, A8, 3X, A12, 3X, A8, 3X, A14, 1X,A6, 2X, A7, 2X, A6, 2X, A13, 2X, A12, 2X, A6)'
    print, '-------------------------------------------------------------------------------------------------------------------------------'


; Create empty arrays

total_files = [file_search(filelist(0)),file_search(filelist(1))]
n_tot = n_elements(total_files)/2

	pu_tag = strarr(n_tot)
	pu_fluxarr = fltarr(2,n_tot)

for j = 0, 1 do begin

	; Directories for CSOs in which to search for SUR peakup files
	
	files_in = filelist[j]
	
	print,''
	print, 'Template files: ',files_in
	
	files=file_search(files_in)
	n_targets= n_elements(files)
	
	print, 'Images found: ',n_targets
		
	; BEGIN LOOPING OVER IMAGES

	for i=0, n_targets-1 do begin
	
		    mos_img=readfits(files(i),hdr_0, /silent)
		    tag = strmid(files(i), 52, 6)
	
		; Get the AORKEY
		    aorkey=sxpar(hdr_0, 'AORKEY')
	
		; Get the Program ID
		    progid='P_'+strn(sxpar(hdr_0, 'PROGID'))
	
		; Get the target name
		    object_long=sxpar(hdr_0, 'OBJECT')
	
		; Get the units
		    bunit=sxpar(hdr_0, 'BUNIT')
		    if strmid(bunit,0,6) ne 'MJy/sr' then begin
			    print, i, ' does not have units of MJy/sr'
			    stop
	    	    endif
	
		; Get the filter used
		    fovname=sxpar(hdr_0, 'FOVNAME')
		   ;filter = STRMID(fovname,4,3)
		   filter = filterlist[j]
		
		; Log the day of observations
		    julianday=sxpar(hdr_0, 'MJD_OBS')+2400000.5
		    daycnv, julianday, year,month,day
		    date=strn(YMD2DATE(year,month,day))
		
	        object=STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(STRJOIN(STRSPLIT(object_long,$
			'[', /extract),''),']', /extract),''),'/', 			/extract),'_'),/extract),''), '?',/extract),'')
	
		nameout=progid+'_'+STRMID(object,0,18)+'_'+STRTRIM(filter,2)+'_'+STRMID(files(i), 52, 6)
	
		; Check if there's WCS in the header

		    is_wcs_in=sxpar(hdr_0, 'CD1_1')
		    
		; Eyeball glance seems to indicate that the image has been dejailbar-ed already
		; Look at image processing in SOM, or think of a quantitative test I could do 
		; to check on the bias levels (mean value of different sky pixel channels?)
		
		;        dejailbar, pu_clean, pu_dj
		
		
		; Photometry
		
				; Starting positions of target in FOV for both filters
		
				pos_x_0 = 32
				pos_y_0 = 92
		
		            if string(filter) eq 'Red' then begin
		
		                aperture=[4]
				ap_correction = 1.412d
		
		            endif else begin
		
		                aperture=[3]
				ap_correction = 1.418d
		
			    endelse
		
		; Sky annulus [8 - 14]
		            skyrad=[12,15]
		            badpix=[-32765,65000]
		            phpadu=1
		;  Centroid is computed using a box of half width equal to 
		; 1.5 sigma = 0.637* FWHM
		            fwhm=2.
		
		; Photometry of ACQ image
	
	            cntrd, mos_img, pos_x_0, pos_y_0, xcen_0, ycen_0, fwhm
	            aper, mos_img, xcen_0, ycen_0, mags_0, errap_0, sky_0, skyerr_0, phpadu, aperture, skyrad, badpix, $
	              /flux;, /silent
	
	      	    if finite(mags_0) eq 1 then begin
				
	
			; Convert mags_0 (MJy/sr) to mJy
	
			mJy_per_MJy = 1d9
			arcsec2_per_pixel = 3.24
			sr_per_arcsec2 = 2.3504d-11
			sr_per_pixel = sr_per_arcsec2 * arcsec2_per_pixel 
			flux_mJy = mags_0 * sr_per_pixel * mJy_per_MJy * ap_correction
	
	                print,  progid, aorkey, date, tag, object, filter, xcen_0, ycen_0, mags_0, sky_0, flux_mJy, $
	                  format='(A7, 4X, I8, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X,  F6.2, 2X, F6.2, 2X, F13.1, 2X, F6.1, 7X, F7.2)'
	
			; Write to file
	
	                printf, 5,  progid, aorkey, date, tag, object, filter, xcen_0, ycen_0, mags_0, sky_0, flux_mJy, $
	                  format='(A7, 4X, I8, 3X, A12, 3X, A8, 3X, A14, 3X,A4, 2X,  F6.2, 2X, F6.2, 2X, F13.1, 2X, F6.1, 7X, F7.2)'
	
	              	printf,5,''
	
			pu_fluxarr(j,i) = flux_mJy
		    endif 
	
	
		pu_tag(i) = tag

	endfor

endfor

printf, 5, '-------------------------------------------------------------------------------------------------------------------------------'
close, 5

print,''

pu16 = transpose(pu_fluxarr(0,*))
pu22 = transpose(pu_fluxarr(1,*))

pu16_err = 0.20 * pu16
pu22_err = 0.20 * pu22

pu16 = pu16 * 1d-3 & pu22 = pu22 * 1d-3				; Convert from mJy to Jy
pu16_err = pu16_err * 1d-3 & pu22_err = pu22_err * 1d-3

pufilelist = pu_tag

; Kluge - add blank entry for cso006 (no dedicated peakup)

pufilelist = [pufilelist,'cso006']
pu16 = [pu16,0.]
pu22 = [pu22,0.]
pu16_err = [pu16_err,0.]
pu22_err = [pu22_err,0.]

save, file='~/Astronomy/Research/Spitzer/cso/data/idl_sav/pu_sur_cso.sav', $
	pufilelist, pu16, pu22, pu16_err, pu22_err


if keyword_set(stop) then stop

end
