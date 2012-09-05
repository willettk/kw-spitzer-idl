;
; NAME:
;   	DEJAILBAR
; PURPOSE:
;    	Removes column-wise additive offsets in Spitzer peakup data by
;	taking arithmetic means of input channels.
;
; CALLING SEQUENCE:
;    	jailbar,inimage,outimage [, /fits]
;
; INPUTS:
;    	inimage = IDL 2D array or string filename of FITS to be reduced.
;
; OUTPUTS:
;
;	outimage = IDL 2D array or string filename of FITS.
;
; KEYWORDS:
;
;	/fits = use to open a FITS file, rather than a 2D array
;
; EXAMPLES:
;
;    	IDL> jailbar,'spitzer.fits', 'spitzer_jb.fits'
;	or
;    	IDL> jailbar,data_array, new_array
;
; MODIFICATION HISTORY
;    	Originally developed by V. Charmandaris, Univ. of Crete/Cornell
; 	FITS option, parametrized, arith. mean - K. Willett, Univ. of Colorado, 2-Feb-07
;
pro dejailbar, inimage, outimage, fits = fits

; If reading a FITS file, open image, find dimensions, and load data array

if keyword_set(fits) then begin
	data = readfits(inimage,header)
	naxis1 = sxpar(header,'NAXIS1')
	naxis2 = sxpar(header,'NAXIS2')
	inimage = data
	cols = naxis1
	rows = naxis2
endif else begin
	insize = size(inimage)
	cols = insize(1)
	rows = insize(2)
endelse

; Create blank array for channel medians; rebin image by the four readout channels

chmed = fltarr(4)
temp=reform(inimage,4, cols / 4 * rows)

; Find the medians of each readout channel

for i = 0,3 do begin
	chmed[i] = median(temp[i,*])
endfor

; Take the arithmetic mean of the channels as the fiducial offest

tieval = mean(chmed)

; Add the difference between the fiducial median and the median for each 
; to normalize the output channels

for j = 0,3 do begin
	temp[j,*] = temp[j,*] + tieval - chmed[j]
endfor

; Rebin back to the original image size

newimage=reform(temp, 128, 128) 

; Write a new FITS file with the jailbars corrected

if keyword_set(fits) then begin
	sxaddhist,'Jailbar corrector applied '+systime(),header
	writefits,outimage,newimage,header
	print,'New FITS file written: '+string(outimage)
endif else begin
	outimage = newimage
endelse

end
