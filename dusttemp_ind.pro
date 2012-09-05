;+
; NAME: 
;       DUSTTEMP_IND
;
; PURPOSE:
;
;	Measure blackbody dust temperature for an individual object
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
;
;
; REQUIRES:
;
;	MPLANCK.pro
;
; EXAMPLE:
;
;	IDL> dusttemp_ind('cso001')
;
; NOTES:
;
;         
; MODIFICATION HISTORY:
;
;	Adapted from DUSTTEMP.pro - KW, Oct 09
;-

; Modified Planck function for determining expected flux; adds a dust emissivity below a critical frequency.
; Equation is from Yun & Carilli (2002); the angular size \theta is a free parameter in this model. 

function mplanck, x, p

prefix = 2.8d-8
betab = 1.35

return, prefix * (x^3 * p[1]^2 / (exp(0.048d * x / p[0]) - 1d)) * ( 1d - exp( -((x/2000d)^betab) ) )

end

; Modified Planck function in which angular size is taken as an input

function mplanck_thetafixed, x, p

prefix = 2.8d-8
theta = 0.1		; cso001 (1146+59) - rough estimate from Perlman et al. (2001)
betab = 1.35

return, prefix * (x^3 * theta^2 / (exp(0.048d * x / p[0]) - 1d)) * ( 1d - exp( -((x/2000d)^betab) ) )

end

; Begin program

pro dusttemp_ind, obj, notheta=notheta, silent = silent

	; Restore spectral information from .sav file

	tag, obj, dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+obj+'.sav'
	
	; Wavelengths of photometric data points
	
	pu1 = 16.
	pu2 = 22.
	
	iras1 = 12.
	iras2 = 25.
	iras3 = 60.
	iras4 = 100.
	
	iso1 = 170.
	iso2 = 150.
	iso3 = 200.
	iso4 = 90.

	c = 3d10
;	waves = [pu1,pu2,iras1,iras2,iras3,iras4,iso1,iso2,iso3,iso4] 		
	waves = [pu1,pu2,iras1,iras2,iras3,iras4] 				; Wavelengths in um
	allfreq = double((c / (1d-4 * waves)) / 1d9)				; Frequencies in GHz

	; Read in photometry from structure
	
;	if keyword_set(cso) or keyword_set(arch) or keyword_set(control) then begin
		sd = [sed.peakup,sed.iras]
		sderr = [sed.peakuperr,sed.iraserr]
		waves = waves(0:5) & allfreq = allfreq(0:5)
		sd_names = ['Blue PU','Red PU','IRAS 12','IRAS 25','IRAS 60','IRAS 100']+replicate(' ',n_elements(sd))
;	endif else begin
;		sd = [sed.peakup,sed.iras,sed.iso]
;		sderr = [sed.peakuperr,sed.iraserr,sed.isoerr]
;		sd_names = ['Blue PU','Red PU','IRAS 12','IRAS 25','IRAS 60','IRAS 100','ISO 170','ISO 150','ISO 200','ISO 90'] $
;			+replicate(' ',n_elements(sd))
;	endelse
	
	isthereflux = where(sd ne 0)
	sd = sd(isthereflux) 
	sderr = sderr(isthereflux)
	freqs = allfreq(isthereflux)
	
	if keyword_set(verbose) then print, 'Fitting '+sed.tag+' to: '+strjoin(sd_names(isthereflux))
	
	; Fit the photometry to the modified blackbody

	if keyword_set(notheta) then begin
		start = [100]
		result = mpfitfun('mplanck_thetafixed',freqs,sd,sderr,start,/quiet,perror=perror)
	endif else begin
		start = [100,1]
		result = mpfitfun('mplanck',freqs,sd,sderr,start,/quiet,perror=perror)
	endelse
	
	xarr = fillarr(1,1d1,1d5)
;	!p.multi = [0,1,1]

	; Plot the blackbody curve along with the data points
	
	if keyword_set(notheta) then begin
		plot,xarr,mplanck_thetafixed(xarr,[50,1]), $
			/xlog,/ylog, $
			yrange = [1d-6,1d3], $
			xrange = [1d1,1d5], $
			xtitle = 'Frequency [GHz]', $
			ytitle = 'Flux density [Jy]', $
			title = sed.obj, $
			/nodata
		oplot,xarr,mplanck_thetafixed(xarr,result),color=fsc_color("Yellow")
	endif else begin
		plot,xarr,mplanck(xarr,[50,1]), $
			/xlog,/ylog, $
			yrange = [1d-6,1d3], $
			xrange = [1d1,1d5], $
			xtitle = 'Frequency [GHz]', $
			ytitle = 'Flux density [Jy]', $
			title = sed.obj, $
			/nodata
		oplot,xarr,mplanck(xarr,result),color=fsc_color("Yellow")
		psize = result(1) / 206265d * sed.dl * 1d6
		psizeerr = perror(1) / 206265d * sed.dl * 1d6
		if not keyword_set(silent) then begin
			xyouts,0.1,0.7,/normal,'!7h!3 = ('+string(result(1),format='(f6.3)')+' +- '$
				+string(perror(1),format='(e7.1)')+') arcsec', charsize=2
			xyouts,0.1,0.6,/normal,'Gal. diameter = '+string(psize,format='(i5)')+' pc', charsize=2
		endif
	endelse
	
		if not keyword_set(silent) then begin
			xyouts,0.1,0.8,/normal,'T!Id!N = ('+string(result(0),format='(f5.1)')+' +- '$
				+string(perror(0),format='(f4.1)')+') K', charsize=2
		endif

	oploterror, freqs, sd, sderr, psym = symcat(14), color=fsc_color("Red"), symsize = 2
	ver, (c / (1d-4 * 16)) * 1d-9,color=fsc_color("Blue")	
	ver, (c / (1d-4 * 22)) * 1d-9,color=fsc_color("Red")	

	ver, 2d3, linestyle=1
	
end
