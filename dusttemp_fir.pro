;+
; NAME: 
;       DUSTTEMP_FIR
;
; PURPOSE:
;
;	Measure the far_IR dust temperature from IRS spectra (20-30 um)
;
; INPUTS:
;	
; OUTPUTS:
;
; KEYWORDS:
;
;	CONTROL - 	if set, computes dust temps for the non-masing control sample. Default is to operate on OHMs. 
;
; REQUIRES:
;
;	MPLANCK.pro
;
; EXAMPLE:
;
;	IDL> dusttemp, /control
;
; NOTES:
;
;         
; MODIFICATION HISTORY:
;
;	Adapted from DUSTTEMP.pro - KW, Nov 08
;-

; Modified Planck function for determining expected flux; adds a dust emissivity below a critical frequency.
; Equation is from Yun & Carilli (2002); the angular size \theta is a free parameter in this model. 

function mplanck, x, p

prefix = 2.8d-8
theta = 0.1d
betab = 1.35

return, prefix * (x^3 * p[1]^2 / (exp(0.048d * x / p[0]) - 1d)) * ( 1d - exp( -((x/2000d)^betab) ) )

end

; Begin program

pro dusttemp_fir, control = control, arch  = arch, cso = cso, verbose = verbose, wait = wait, noplot = noplot

for k= 0,1 do begin

	if k eq 0 then fnames = [transpose(ohmdat('tag')),transpose(archdat('tag'))] $
		else fnames = condat('tag')
	
	nf = n_elements(fnames)
	
	; Empty arrays for the dust temperature, angular size, and physical size
	
	dtemp = fltarr(nf,2)
	angsize = fltarr(nf,2)
	physize = fltarr(nf,2)
	
	; Loop over objects
	
	for i = 0, nf - 1 do begin
	
		; Restore spectral information from .sav file
	
		tag,fnames[i],dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fnames(i)+'.sav'
		
		; Use spectra along sampled points instead of photometry
		
		c = 3d10
		waves = fillarr(1,20,30) 					; Wavelengths in um
		freqs = double((c / (1d-4 * waves)) / 1d9)			; Frequencies in GHz
	
		; Read in photometry from structure
		
		sd = fltarr(n_elements(waves))
		sderr = fltarr(n_elements(waves))
		for j = 0, n_elements(waves) - 1 do begin
			tempphot = cmw(sed.wave_lr,sed.flux_lr,sed.err_lr,waves[j])
			sd[j] = tempphot
			bind = closeto(sed.wave_lr,waves[j]-1)
			eind = closeto(sed.wave_lr,waves[j]+1)
			temperr = mean(sed.err_lr[bind:eind])
			sderr[j] = temperr
		endfor
	
		if keyword_set(verbose) then print, 'Fitting '+sed.tag+' to: '+strjoin(sd_names(isthereflux))
		
		; Fit the photometry to the modified blackbody
	
		start = [100,1]
		result = mpfitfun('mplanck',freqs,sd,sderr,start,/quiet,perror=perror)
		
		xarr = fillarr(1,1d1,1d5)
		!p.multi = [0,1,1]
	
		psize = result(1) / 206265d * sed.dl * 1d6
		psizeerr = perror(1) / 206265d * sed.dl * 1d6

		; Plot the blackbody curve along with the data points
		
		if not keyword_set(noplot) then begin

			plot,xarr,mplanck(xarr,[50,1]), $
				/xlog,/ylog, $
				yrange = [1d-6,1d3], $
				xrange = [1d1,1d5], $
				xtitle = 'Frequency [GHz]', $
				ytitle = 'Flux density [Jy]', $
				title = sed.obj, $
				/nodata
			
			oploterror, freqs, sd, sderr, psym = symcat(14), color=fsc_color("Red"), symsize = 2
			ver, (c / (1d-4 * 16)) * 1d-9,color=fsc_color("Blue")	
			ver, (c / (1d-4 * 22)) * 1d-9,color=fsc_color("Red")	
				
			oplot,xarr,mplanck(xarr,result),color=fsc_color("Yellow")
			
			xyouts,0.1,0.8,/normal,'T!Id!N = ('+string(result(0),format='(f5.1)')+' +- '$
				+string(perror(0),format='(f4.1)')+') K', charsize=2
			xyouts,0.1,0.7,/normal,'!7h!3 = ('+string(result(1),format='(f6.3)')+' +- '$
				+string(perror(1),format='(e7.1)')+') arcsec', charsize=2
			xyouts,0.1,0.6,/normal,'Gal. diameter = '+string(psize,format='(i5)')+' pc', charsize=2

		endif
			
		; Place results in array
	
		dtemp(i,*) = [result(0),perror(0)]
		angsize(i,*) = [result(1),perror(1)]
		physize(i,*) = [psize,psizeerr]
	
		if keyword_set(wait) then wait, 1
	endfor
	
	if k eq 0 then begin
;		print,'OHMs: ',mean(dtemp[*,0]),' +- ',stddev(dtemp[*,0])
		otemp = dtemp[*,0]
	endif else begin
;		print,'Control: ',mean(dtemp[*,0]),' +- ',stddev(dtemp[*,0]) 
		ctemp = dtemp[*,0]
	endelse
endfor

kstwo, otemp, ctemp, D_temp, prob_temp
gauss_temp = sqrt(2d) * inverf(1d - prob_temp)

print,''
print,'D_KS    for temp_fir: '+string(D_temp,format='(f7.3)')
print,'KS-prob for temp_fir: '+string(prob_temp,format='(e9.2)')
print,'Gaussian probability:  '+string(gauss_temp,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(otemp),format='(f7.2)')+' +- '+string(stddev(otemp),format='(f7.2)')
print,'Average value (con) :  '+string(mean(ctemp),format='(f7.2)')+' +- '+string(stddev(ctemp),format='(f7.2)')

end
