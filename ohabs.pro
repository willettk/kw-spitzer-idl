;pro ohabs, stop = stop, arr, tau, intflux
;+
; NAME:
;       
;	OHABS
;
; PURPOSE:
;
;	Measure the OH- absorption at 34.6 um
;
; INPUTS:
;
;	
;
; OUTPUTS:
;
;	
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	For a few individual objects, spline and boundaries must be manually adjusted to make an accurate flux measurement.
;
; REVISION HISTORY
;
;       Written by K. Willett                Dec 08
;-

device, window_state = state

limfac = 3.			; Factor by which to multiply positive flux (for a conservative upper limit)
c = 299792.458d * 1d9		; um / s

; Read in the full LR spectrum

arr = 'arch'+string([3,8,24,29],format='(i03)')
arr = 'arch'+string([29],format='(i03)')
arr = 'control'+string([25],format='(i03)')
tau = dblarr(2,n_elements(arr))
intflux = dblarr(2,n_elements(arr))

for k = 0, n_elements(arr)-1 do begin
	
	fname = arr[k]
	
	tag,fname,dirtag
	
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
	
	wave = sed.wave_lh
	flux = sed.flux_lh
	err = sed.err_lh
	
	pivarr = [33.8,34.1,35.0]
	npiv = n_elements(pivarr)
	piv = intarr(npiv)
	for j=0,npiv - 1 do piv[j] = closeto(wave,pivarr[j])
	
	; OH

	bindoh = closeto(wave,34.5373)
	eindoh = closetomed(wave,flux,34.688)
	
	; Spline fit to continuum
	
	wp = wave(piv)
	fp = flux(piv)
	pspl1 = spl_init(wp,fp)
	splfit = spl_interp(wp,fp,pspl1,wave)
	
	; Add the flux over designated area
	
	spec_sub = flux - splfit
	addfeat_oh = fltarr(2,eindoh - bindoh + 1)
	for i = bindoh, eindoh do begin
		addfeat_oh(0,i-bindoh) = spec_sub(i)
		addfeat_oh(1,i-bindoh) = (wave(i+1) - wave(i)) * c / (wave(i))^2
	endfor
	
	hacflux_oharr = addfeat_oh(0,*) * addfeat_oh(1,*) 					; Flux density in Jy, d_nu in Hz
	hacflux_oh = total(hacflux_oharr) * 1d-30						; Flux in W/cm^2
	contoh = splfit(closetomed(wave,flux,34.6))
	hacewoh = hacflux_oh / contoh * (34.6)^2 / c * 1d30				; EW in um
	baseerroh = stddev(flux(closetomed(wave,flux,33.8):closetomed(wave,flux,34.2)))	; Estimate error in baseline flux density [Jy]
	hacewoherr = baseerroh * hacflux_oh / contoh^2 * (34.6)^2/ c*1d30		; Error in equivalent width
	
	; Optical depth
	
	optarr = alog(flux / splfit)
	
	; Stats

	minopt_oh = min(optarr[bindoh:eindoh])
	print,''
	print,'OH int.  flux for ',fname,'     : ',string(hacflux_oh,format='(e9.2)')
	print,'OH opt. depth for ',fname,'     : ',string(-1d * minopt_oh,format='(f9.2)')

	tau[0,k] = minopt_oh

	intflux[0,k] = hacflux_oh

	; Plot the results
	
	if not keyword_set(noplot) then begin
		
		;junk = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/plots/hac/')
		;if junk eq '' then spawn, 'mkdir '+'~/Astronomy/Research/Spitzer/'+dirtag+'/plots/hac/'
		plotname = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/hac/'+fname+'_splfitnu.ps'
		!p.multi = [0,1,1]
		if keyword_set(ps) then begin
			set_plot,'ps'
			device,filename = plotname, /color
			cs = 1
			ls = 2
		endif else begin
			cs = 2
			ls = 1
		endelse
		
		xroh = [33,36]
		yroh = [0,3]
		
		; 6.2 um PAH
		
		plot, wave, flux, $
			xtitle = 'Wavelength [!7l!3m]', $
			ytitle = 'Flux [Jy]', $
		;	/xlog, $
;			/ylog, $
			xrange = xroh, $
			yrange = yroh, $
			/xstyle, $
			/ystyle, $
			charsize = cs, $
			thick = ls, $
			charthick = ls, $
			psym = 10, $
			title = sed.obj
		
		oplot, wp, fp, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
		oplot, wave, splfit, color=fsc_color("Green"), thick = ls
		
		oplot, wave[bindoh:eindoh], flux[bindoh:eindoh], psym=10, thick = ls, color=fsc_color("Red")

		xyouts, 0.8,0.50, fname, /normal
		
		if keyword_set(ps) then begin
			device,/close
			set_plot,'x'
		endif
		
	endif	; noplot
	
;	if fname eq 'mega008' then stop

stop

wait,10

endfor

print,''

!p.multi=[0,1,1]

if keyword_set(stop) then stop

end
