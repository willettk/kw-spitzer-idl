pro hac_norm, stop = stop, arr, tau, intflux
;+
; NAME:
;       
;	HAC_NORM
;
; PURPOSE:
;
;	Measure the aliphatic hydrocarbons in Spitzer IRS data using normalized flux units
;
;	- peak depth in normalized flux (e^tau, essentially)
;
;	- integrated normalized flux [cm^-1]
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
;       Written by K. Willett                Feb 09
;-

device, window_state = state

limfac = 3.			; Factor by which to multiply positive flux (for a conservative upper limit)
c = 299792.458d * 1d9		; um / s

; Read in the full LR spectrum

;arr = 'mega'+string([2,4,5,8,10,13,18,32],format='(i03)')
arr = 'arch'+string([3,4,5,7,9,10,13,14,16,17,18,24,26,29,30,32,33,35,36,48],format='(i03)')
tau = dblarr(2,n_elements(arr))
intflux = dblarr(2,n_elements(arr))

for k = 0, n_elements(arr)-1 do begin
	
	fname = arr[k]
	
	tag,fname,dirtag
	
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
	
	wave = sed.wave_lr
	flux = sed.flux_lr
	
	pivarr = [5.2,5.6,6.4,8.1,14.0,26.0]
	npiv = n_elements(pivarr)
	piv = fltarr(npiv)
	for j=0,npiv - 1 do piv[j] = closetomed(wave,flux,pivarr[j])
	
	; 6.85 um feature

	bind685 = closetomed(wave,flux,6.7)
	eind685 = closetomed(wave,flux,7.0)
	
	; Spline fit to continuum
	
	wp = wave(piv)
	fp = flux(piv)
	pspl1 = spl_init(wp,fp)
	splfit = spl_interp(wp,fp,pspl1,wave)
	
	; Add the flux over designated area
	
	spec_sub = flux - splfit
	addfeat_685 = fltarr(2,eind685 - bind685 + 1)
	for i = bind685, eind685 do begin
		addfeat_685(0,i-bind685) = spec_sub(i)
		addfeat_685(1,i-bind685) = (wave(i+1) - wave(i)) * c / (wave(i))^2
	endfor
	
	hacflux_685arr = addfeat_685(0,*) * addfeat_685(1,*) 					; Flux density in Jy, d_nu in Hz
	hacflux_685 = total(hacflux_685arr) * 1d-30						; Flux in W/cm^2
	cont685 = splfit(closetomed(wave,flux,6.85))
	hacew685 = hacflux_685 / cont685 * (6.85)^2 / c * 1d30				; EW in um
	baseerr685 = stddev(flux(closetomed(wave,flux,6.7):closetomed(wave,flux,7.0)))	; Estimate error in baseline flux density [Jy]
	hacew685err = baseerr685 * hacflux_685 / cont685^2 * (6.85)^2/ c*1d30		; Error in equivalent width
	
	; 7.25 um feature

	bind725 = closetomed(wave,flux,7.15)
	eind725 = closetomed(wave,flux,7.3)
	
	; Add the flux over designated area
	
	addfeat_725 = fltarr(2,eind725 - bind725 + 1)
	for i = bind725, eind725 do begin
		addfeat_725(0,i-bind725) = spec_sub(i)
		addfeat_725(1,i-bind725) = (wave(i+1) - wave(i)) * c / (wave(i))^2
	endfor
	
	hacflux_725arr = addfeat_725(0,*) * addfeat_725(1,*) 					; Flux density in Jy, d_nu in Hz
	hacflux_725 = total(hacflux_725arr) * 1d-30						; Flux in W/cm^2
	cont725 = splfit(closetomed(wave,flux,7.25))
	hacew725 = hacflux_725 / cont725 * (7.25)^2 / c * 1d30				; EW in um
	baseerr725 = stddev(flux[closetomed(wave,flux,7.1):closetomed(wave,flux,7.4)])	; Estimate error in baseline flux density [Jy]
	hacew725err = baseerr725 * hacflux_725 / cont725^2 * (7.25)^2/ c*1d30		; Error in equivalent width
	
	; Optical depth
	
	optarr = alog(flux / splfit)
	
	; Stats

	minopt_685 = min(optarr[bind685:eind685])
	minopt_725 = min(optarr[bind725:eind725])
	print,''
	print,'6.85 optical depth for ',fname,' : ',string(minopt_685,format='(f6.3)')
	print,'6.85 int. flux for ',fname,'     : ',string(hacflux_685,format='(e9.2)')
	print,'7.25 optical depth for ',fname,' : ',string(minopt_725,format='(f6.3)')
	print,'7.25 int. flux for ',fname,'     : ',string(hacflux_725,format='(e9.2)')

	tau[0,k] = minopt_685
	tau[1,k] = minopt_725

	intflux[0,k] = hacflux_685
	intflux[1,k] = hacflux_725

	; Plot the results
	
	if not keyword_set(noplot) then begin
		
		;junk = findfile('~/Astronomy/Research/Spitzer/'+dirtag+'/plots/hac/')
		;if junk eq '' then spawn, 'mkdir '+'~/Astronomy/Research/Spitzer/'+dirtag+'/plots/hac/'
		plotname = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/hac/'+fname+'_splfitnu.ps'
		!p.multi = [0,1,2]
		if keyword_set(ps) then begin
			set_plot,'ps'
			device,filename = plotname, /color
			cs = 1
			ls = 2
		endif else begin
			cs = 2
			ls = 1
		endelse
		
		xr685 = [5,10]
		yr685 = [1d-4,1d0]
		
		; 6.2 um PAH
		
		plot, wave, flux, $
			xtitle = 'Wavelength [!7l!3m]', $
			ytitle = 'Flux [Jy]', $
		;	/xlog, $
			/ylog, $
			xrange = xr685, $
			yrange = yr685, $
			/xstyle, $
			/ystyle, $
			charsize = cs, $
			thick = ls, $
			charthick = ls, $
			psym = 10, $
			title = sed.obj
		
		oplot, wp, fp, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
		oplot, wave, splfit, color=fsc_color("Green"), thick = ls
		ver, wave(bind685), linestyle = 1
		ver, wave(eind685), linestyle = 1
		xyouts, 0.8,0.50, fname, /normal
		
		; Optical depth
		
		xr11 = xr685
		yr11 = [0,2]
		
		plot, wave, optarr, $
			psym = 10, $
			xtitle = 'Wavelength [!7l!3m]', $
			ytitle = 'Optical depth !7s!3', $
		;	/xlog, $
	;		/ylog, $
			xrange = xr11, $
			/xstyle, $
	;		yrange = [0,2], /ystyle = 1, $
			charsize = cs, $
			thick = ls, $
			charthick = ls, $
			title = sed.obj
		
		ver, wave(bind685), linestyle = 1
		ver, wave(eind685), linestyle = 1
	
		hor, 0, linestyle = 2
		
		if keyword_set(ps) then begin
			device,/close
			set_plot,'x'
		endif
		
	endif	; noplot
	
;	if fname eq 'mega008' then stop
;wait,10

endfor

print,''

!p.multi=[0,1,1]

if keyword_set(stop) then stop

end
