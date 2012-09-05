;pro gas_norm, stop = stop, arr, tau, intflux
;+
; NAME:
;       
;	GAS_NORM
;
; PURPOSE:
;
;	Measure the HR gas-phase absorption in units of normalized flux
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
;	Make table of splines and boundaries so all results can be easily replicated
;
;		Not sure how this is different - discuss diff. between normalized flux and optical depth with Jeremy/Henrik
;
; REVISION HISTORY
;       Written by K. Willett                Mar 09
;-

device, window_state = state

limfac = 3.			; Factor by which to multiply positive flux (for a conservative upper limit)
c = 299792.458d * 1d9		; um / s

; Read in the full LR spectrum

;arr = 'mega'+string([2,4,5,8,10,13,18,32],format='(i03)')
;arr = 'arch'+string([3,4,5,7,9,10,13,14,16,17,18,24,26,29,30,32,33,35,36,48],format='(i03)')
arr = ['arch029']
tau = dblarr(2,n_elements(arr))
intflux = dblarr(2,n_elements(arr))

for k = 0, n_elements(arr)-1 do begin
	
	fname = arr[k]
	
	tag,fname,dirtag
	
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
	
	wave = sed.wave_sh
	flux = sed.flux_sh
	
	pivarr = [13.3,13.5,13.8,14.2]
	npiv = n_elements(pivarr)
	piv = fltarr(npiv)
	for j=0,npiv - 1 do piv[j] = closetomed(wave,flux,pivarr[j])
	
	; 13.7 um feature

	linecen = 13.7

	bind_c2h2 = closetomed(wave,flux,13.67)
	eind_c2h2 = closetomed(wave,flux,13.74)
	
	; Spline fit to continuum
	
	wp = wave[piv]
	fp = flux[piv]
	pspl1  = spl_init(wp,fp)
	splfit = spl_interp(wp,fp,pspl1,wave)
	
	; Add the flux over designated area
	
	spec_sub = flux - splfit
	addfeat_c2h2 = fltarr(2,eind_c2h2 - bind_c2h2 + 1)
	for i = bind_c2h2, eind_c2h2 do begin
		addfeat_c2h2[0,i-bind_c2h2] = spec_sub[i]
		addfeat_c2h2[1,i-bind_c2h2] = (wave[i+1] - wave[i]) * c / (wave[i])^2
	endfor
	
	hacflux_c2h2arr = addfeat_c2h2(0,*) * addfeat_c2h2(1,*) 				; Flux density in Jy, d_nu in Hz
	hacflux_c2h2 = total(hacflux_c2h2arr) * 1d-30						; Flux in W/cm^2
	cont_c2h2 = splfit[closetomed(wave,flux,linecen)]
	hacew_c2h2 = hacflux_c2h2 / cont_c2h2 * (linecen)^2 / c * 1d30				; EW in um
	baseerr_c2h2 = stddev(flux(closetomed(wave,flux,13.5):closetomed(wave,flux,13.6)))	; Estimate error in baseline flux density [Jy]
	hacew_c2h2err = baseerr_c2h2 * hacflux_c2h2 / cont_c2h2^2 * (linecen)^2/ c*1d30		; Error in equivalent width
	
	; Normalized flux
	
	normarr = flux / splfit
	
	; Stats

	minnorm_c2h2 = min(normarr[bind_c2h2:eind_c2h2])
	print,''
	print,'C2H2 13.7 um peak depth for ',fname,' : ',string(minnorm_c2h2,format='(f6.3)')
	print,'C2H2 13.7 um int. flux  for ',fname,' : ',string(hacflux_c2h2,format='(e9.2)')

	tau[0,k] = minnorm_c2h2

	intflux[0,k] = hacflux_c2h2

	; Plot the results
	
	if not keyword_set(noplot) then begin
		
		plotname = '~/Astronomy/Research/Spitzer/'+dirtag+'/plots/gas/'+fname+'_splfitnu.ps'
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
		
		xrc2h2 = [13.2,14.3]
		yrc2h2 = [0.5,1.5]
		
		; HR spectrum
		
		plot, wave, flux, $
			xtitle = 'Wavelength [!7l!3m]', $
			ytitle = 'Flux [Jy]', $
		;	/xlog, $
;			/ylog, $
			xrange = xrc2h2, $
;			yrange = yrc2h2, /ystyle, $
			/xstyle, $
			charsize = cs, $
			thick = ls, $
			charthick = ls, $
			psym = 10, $
			title = sed.obj
		
		oplot, wp, fp, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
		oplot, wave, splfit, color=fsc_color("Green"), thick = ls
		ver, wave(bind_c2h2), linestyle = 1
		ver, wave(eind_c2h2), linestyle = 1
		xyouts, 0.8,0.50, fname, /normal
		
		; Normalized flux
		
		xr_tau = xrc2h2
		yr_tau = [0,1.5]
		
		plot, wave, normarr, $
			psym = 10, $
			xtitle = 'Wavelength [!7l!3m]', $
			ytitle = 'Normalized flux', $
		;	/xlog, $
	;		/ylog, $
			xrange = xr_tau, $
			/xstyle, $
			yrange = yr_tau, $
			/ystyle, $
			charsize = cs, $
			thick = ls, $
			charthick = ls, $
			title = sed.obj
		
		ver, wave(bind_c2h2), linestyle = 1
		ver, wave(eind_c2h2), linestyle = 1
	
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
