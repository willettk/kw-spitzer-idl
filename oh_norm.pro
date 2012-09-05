;pro oh_norm, stop = stop, arr, tau, intflux
;+
; NAME:
;       
;	OH_NORM
;
; PURPOSE:
;
;	Measure the OH 34.6 um absorption in units of normalized flux
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
;
; REVISION HISTORY
;       Written by K. Willett                Sep 09
;-

device, window_state = state


c = 299792.458d * 1d9		; um / s
c_cms = 299792458. * 100
h = 6.626d-27			; erg sec
lsun = 3.862d33			; erg/s
mpc2cm = 3.086d24

; Read in the spectrum

tau = dblarr(2)
intflux = dblarr(2)

;for k = 0, n_elements(arr)-1 do begin
	
	; Specific data for 04454

	dl = 235. * mpc2cm		; cm
	log_L_ohm = 2.95

	restore,'~/Astronomy/Research/Spitzer/henrik/04454-4838.irshigh.xdr'
	
	wave = sed.wave
	flux = sed.flux.jy
	
	restore, '~/Astronomy/Research/Spitzer/archived/data/structures/arch029.sav'
	wave = sed.wave_lh
	flux = sed.flux_lh
	dl = sed.dl * mpc2cm
	log_L_ohm = sed.logoh
	pivarr = [34.0,34.3,34.9,35.3]

	;pivarr = [34.0,34.3,34.8,35.3]
	npiv = n_elements(pivarr)
	piv = fltarr(npiv)
	for j=0,npiv - 1 do piv[j] = closetomed(wave,flux,pivarr[j])
	
	linecen = 34.6

	bind_oh = closetomed(wave,flux,34.5)
	eind_oh = closetomed(wave,flux,34.7)
	
	; Spline fit to continuum
	
	wp = wave[piv]
	fp = flux[piv]
	pspl1  = spl_init(wp,fp)
	splfit = spl_interp(wp,fp,pspl1,wave)
	
	; Add the flux over designated area
	
	spec_sub = flux - splfit
	addfeat_oh = fltarr(2,eind_oh - bind_oh + 1)
	for i = bind_oh, eind_oh do begin
		addfeat_oh[0,i-bind_oh] = spec_sub[i]
		addfeat_oh[1,i-bind_oh] = (wave[i+1] - wave[i]) * c / (wave[i])^2
	endfor
	
	flux_oharr = addfeat_oh(0,*) * addfeat_oh(1,*) 				; Flux density in Jy, d_nu in Hz
	flux_oh = total(flux_oharr) * 1d-30						; Flux in W/cm^2
	cont_oh = splfit[closetomed(wave,flux,linecen)]
	ew_oh = flux_oh / cont_oh * (linecen)^2 / c * 1d30				; EW in um
	baseerr_oh = stddev(flux(closetomed(wave,flux,13.5):closetomed(wave,flux,13.6)))	; Estimate error in baseline flux density [Jy]
	ew_oherr = baseerr_oh * flux_oh / cont_oh^2 * (linecen)^2/ c*1d30		; Error in equivalent width

	gl = 4
	gu = 6
	a_ul = 1.74d-2

	noh = -1d * (ew_oh * 1d-4) * 8d * !dpi * c_cms * gl / (a_ul * (linecen * 1d-4) ^4 * gu)

	gamma_abs = -1d * flux_oh * 1d7 * 4d * !dpi * dl^2 / (h * c_cms / (linecen * 1d-4))

	lambda_ohm = 17.98
	L_ohm = 10.^log_L_ohm * lsun
	gamma_ohm = L_ohm / (h * c_cms / (lambda_ohm))
	
	phi_pump = gamma_ohm / gamma_abs * 100.

	; Normalized flux
	
	normarr = flux / splfit
	
	; Stats

	minnorm_oh = min(normarr[bind_oh:eind_oh])
	print,''
	print,'OH 34.6 um peak depth for IRAS 04454-4838: ',string(minnorm_oh,format='(f6.3)')
	print,'OH 34.6 um peah tau   for IRAS 04454-4838: ',string(-1 * alog(minnorm_oh),format='(f6.3)')
	print,'OH 34.6 um int. flux  for IRAS 04454-4838: ',string(flux_oh,format='(e9.2)')
	print,'OH EW [10^-3 um]: ',string((-1d * ew_oh / 1d-3),format='(f6.2)')
	print,'N_OH [cm^-2]: ',string(noh,format='(e9.1)')
	print,'gamma_abs [ph/s]: ',string(gamma_abs,format='(e9.1)')
	print,'gamma_OHM [ph/s]: ',string(gamma_ohm,format='(e9.1)')
	print,'Pumping efficiency [%]: ',string(phi_pump,format='(e9.1)')

	tau[0] = minnorm_oh

	intflux[0] = flux_oh

	; Plot the results
	
	if not keyword_set(noplot) then begin
		
		plotname = '~/Astronomy/Research/Spitzer/archived/plots/gas/arch005_ohabs.ps'
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
		
		xroh = [33,36]
		yroh = [2,4]
		
		; HR spectrum
		
		plot, wave, flux, $
			xtitle = 'Wavelength [!7l!3m]', $
			ytitle = 'Flux [Jy]', $
		;	/xlog, $
;			/ylog, $
			xrange = xroh, /xstyle, $
;			yrange = yroh, /ystyle, $
			charsize = cs, $
			thick = ls, $
			charthick = ls, $
			psym = 10, $
			title = 'IRAS 04454-4838'
		
		oplot, wp, fp, psym = 4, symsize = 2, color=fsc_color("Orange"), thick = ls
		oplot, wave, splfit, color=fsc_color("Green"), thick = ls
		ver, wave(bind_oh), linestyle = 1
		ver, wave(eind_oh), linestyle = 1
		xyouts, 0.8,0.50, 'IRAS 04454-4838', /normal
		
		; Normalized flux
		
		xr_tau = xroh
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
			title = 'IRAS 04454-4838'
		
		ver, wave(bind_oh), linestyle = 1
		ver, wave(eind_oh), linestyle = 1
	
		hor, 0, linestyle = 2
		
		if keyword_set(ps) then begin
			device,/close
			set_plot,'x'
		endif
		
	endif	; noplot
	
;	if fname eq 'mega008' then stop
;wait,10

;endfor

print,''

!p.multi=[0,1,1]

if keyword_set(stop) then stop

end
