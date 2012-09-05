;pro signaltonoise, noplot = noplot
;+
; NAME:
;       
;	SIGNALTONOISE
;
; PURPOSE:
;
;	Measure the signal-to-noise figure of merit for Spitzer IRS spectra
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
;	IDL> .r signaltonoise
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Nov 08
;-

; Load data

for j = 0,3 do begin

	case j of
		0: fname = ohmdat('tag')
		1: fname = archdat('tag')
		2: fname = condat('tag')
		3: fname = csodat('tag')
	endcase

	nf = n_elements(fname)
	
	sn = fltarr(nf)
	chi2 = fltarr(nf)
	
	for i = 0, nf - 1 do begin
	
		tag,fname[i],dirtag
		restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname[i]+'.sav'
		
		wave = sed.wave_lr
		flux = sed.flux_lr
		err = sed.err_lr

		bwave = 19.2
		ewave = 23.0

		bind = closeto(wave, bwave)
		eind = closeto(wave, ewave)
		
		; Fit baseline

;		expr = 'p[0] + p[1] * x'
;		start = [0., 0.]
		expr = 'p[0] + p[1] * x + p[2] * x^2'
		start = [0., 0., 0.]

		result = mpfitexpr(expr,wave[bind:eind],flux[bind:eind],err[bind:eind], start, bestnorm = bestnorm, dof = dof, /quiet)

;		fluxfit = result[0] + result[1] * wave
		fluxfit = result[0] + result[1] * wave + result[2] * wave * wave
		bsubtracted = flux - fluxfit

		signal = cmw(wave, flux, err, mean([bwave,ewave]))
		noise = stddev(bsubtracted[bind:eind])
		
		if not keyword_set(noplot) then begin

			!p.multi=[0,2,1]
	
			plot,wave, flux, $
				/xlog, /ylog, $
				xrange = [bwave-1,ewave+1], /xstyle, $
				xtitle = 'Wavelength [um]', $
				ytitle = 'Flux [Jy]', $
				title = sed.obj, $
				psym = 10

			oplot, wave[bind:eind], flux[bind:eind], thick = 2, color=fsc_color("Red"), psym = 10
			oplot, wave, fluxfit, linestyle = 2
				
			plot, wave, bsubtracted, $
				/xlog, $
;				/ylog, $
				xrange = [bwave-1,ewave+1], /xstyle, $
				xtitle = 'Wavelength [um]', $
				ytitle = 'Flux [Jy]', $
				title = sed.obj, $
				psym = 10

			oplot, wave[bind:eind], bsubtracted[bind:eind], thick = 2, color=fsc_color("Red"), psym = 10

			hor, 0.0, linestyle=1
				
		endif
		
		print,''
		print,'Reduced chi^2 of fit', bestnorm/dof
		print,'S/N = ',signal/noise
		print,''

		sn[i] = signal/noise
		chi2[i] = bestnorm/dof

	endfor

	fname_sn = fname

	case j of
		0: save, sn, chi2, fname_sn, file='~/Astronomy/Research/Spitzer/ohm/data/idl_sav/sn.sav' 
		1: save, sn, chi2, fname_sn, file='~/Astronomy/Research/Spitzer/archived/data/idl_sav/sn.sav' 
		2: save, sn, chi2, fname_sn, file='~/Astronomy/Research/Spitzer/control/data/idl_sav/sn.sav' 
		3: save, sn, chi2, fname_sn, file='~/Astronomy/Research/Spitzer/cso/data/idl_sav/sn.sav' 
	endcase

endfor

restore, '~/Astronomy/Research/Spitzer/ohm/data/idl_sav/sn.sav'
sn_ohm = sn & fname_ohm = transpose(fname_sn)
restore, '~/Astronomy/Research/Spitzer/archived/data/idl_sav/sn.sav'
sn_arch = sn & fname_arch = transpose(fname_sn)

f = [fname_ohm,fname_arch]
s = [sn_ohm,sn_arch]

obj = strarr(n_elements(f))

for i = 0, n_elements(f)-1 do begin
	targets,f[i],r,o
	obj[i] = o
endfor

sind = sort(obj)

f = f[sind]
s = s[sind]

;print,[transpose(string(s,format='(i5)'))]

restore, '~/Astronomy/Research/Spitzer/control/data/idl_sav/sn.sav'
sn_con = sn & fname_con = transpose(fname_sn)

obj = strarr(n_elements(sn_con))

for i = 0, n_elements(sn_con)-1 do begin
	targets,fname_con[i],r,o
	obj[i] = o
endfor

sind = sort(obj)

fname_con = fname_con[sind]
sn_con = sn_con[sind]

;print,[transpose(string(sn_con,format='(i5)'))]

;print,median([s,sn_con])

end
