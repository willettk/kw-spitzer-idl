pro h2_lorescon, ps = ps
;+
; NAME: 
;       H2_LORESCON
;
; PURPOSE:
;
;	Measure the masses and excitation temperatures from molecular H2 lines for control sample
;
; CATEGORY:
;	ASTRONOMY; DATA REDUCTION
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
;	IDL> h2_lorescon
;         
; MODIFICATION HISTORY:
;
;	Written by KW - Sep 07
;-

; Window availability
if not keyword_set(ps) then device, decomposed = 1, window_state = state

; Load in the line fluxes and errors from Excel -> csv files


h2s7 = [0.57,0,0,0,0.46,0.15,0,0,0.82,0,0,0] * 1d-21
h2s5 = [4.17,2.55,2.31,0,2.01,2.08,1.54,0,4.89,2.06,3.14,1.15] * 1d-21
h2s3 = [0.70,0.66,0.39,0,1.22,0.52,0.59,0,0.66,0.83,1.06,0.42] * 1d-21
h2s1 = [1.69,1.73,0.70,0,4.10,1.59,0,0,2.17,2.09,2.81,1.07] * 1d-21

h2s7_err = [0.12,0,0,0,0.1,0.1,0,0,0.17,0,0,0] * 1d-21
h2s5_err = [0.48,0.43,0.22,0,0.20,0.38,0.60,0,0.48,0.45,0.19,0.80] * 1d-21
h2s3_err = [0.15,0.14,0.07,0,0.06,0.06,0.61,0,0.09,0.25,0.11,0.06] * 1d-21
h2s1_err = [0.27,0.13,0.40,0,1.87,0.25,0,0,0.25,0.21,0.19,0.18] * 1d-21

; Physical constants for H2 transitions taken from 
; http://www.not.iac.es/instruments/notcam/ReferenceInfo/h2_lines.html

ej = [1015, 2504, 4586, 7197]				; Energy of upper state
aj = [4.8d-3, 9.84d-2, 0.588, 2.00]			; Einstein-A coefficient
lambda = [17.035,9.6649,6.9091,5.5115]			; Rest wavelength of transition
gj = [21,33,45,57]					; Statistical weight of upper level

plancks = 6.6261d-27					; Planck's constant
c = 3d10						; Speed of light [cm/s]
del_ej = plancks * c / lambda				; Energy of transition

; Linear fit for systems with a single excitation temperature

expr = 'p[0] + x*p[1]'
start = [10, 3d-3]
xarr = fillarr(1, 0, 8000)

; Read in the systems with H2 detections

allh2 = indgen(n_elements(h2s7))
noh2 = [3,7]
isthereh2 = setdifference(allh2,noh2)
restore,'~/Astronomy/Research/Spitzer/control/conname.sav'
fnames = conname
nf = n_elements(fnames)
fobj=strarr(nf)
for i = 0,nf-1 do begin
	targets,fnames(i),r,o
	fobj(i)=o
endfor
fnames = fnames(isthereh2)
fobj = fobj(isthereh2)

h2s1 = h2s1(isthereh2) & h2s3 = h2s3(isthereh2) & h2s5 = h2s5(isthereh2) & h2s7 = h2s7(isthereh2) 
h2s1_err = h2s1_err(isthereh2) & h2s3_err = h2s3_err(isthereh2) & h2s5_err = h2s5_err(isthereh2) & h2s7_err = h2s7_err(isthereh2) 


; Plot the excitation diagrams in 2x2 windows

!p.multi = [0,2,2]
pspath = '~/Astronomy/Research/Spitzer/control/plots/'
hlines = ['H!I2!N S(1)', 'H!I2!N S(3)', 'H!I2!N S(5)', 'H!I2!N S(7)']
ncycles = (n_elements(h2s1) / 4) + 1
temparr = fltarr(2,n_elements(h2s1))

for j = 0, ncycles-1 do begin
	if keyword_set(ps) then begin
		set_plot, 'ps'
		device, filename = pspath+'h2_lores_plot'+strtrim(j,2)+'_con.ps', /landscape
		cs = 1
		lthick = 2
	endif else begin
		cs = 2
		lthick = 1
	endelse
	if j lt ncycles-1 then endind = 3 else endind = (n_elements(h2s1) mod 4) - 1
	for i = 0, endind do begin
		ind = 4*j + i
		if not keyword_set(ps) then if state(j) ne 1 then window, j else wset, j
		h = [h2s1(ind),h2s3(ind),h2s5(ind),h2s7(ind)]
		herr = [h2s1_err(ind),h2s3_err(ind),h2s5_err(ind),h2s7_err(ind)]
		hyes = where(h ne 0.)
		nj     = h(hyes)     / (aj(hyes) * del_ej(hyes))
		nj_err = herr(hyes) / (aj(hyes) * del_ej(hyes))

		; Multiple-temperature fit if the H2 S(7) line is present

		if where(hyes eq 3) ne -1 then begin
			
			plot, ej(hyes), alog(nj / gj(hyes)), $
				xtitle = 'Upper energy level E!IJ!N/k!Ib!N [K]', $
				ytitle = 'ln(N!IJ!N / g!IJ!N) + constant', $
				psym = 4, $
				charsize = cs, $
				thick = lthick, $
				title = fobj(ind), $
				xr = [0,9000], $
				yr = [-20, 2]

			oploterror, ej(hyes), alog(nj / gj(hyes)), alog(nj / gj(hyes)) - alog((nj - nj_err) / gj(hyes)), psym = 3, /lobar
			oploterror, ej(hyes), alog(nj / gj(hyes)), -alog(nj / gj(hyes)) + alog((nj + nj_err) / gj(hyes)), psym = 3, /hibar

			plotsym2,1,3
			oplot, ej([2]), alog(nj[2] / gj([2])), psym = 8
	
			hotind = [1,3]
			hotfit = mpfitexpr(expr, ej(hotind), alog(nj(hotind) / gj(hotind)), alog(nj_err / gj(hyes)), $
				perror = perror, start, /quiet)

			print,fnames(ind), 1d/perror(1)
	
			hotarr = fillarr(1,ej(1),ej(3))
			warmarr = fillarr(1,ej(0),ej(1))
			oplot, hotarr, hotfit(0) + hotfit(1) * hotarr, linestyle = 1, thick = lthick
		
			fourarr = [0,1]
			sub_yvalues = ej(fourarr) * hotfit(1) + hotfit(0)
			sub_yvalues = sub_yvalues - sub_yvalues(0)
			new_yvalues = alog(nj(fourarr) / gj(fourarr)) - sub_yvalues
		
			oplot, ej(fourarr), new_yvalues, psym = 2
		
			warmfit = mpfitexpr(expr, ej(fourarr), new_yvalues, replicate(stddev(0.1*new_yvalues),n_elements(new_yvalues)), $
				start, /quiet)
			oplot, warmarr, warmfit(0) + warmarr*warmfit(1), linestyle = 2, thick = lthick
			
			xyouts, ej(hyes)+100, alog(nj / gj(hyes)), hlines(hyes), charsize = cs, charthick = lthick
			
			xyouts, 2000, -2, 'T!Iwarm!N = '+string(-1d/warmfit(1), format = '(i4)')+' K', charsize = cs, charthick = lthick
			xyouts, 5000, -8, 'T!Ihot!N  = '+string(-1d/hotfit(1), format = '(i4)')+' K', charsize = cs, charthick = lthick

			print,'Multiple fit for ',fnames(ind)

			temparr(0,ind) = -1d/warmfit(1)
			temparr(1,ind) = -1d/hotfit(1)

		endif else begin

		; If not, single fit to the H2 S(1) through H2 S(3) lines is done

		plot, ej(hyes), alog(nj / gj(hyes)), $
			xtitle = 'Upper energy level E!IJ!N/k!Ib!N [K]', $
			ytitle = 'ln(N!IJ!N / g!IJ!N)', $
			psym = 4, $
			charsize = cs, $
			thick = lthick, $
			title = fobj(ind), $
			xr = [0,9000], $
			yr = [-20,0]

		oploterror, ej(hyes), alog(nj / gj(hyes)), alog(nj / gj(hyes)) - alog((nj - nj_err) / gj(hyes)), psym = 3, /lobar
		oploterror, ej(hyes), alog(nj / gj(hyes)), -alog(nj / gj(hyes)) + alog((nj + nj_err) / gj(hyes)), psym = 3, /hibar

		plotsym2,1,3
		if (where(hyes eq 2) ne -1) then begin
			fiveind = hyes(where(hyes eq 2))
			oplot, ej(fiveind), alog(nj(fiveind) / gj(fiveind)), psym = 8
		endif

		xyouts, ej(hyes)+100, alog(nj / gj(hyes)), hlines(hyes), $
			charsize = cs

		; Linear fit to the points in log space gives the excitation temperature

		if where(hyes eq 0) ne -1 and where(hyes eq 1) ne -1 then begin
			onethree = [0,1]
			result = mpfitexpr(expr, ej(onethree), alog(nj / gj(onethree)), 0.1 * alog(nj / gj(onethree)), start, /quiet)
			oplot, xarr, (result(0) + xarr*result(1)), linestyle = 2, thick = lthick
			xyouts, 4000, -3, 'T!Iex!N = '+string(-1d/result(1), format = '(i4)')+' K', charsize = cs
			temparr(0,ind) = -1d/result(1)

		endif else begin if where(hyes eq 1) ne -1 then begin

			; Plot the median temperature of the measured sample through the H2 S(3) flux

			meanslope = -1d/349
			oplot, xarr, meanslope*xarr + (alog(nj(0)/gj(hyes(0))) - meanslope * ej(hyes(0))), linestyle = 4, thick = lthick
			xyouts, 4000, -3, 'T!Iavg!N = 349 K', charsize = cs
		endif
	endelse
		endelse
	endfor
	if keyword_set(ps) then begin
		device, /close
		set_plot,'x'
	endif
endfor

!p.multi = [0,1,1]

stop
end


