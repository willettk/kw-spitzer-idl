pro opr_hight, ps=ps
;+
; NAME:
;       
;	OPR_HIGHT
;
; PURPOSE:
;
;	Compute ortho-to-para ratio of H2 in high-temperature limit for objects in which the S(0)-S(3) lines are all detected
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
;	Based off Roussel et al. (2007)
;
; REVISION HISTORY
;       Written by K. Willett                Feb 10
;-

; Test data - IRAS 04454-4838 (OHM)

fluxarr = [[0.76,3.05,1.27,1.04], $
	[1.27,2.91,1.46,2.21], $
	[2.85,3.73,1.85,1.93], $
	[2.21,1.28,0.36,0.47], $
	[1.30,3.71,1.66,1.87], $
	[0.94,0.92,0.35,0.78], $
	[0.59,2.21,0.79,0.77], $
	[1.02,0.94,0.46,0.65], $
	[1.79,5.24,2.44,3.15]]

ej = [510, 1015, 1682, 2504]			; Energy of upper state
aj = [3d-4, 4.8d-3, 2.76d-2, 9.84d-2]		; Einstein-A coefficient
lambda = [28.221,17.035,12.279,9.6649]		; Rest wavelength of transition
gj = [5,21,9,33]				; Statistical weight of upper level

red  = fsc_color("Red")
blue  = fsc_color("Blue")
green  = fsc_color("Green")

;fluxarr = [[13.2,32.1,10.2,18.7], $
;	[58.2,123.6,60.4,46.8]]
;
!p.multi=[0,3,3]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/plots/opr_hight.ps', /color,/landscape
	defcolor=fsc_color("Black")
	thick=3
	lthick=2 
	cthick = 3
endif else begin
	defcolor=fsc_color("White")
	cthick = 1
endelse

xarr = fillarr(0.01,0,3)

arrsize = size(fluxarr,/dim)
for i=0, arrsize[1] - 1 do begin
	
	flux = fluxarr[*,i]

	t_s0_s2 = (ej[2] - ej[0]) / alog(flux[0] / flux[2] * aj[2]/aj[0] * lambda[0]/lambda[2] * gj[2]/gj[0])
	t_s1_s3 = (ej[3] - ej[1]) / alog(flux[1] / flux[3] * aj[3]/aj[1] * lambda[1]/lambda[3] * gj[3]/gj[1])

	; J is in the upper state, so J for S(0) is 2

	r_s2_s3 = flux[2] / flux[3] * aj[3]/aj[2] * lambda[2]/lambda[3] * (2*(3+2) + 1) / (2*(2+2) + 1) 
	r_s1_s2 = flux[2] / flux[1] * aj[1]/aj[2] * lambda[2]/lambda[1] * (2*(1+2) + 1) / (2*(2+2) + 1)
	r_s0_s1 = flux[0] / flux[1] * aj[1]/aj[0] * lambda[0]/lambda[1] * (2*(1+2) + 1) / (2*(0+2) + 1)

	plot, indgen(10), $
		/nodata, $
		xr=[0,5], $
		yr=[0,800], $
		/xstyle, /ystyle, $
		thick = thick, $
		xthick = thick, $
		ythick = thick, $
		xtitle='OPR!Ihigh T!N', $
		ytitle='T [K]'
	
	;oplot, (euo - eup) / alog(opr * r)
	yred   =  (ej[1] - ej[0]) / alog(r_s0_s1 * xarr)
	ygreen =  (ej[1] - ej[2]) / alog(r_s1_s2 * xarr)
	yblue  =  (ej[3] - ej[2]) / alog(r_s2_s3 * xarr)
	
	redgreen = min(where(ygreen gt yred and xarr gt 0.2))
	greenblue = max(where(yblue gt ygreen and xarr gt 0.5))

	polyfill, [xarr[redgreen], xarr[greenblue], xarr[greenblue], xarr[redgreen]], $
		[ygreen[redgreen], ygreen[redgreen], ygreen[greenblue] > t_s1_s3, ygreen[greenblue] > t_s1_s3], $
		color=fsc_color("Grey"),/data
	;print, [xarr[redgreen], xarr[greenblue], xarr[greenblue], xarr[redgreen], xarr[redgreen]]
	;print,	[ygreen[redgreen], ygreen[greenblue], ygreen[greenblue], ygreen[redgreen], ygreen[redgreen]]
	

	plots, [0,3], [t_s0_s2, t_s0_s2], linestyle=1
	plots, [0,3], [t_s1_s3, t_s1_s3], linestyle=1
	
	xyouts, 3.6, t_s0_s2, 'S(0)-S(2)'
	xyouts, 3.6, t_s1_s3, 'S(1)-S(3)'
	
	xyouts, 2.0, 700, 'T(S1-S3): '+string(t_s1_s3,format='(i4)')+' K', charsize=1, charthick = cthick
	xyouts, 2.0, 630, 'T(S0-S2): '+string(t_s0_s2,format='(i4)')+' K', charsize=1, charthick = cthick
	
	oplot, xarr, yred,   thick = 2, color=red
	oplot, xarr, ygreen, thick = 2, color=green
	oplot, xarr, yblue,  thick = 2, color=blue
	
	xyouts, 3.2, (ej[1] - ej[0]) / alog(r_s0_s1 * 3), 'S(1)-S(0)', color=red, charthick = cthick
	xyouts, 3.2, (ej[1] - ej[2]) / alog(r_s1_s2 * 3), 'S(2)-S(1)', color=green, charthick = cthick
	xyouts, 3.2, (ej[3] - ej[2]) / alog(r_s2_s3 * 3), 'S(3)-S(2)', color=blue, charthick = cthick

	ver, (xarr[redgreen]+xarr[greenblue])/2., linestyle=2
	xyouts,(xarr[redgreen]+xarr[greenblue])/2., 550, charsize=1, string(xarr[redgreen],format='(f3.1)')+'  '+string(xarr[greenblue],format='(f3.1)'), charthick = cthick
	xyouts,(xarr[redgreen]+xarr[greenblue])/2.+0.3, 480, charsize=1, '('+string((xarr[redgreen]+xarr[greenblue])/2.,format='(f3.1)')+')', charthick = cthick

endfor

if keyword_set(ps) then begin
	device, /close
	set_plot,'x'
endif


end
