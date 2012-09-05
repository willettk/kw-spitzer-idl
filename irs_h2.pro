pro irs_h2, ps = ps
;+
; NAME: 
;       IRS_H2 
;
; PURPOSE:
;
;	Measure the masses and excitation temperatures from molecular H2 lines 
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
;	IDL> irs_h2
;         
; MODIFICATION HISTORY:
;
;	Written by KW - Aug 07
;-

; Window availability
if not keyword_set(ps) then device, decomposed = 1, window_state = state

; Load in the line fluxes and errors from Excel -> csv files

fpath = '~/Astronomy/Research/Spitzer/'
fluxfile = 'mega_lines_aug29.csv'
errfile = 'mega_lines_aug29_errs.csv'

readcol, fpath+fluxfile, fname, ftag, h2s7, h2s3, h2s2, h2s1_sh, h2s1_lh, h2s0, $
	format = 'a,a,f,f,f,f,f', skipline = 4, numline = 22;, /silent

readcol, fpath+errfile, fname, ftag, h2s7_err, h2s3_err, h2s2_err, h2s1_sh_err, h2s1_lh_err, h2s0_err, $
	format = 'a,a,f,f,f,f,f', skipline = 4, numline = 22;, /silent

; Physical constants for H2 transitions taken from 
; http://www.not.iac.es/instruments/notcam/ReferenceInfo/h2_lines.html

ej = [510, 1015, 1015, 1682, 2504, 7197]			; Energy of upper state
aj = [3d-4, 4.8d-3,4.8d-3, 2.76d-2, 9.84d-2, 2.00]		; Einstein-A coefficient
lambda = [28.221, 17.035,17.035, 12.279, 9.6649, 5.5115]	; Wavelength of transition
plancks = 6.6261d-27					; Planck's constant
c = 3d10						; Speed of light [cm/s]
del_ej = plancks * c / lambda				; Energy of transition
gj = [5,21,21,9,33,57]					; Statistical weight of upper level

; Linear fit for systems with a single excitation temperature

expr = 'p[0] + x*p[1]'
start = [10, 3d-3]
xarr = fillarr(1, 0, 8000)

; Read in the systems with H2 detections

allh2 = indgen(22)
noh2 = [1,2,6,7]			; No H2 seen in mega006, mega007, mega014, mega016
isthereh2 = setdifference(allh2,noh2)
fname = fname(isthereh2)
ftag = ftag(isthereh2)

h2s3 = h2s3(isthereh2) & h2s2 = h2s2(isthereh2) & h2s1_lh = h2s1_lh(isthereh2) 
h2s1_sh = h2s1_sh(isthereh2) & h2s0 = h2s0(isthereh2) & h2s7 = h2s7(isthereh2) 

h2s3_err = h2s3_err(isthereh2) & h2s2_err = h2s2_err(isthereh2) & h2s1_lh_err = h2s1_lh_err(isthereh2) 
h2s1_sh_err = h2s1_sh_err(isthereh2) & h2s0_err = h2s0_err(isthereh2) & h2s7_err = h2s7_err(isthereh2) 

; Plot the excitation diagrams in 2x2 windows

!p.multi = [0,2,2]
pspath = '~/Astronomy/Research/Spitzer/ohm/lines/plots/'
hlines = ['H!I2!N S(0)', 'H!I2!N S(1) [LH]', 'H!I2!N S(1) [SH]', 'H!I2!N S(2)', 'H!I2!N S(3)', 'H!I2!N S(7)']

for j = 0, 4 do begin
	if keyword_set(ps) then begin
		set_plot, 'ps'
		device, filename = pspath+'h2_plot'+strtrim(j,2)+'.ps', /landscape
		cs = 1
	endif else cs = 2
	if j lt 4 then endind = 3 else endind = 0
	for i = 0, endind do begin
		ind = 4*j + i
		if not keyword_set(ps) then if state(j) ne 1 then window, j else wset, j
		h = [h2s0(ind),h2s1_lh(ind),h2s1_sh(ind),h2s2(ind),h2s3(ind),h2s7(ind)]
		h_err = [h2s0_err(ind),h2s1_lh_err(ind),h2s1_sh_err(ind),h2s2_err(ind),h2s3_err(ind),h2s7_err(ind)]
		hyes = where(h ne 0.)
		nj     = h(hyes)     / (aj(hyes) * del_ej(hyes))
		nj_err = h_err(hyes) / (aj(hyes) * del_ej(hyes))
		plot, ej(hyes), alog(nj / gj(hyes)), $
			xtitle = 'Upper energy level E!IJ!N/k!Ib!N [K]', $
			ytitle = 'ln(N!IJ!N / g!IJ!N)', $
			psym = 4, $
			charsize = cs, $
			title = fname(ind)+' ('+ftag(ind)+')', $
			xr = [0,8000], $
			yr = [-20,0]
		oploterror, ej(hyes), alog(nj / gj(hyes)), alog(nj / gj(hyes)) - alog((nj - nj_err) / gj(hyes)), psym = 3, /lobar
		oploterror, ej(hyes), alog(nj / gj(hyes)), -alog(nj / gj(hyes)) + alog((nj + nj_err) / gj(hyes)), psym = 3, /hibar
		xyouts, ej(hyes)+100, alog(nj / gj(hyes)), hlines(hyes), $
			charsize = cs

		; Linear fit to the points in log space gives the excitation temperature

		if n_elements(hyes) ge 2 then begin
			result = mpfitexpr(expr, ej(hyes), alog(nj / gj(hyes)), alog(nj_err / gj(hyes)), start, /quiet)
			oplot, xarr, (result(0) + xarr*result(1)), linestyle = 2
			xyouts, 4000, -3, 'T!Iex!N = '+string(-1d/result(1), format = '(i4)')+' K', charsize = cs
		endif
	endfor
	if keyword_set(ps) then begin
		device, /close
		set_plot,'x'
	endif
endfor

!p.multi = [0,1,1]

stop
end
