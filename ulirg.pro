pro ulirg, ps = ps
;+
; NAME:
;       ULIRG
;
; PURPOSE:
; 	Read and plot ULIRG data to compare my measurements with published data
;
; INPUTS:
;
; OUTPUTS:
;
; 	Several scatterplots of data for three targets, both nods and coadded data
;
; KEYWORDS:
;
;	PS - hard copies of data
;
; EXAMPLE:
;	IDL> ulirg, /ps
;
; REQUIRES:
;
; NOTES:
;
;	- Data ("true" values) are published in Armus et al (2007), ApJ, 656:148-167
;
; REVISION HISTORY
;       Written by K. Willett                Jun 2007
;	Added analysis of IRSCLEANed data - KW, Jul 2007
;-

device, decomposed = 1

fpath = '~/Astronomy/Research/Spitzer/ulirg/'

readcol, fpath+'ulirg_fits.txt', obj, iras_flux, iras_err, iras_ew, ugc_flux, ugc_err, ugc_ew, $
	mrk_flux, mrk_err, mrk_ew, format = 'a,f,f,f,f,f,f,f,f,f', /silent

readcol, fpath+'ulirg_true.txt', obj, iras_flux_true, iras_err_true, iras_ew_true, ugc_flux_true, $
	ugc_err_true, ugc_ew_true, mrk_flux_true, mrk_err_true, mrk_ew_true, format = 'a,f,f,f,f,f,f,f,f,f', $
	/silent

readcol, fpath+'ulirg_nods.txt', obj, wave, $
	iras_flux_1p, iras_err_1p, iras_ew_1p, iras_flux_2p, iras_err_2p, iras_ew_2p, $
	ugc_flux_1p, ugc_err_1p, ugc_ew_1p, ugc_flux_2p, ugc_err_2p, ugc_ew_2p, $
	mrk_flux_1p, mrk_err_1p, mrk_ew_1p, mrk_flux_2p, mrk_err_2p, mrk_ew_2p, $
	format = 'a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f', /silent

readcol, fpath+'ulirg_clean.txt', $
	iras_clean_flux, iras_clean_err, iras_clean_ew, $
	ugc_clean_flux, ugc_clean_err, ugc_clean_ew, $
	mrk_clean_flux, mrk_clean_err, mrk_clean_ew, /silent

readcol, fpath+'ulirg_clean_nods.txt', $
	iras_clean_flux_1p, iras_clean_err_1p, iras_clean_ew_1p, $
	iras_clean_flux_2p, iras_clean_err_2p, iras_clean_ew_2p, $
	ugc_clean_flux_1p, ugc_clean_err_1p, ugc_clean_ew_1p, $
	ugc_clean_flux_2p, ugc_clean_err_2p, ugc_clean_ew_2p, $
	mrk_clean_flux_1p, mrk_clean_err_1p, mrk_clean_ew_1p, $
	mrk_clean_flux_2p, mrk_clean_err_2p, mrk_clean_ew_2p, /silent

readcol, fpath+'ulirg_lores.txt', $
	mrk_flux_lores, mrk_err_lores, mrk_ew_lores, /silent

set = [0,1,2,4,5]

obj = obj[set]
wave = wave[set]
iras_flux = iras_flux[set] & iras_err = iras_err[set] & iras_ew = iras_ew[set]

iras_flux_true = iras_flux_true[set] & iras_err_true = iras_err_true[set] & iras_ew_true = iras_ew_true[set]
ugc_flux = ugc_flux[set] & ugc_err = ugc_err[set] & ugc_ew = ugc_ew[set]
ugc_flux_true = ugc_flux_true[set] & ugc_err_true = ugc_err_true[set] & ugc_ew_true = ugc_ew_true[set]
mrk_flux = mrk_flux[set] & mrk_err = mrk_err[set] & mrk_ew = mrk_ew[set]
mrk_flux_true = mrk_flux_true[set] & mrk_err_true = mrk_err_true[set] & mrk_ew_true = mrk_ew_true[set]

iras_flux_1p = iras_flux_1p[set] & iras_err_1p = iras_err_1p[set] & iras_ew_1p = iras_ew_1p[set]
ugc_flux_1p = ugc_flux_1p[set] & ugc_err_1p = ugc_err_1p[set] & ugc_ew_1p = ugc_ew_1p[set]
mrk_flux_1p = mrk_flux_1p[set] & mrk_err_1p = mrk_err_1p[set] & mrk_ew_1p = mrk_ew_1p[set]

iras_flux_2p = iras_flux_2p[set] & iras_err_2p = iras_err_2p[set] & iras_ew_2p = iras_ew_2p[set]
ugc_flux_2p = ugc_flux_2p[set] & ugc_err_2p = ugc_err_2p[set] & ugc_ew_2p = ugc_ew_2p[set]
mrk_flux_2p = mrk_flux_2p[set] & mrk_err_2p = mrk_err_2p[set] & mrk_ew_2p = mrk_ew_2p[set]

if keyword_set(ps) then cs = 1.0 else cs = 2.0
if keyword_set(ps) then ss = 1.0 else ss = 2.0
if keyword_set(ps) then ct = 2.0 else ct = 1.0

expr = 'p[0]*x'
expr2 = 'p[1]*x + p[0]'

iras_start = [1]
ugc_start = [1]
mrk_start = [1]

; Uncleaned fits

iras_result = mpfitexpr(expr, iras_flux, iras_flux_true, iras_err_true, iras_start, $
	perror = iras_perror, /quiet, bestnorm = iras_bestnorm, dof = iras_dof)
ugc_result = mpfitexpr(expr, ugc_flux, ugc_flux_true, ugc_err_true, ugc_start, $
	perror = ugc_perror, /quiet, bestnorm = ugc_bestnorm, dof = ugc_dof)
mrk_result = mpfitexpr(expr, mrk_flux, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_perror, /quiet, bestnorm = mrk_bestnorm, dof = mrk_dof)

iras_result_1p = mpfitexpr(expr, iras_flux_1p, iras_flux_true, iras_err_true, iras_start, $
	perror = iras_perror_1p, /quiet, bestnorm = iras_bestnorm_1p, dof = iras_dof_1p)
ugc_result_1p = mpfitexpr(expr, ugc_flux_1p, ugc_flux_true, ugc_err_true, ugc_start, $
	perror = ugc_perror_1p, /quiet, bestnorm = ugc_bestnorm_1p, dof = ugc_dof_1p)
mrk_result_1p = mpfitexpr(expr, mrk_flux_1p, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_perror_1p, /quiet, bestnorm = mrk_bestnorm_1p, dof = mrk_dof_1p)

iras_result_2p = mpfitexpr(expr, iras_flux_2p, iras_flux_true, iras_err_true, iras_start, $
	perror = iras_perror_2p, /quiet, bestnorm = iras_bestnorm_2p, dof = iras_dof_2p)
ugc_result_2p = mpfitexpr(expr, ugc_flux_2p, ugc_flux_true, ugc_err_true, ugc_start, $
	perror = ugc_perror_2p, /quiet, bestnorm = ugc_bestnorm_2p, dof = ugc_dof_2p)
mrk_result_2p = mpfitexpr(expr, mrk_flux_2p, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_perror_2p, /quiet, bestnorm = mrk_bestnorm_2p, dof = mrk_dof_2p)

xarr = fillarr(1, -100, 500)

iras_ew_result = mpfitexpr(expr, iras_ew, iras_ew_true, iras_err_true, iras_start, $
	perror = iras_ew_perror, /quiet, bestnorm = iras_ew_bestnorm, dof = iras_ew_dof)
ugc_ew_result = mpfitexpr(expr, ugc_ew, ugc_ew_true, ugc_err_true, ugc_start, $
	perror = ugc_ew_perror, /quiet, bestnorm = ugc_ew_bestnorm, dof = ugc_ew_dof)
mrk_ew_result = mpfitexpr(expr, mrk_ew, mrk_ew_true, mrk_err_true, mrk_start, $
	perror = mrk_ew_perror, /quiet, bestnorm = mrk_ew_bestnorm, dof = mrk_ew_dof)

iras_ew_result_1p = mpfitexpr(expr, iras_ew_1p, iras_ew_true, iras_err_true, iras_start, $
	perror = iras_ew_perror_1p, /quiet, bestnorm = iras_ew_bestnorm_1p, dof = iras_ew_dof_1p)
ugc_ew_result_1p = mpfitexpr(expr, ugc_ew_1p, ugc_ew_true, ugc_err_true, ugc_start, $
	perror = ugc_ew_perror_1p, /quiet, bestnorm = ugc_ew_bestnorm_1p, dof = ugc_ew_dof_1p)
mrk_ew_result_1p = mpfitexpr(expr, mrk_ew_1p, mrk_ew_true, mrk_err_true, mrk_start, $
	perror = mrk_ew_perror_1p, /quiet, bestnorm = mrk_ew_bestnorm_1p, dof = mrk_ew_dof_1p)

iras_ew_result_2p = mpfitexpr(expr, iras_ew_2p, iras_ew_true, iras_err_true, iras_start, $
	perror = iras_ew_perror_2p, /quiet, bestnorm = iras_ew_bestnorm_2p, dof = iras_ew_dof_2p)
ugc_ew_result_2p = mpfitexpr(expr, ugc_ew_2p, ugc_ew_true, ugc_err_true, ugc_start, $
	perror = ugc_ew_perror_2p, /quiet, bestnorm = ugc_ew_bestnorm_2p, dof = ugc_ew_dof_2p)
mrk_ew_result_2p = mpfitexpr(expr, mrk_ew_2p, mrk_ew_true, mrk_err_true, mrk_start, $
	perror = mrk_ew_perror_2p, /quiet, bestnorm = mrk_ew_bestnorm_2p, dof = mrk_ew_dof_2p)

; Cleaned fits

iras_clean_result = mpfitexpr(expr, iras_clean_flux, iras_flux_true, iras_err_true, iras_start, $
	perror = iras_clean_perror, /quiet, bestnorm = iras_clean_bestnorm, dof = iras_clean_dof)
ugc_clean_result = mpfitexpr(expr, ugc_clean_flux, ugc_flux_true, ugc_err_true, ugc_start, $
	perror = ugc_clean_perror, /quiet, bestnorm = ugc_clean_bestnorm, dof = ugc_clean_dof)
mrk_clean_result = mpfitexpr(expr, mrk_clean_flux, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_clean_perror, /quiet, bestnorm = mrk_clean_bestnorm, dof = mrk_clean_dof)

iras_clean_result_1p = mpfitexpr(expr, iras_clean_flux_1p, iras_flux_true, iras_err_true, iras_start, $
	perror = iras_clean_perror_1p, /quiet, bestnorm = iras_clean_bestnorm_1p, dof = iras_clean_dof_1p)
ugc_clean_result_1p = mpfitexpr(expr, ugc_clean_flux_1p, ugc_flux_true, ugc_err_true, ugc_start, $
	perror = ugc_clean_perror_1p, /quiet, bestnorm = ugc_clean_bestnorm_1p, dof = ugc_clean_dof_1p)
mrk_clean_result_1p = mpfitexpr(expr, mrk_clean_flux_1p, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_clean_perror_1p, /quiet, bestnorm = mrk_clean_bestnorm_1p, dof = mrk_clean_dof_1p)

iras_clean_result_2p = mpfitexpr(expr, iras_clean_flux_2p, iras_flux_true, iras_err_true, iras_start, $
	perror = iras_clean_perror_2p, /quiet, bestnorm = iras_clean_bestnorm_2p, dof = iras_clean_dof_2p)
ugc_clean_result_2p = mpfitexpr(expr, ugc_clean_flux_2p, ugc_flux_true, ugc_err_true, ugc_start, $
	perror = ugc_clean_perror_2p, /quiet, bestnorm = ugc_clean_bestnorm_2p, dof = ugc_clean_dof_2p)
mrk_clean_result_2p = mpfitexpr(expr, mrk_clean_flux_2p, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_clean_perror_2p, /quiet, bestnorm = mrk_clean_bestnorm_2p, dof = mrk_clean_dof_2p)

xarr = fillarr(1, -100, 500)

iras_clean_ew_result = mpfitexpr(expr, iras_clean_ew, iras_ew_true, iras_err_true, iras_start, $
	perror = iras_clean_ew_perror, /quiet, bestnorm = iras_clean_ew_bestnorm, dof = iras_clean_ew_dof)
ugc_clean_ew_result = mpfitexpr(expr, ugc_clean_ew, ugc_ew_true, ugc_err_true, ugc_start, $
	perror = ugc_clean_ew_perror, /quiet, bestnorm = ugc_clean_ew_bestnorm, dof = ugc_clean_ew_dof)
mrk_clean_ew_result = mpfitexpr(expr, mrk_clean_ew, mrk_ew_true, mrk_err_true, mrk_start, $
	perror = mrk_clean_ew_perror, /quiet, bestnorm = mrk_clean_ew_bestnorm, dof = mrk_clean_ew_dof)

iras_clean_ew_result_1p = mpfitexpr(expr, iras_clean_ew_1p, iras_ew_true, iras_err_true, iras_start, $
	perror = iras_clean_ew_perror_1p, /quiet, bestnorm = iras_clean_ew_bestnorm_1p, dof = iras_clean_ew_dof_1p)
ugc_clean_ew_result_1p = mpfitexpr(expr, ugc_clean_ew_1p, ugc_ew_true, ugc_err_true, ugc_start, $
	perror = ugc_clean_ew_perror_1p, /quiet, bestnorm = ugc_clean_ew_bestnorm_1p, dof = ugc_clean_ew_dof_1p)
mrk_clean_ew_result_1p = mpfitexpr(expr, mrk_clean_ew_1p, mrk_ew_true, mrk_err_true, mrk_start, $
	perror = mrk_clean_ew_perror_1p, /quiet, bestnorm = mrk_clean_ew_bestnorm_1p, dof = mrk_clean_ew_dof_1p)

iras_clean_ew_result_2p = mpfitexpr(expr, iras_clean_ew_2p, iras_ew_true, iras_err_true, iras_start, $
	perror = iras_clean_ew_perror_2p, /quiet, bestnorm = iras_clean_ew_bestnorm_2p, dof = iras_clean_ew_dof_2p)
ugc_clean_ew_result_2p = mpfitexpr(expr, ugc_clean_ew_2p, ugc_ew_true, ugc_err_true, ugc_start, $
	perror = ugc_clean_ew_perror_2p, /quiet, bestnorm = ugc_clean_ew_bestnorm_2p, dof = ugc_clean_ew_dof_2p)
mrk_clean_ew_result_2p = mpfitexpr(expr, mrk_clean_ew_2p, mrk_ew_true, mrk_err_true, mrk_start, $
	perror = mrk_clean_ew_perror_2p, /quiet, bestnorm = mrk_clean_ew_bestnorm_2p, dof = mrk_clean_ew_dof_2p)

; Lo-res fits

mrk_lores_result = mpfitexpr(expr, mrk_flux_lores, mrk_flux_true, mrk_err_true, mrk_start, $
	perror = mrk_lores_perror, bestnorm = mrk_lores_bestnorm, dof = mrk_lores_dof)

mrk_lohi_result = mpfitexpr(expr2, mrk_flux_lores, mrk_clean_flux, mrk_clean_err, [1,1], $
	perror = mrk_lohi_perror, bestnorm = mrk_lohi_bestnorm, dof = mrk_lohi_dof)


!p.multi = [0,2,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg028.ps', /landscape
endif

; Fluxes

ploterror, iras_flux, iras_flux_true, iras_err, iras_err_true, $
	psym = 4, $
	title = 'IRAS 05189-2524', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 30], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oploterror, iras_flux_1p, iras_flux_true, iras_err_1p, iras_err_true, psym = 5, symsize = ss, thick = ct
oploterror, iras_flux_2p, iras_flux_true, iras_err_2p, iras_err_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*iras_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*iras_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*iras_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 3, 18, /data, 'Coadded slope: '+string(iras_result(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 3, 16, /data, 'Nod 1 slope: '+string(iras_result_1p(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 3, 14, /data, 'Nod 2 slope: '+string(iras_result_2p(0), format = '(f5.3)'), charsize = cs, charthick = ct

xyouts, fltarr(n_elements(iras_flux))+27, iras_flux_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

; EW

plot, iras_ew, iras_ew_true, $
	psym = 4, $
	title = 'IRAS 05189-2524', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 16], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, iras_ew_1p, iras_ew_true, psym = 5, symsize = ss, thick = ct
oplot, iras_ew_2p, iras_ew_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*iras_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*iras_ew_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*iras_ew_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 1, 11, /data, 'Coadded slope: '+string(iras_ew_result(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 1, 9, /data, 'Nod 1 slope: '+string(iras_ew_result_1p(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 1, 7, /data, 'Nod 2 slope: '+string(iras_ew_result_2p(0), format = '(f5.3)'), charsize = cs, charthick = ct

xyouts, fltarr(n_elements(iras_ew))+17, iras_ew_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif


; UGC 5101

; Fluxes

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg041.ps', /landscape
endif

ploterror, ugc_flux, ugc_flux_true, ugc_err, ugc_err_true, $
	psym = 4, $
	title = 'UGC 5101', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 55], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oploterror, ugc_flux_1p, ugc_flux_true, ugc_err_1p, ugc_err_true, psym = 5, symsize = ss, thick = ct
oploterror, ugc_flux_2p, ugc_flux_true, ugc_err_2p, ugc_err_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*ugc_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*ugc_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*ugc_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 3, 40, /data, 'Coadded slope: '+string(ugc_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 36, /data, 'Nod 1 slope: '+string(ugc_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 32, /data, 'Nod 2 slope: '+string(ugc_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(ugc_flux))+45, ugc_flux_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

; EW

plot, ugc_ew, ugc_ew_true, $
	psym = 4, $
	title = 'UGC 5101', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 120], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, ugc_ew_1p, ugc_ew_true, psym = 5, symsize = ss, thick = ct
oplot, ugc_ew_2p, ugc_ew_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*ugc_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*ugc_ew_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*ugc_ew_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 1, 65, /data, 'Coadded slope: '+string(ugc_ew_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 1, 60, /data, 'Nod 1 slope: '+string(ugc_ew_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 1, 55, /data, 'Nod 2 slope: '+string(ugc_ew_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(ugc_ew))+100, ugc_ew_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

; Mrk 273

; Fluxes

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg065.ps', /landscape
endif


ploterror, mrk_flux, mrk_flux_true, mrk_err, mrk_err_true, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 50], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oploterror, mrk_flux_1p, mrk_flux_true, mrk_err_1p, mrk_err_true, psym = 5, symsize = ss, thick = ct
oploterror, mrk_flux_2p, mrk_flux_true, mrk_err_2p, mrk_err_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*mrk_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*mrk_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*mrk_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 3, 50, /data, 'Coadded slope: '+string(mrk_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 46, /data, 'Nod 1 slope: '+string(mrk_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 42, /data, 'Nod 2 slope: '+string(mrk_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(mrk_flux))+45, mrk_flux_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

; EW

plot, mrk_ew, mrk_ew_true, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 120], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, mrk_ew_1p, mrk_ew_true, psym = 5, symsize = ss, thick = ct
oplot, mrk_ew_2p, mrk_ew_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*mrk_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*mrk_ew_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*mrk_ew_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 5, 120, /data, 'Coadded slope: '+string(mrk_ew_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 5, 110, /data, 'Nod 1 slope: '+string(mrk_ew_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 5, 100, /data, 'Nod 2 slope: '+string(mrk_ew_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(mrk_ew))+110, mrk_ew_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

print,''
print,'Uncleaned data fits'
print,''
print, 'IRAS coadded best fit: ', iras_result(0), ' +- ', iras_perror(0)
print, 'IRAS chi^2: ', iras_bestnorm / iras_dof
print, 'UGC coadded best fit: ', ugc_result(0), ' +- ', ugc_perror(0)
print, 'UGC chi^2: ', ugc_bestnorm / ugc_dof
print, 'Mrk coadded best fit: ', mrk_result(0), ' +- ', mrk_perror(0)
print, 'Mrk chi^2: ', mrk_bestnorm / mrk_dof
print,''

; Check to see if flux scaling has any wavelength dependence for coadded spectra

!p.multi = [0,1,1]
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg_ratio_wave.ps', /landscape
endif
plot, iras_flux / iras_flux_true, wave, /nodata,  $
	xtitle = 'Flux!Imeas!N / Flux!Ipub!N', $
	ytitle = 'Wavelength of line [!7l!3m]', $
	xrange = [0,2.5], $
	charsize = 1.5, charthick = ct, thick = ct
oplot, iras_flux / iras_flux_true, wave, psym = -4, linestyle = 0, symsize = ss, thick = ct
oplot, ugc_flux / ugc_flux_true,  wave, psym = -5, linestyle = 1, symsize = ss, thick = ct
oplot, mrk_flux / mrk_flux_true,  wave, psym = -6, linestyle = 2, symsize = ss, thick = ct
legend, ['IRAS','UGC','Mrk'], linestyle = [0,1,2], /top, /right, thick = ct, charthick = ct
if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

; ########################
; ##	Cleaned data	##
; ########################

!p.multi = [0,2,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg028_clean.ps', /landscape
endif

; Fluxes

ploterror, iras_clean_flux, iras_flux_true, iras_clean_err, iras_err_true, $
	psym = 4, $
	title = 'IRAS 05189-2524', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 30], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oploterror, iras_clean_flux_1p, iras_flux_true, iras_clean_err_1p, iras_err_true, psym = 5, symsize = ss, thick = ct
oploterror, iras_clean_flux_2p, iras_flux_true, iras_clean_err_2p, iras_err_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*iras_clean_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*iras_clean_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*iras_clean_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 3, 18, /data, 'Coadded slope: '+string(iras_result(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 3, 16, /data, 'Nod 1 slope: '+string(iras_result_1p(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 3, 14, /data, 'Nod 2 slope: '+string(iras_result_2p(0), format = '(f5.3)'), charsize = cs, charthick = ct

xyouts, fltarr(n_elements(iras_clean_flux))+27, iras_flux_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

; EW

plot, iras_clean_ew, iras_ew_true, $
	psym = 4, $
	title = 'IRAS 05189-2524', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 16], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, iras_clean_ew_1p, iras_ew_true, psym = 5, symsize = ss, thick = ct
oplot, iras_clean_ew_2p, iras_ew_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*iras_clean_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*iras_clean_ew_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*iras_clean_ew_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 1, 11, /data, 'Coadded slope: '+string(iras_ew_result(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 1, 9, /data, 'Nod 1 slope: '+string(iras_ew_result_1p(0), format = '(f5.3)'), charsize = cs, charthick = ct
;xyouts, 1, 7, /data, 'Nod 2 slope: '+string(iras_ew_result_2p(0), format = '(f5.3)'), charsize = cs, charthick = ct

xyouts, fltarr(n_elements(iras_clean_ew))+17, iras_ew_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif


; UGC 5101

; Fluxes

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg041_clean.ps', /landscape
endif

ploterror, ugc_clean_flux, ugc_flux_true, ugc_clean_err, ugc_err_true, $
	psym = 4, $
	title = 'UGC 5101', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 55], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oploterror, ugc_clean_flux_1p, ugc_flux_true, ugc_clean_err_1p, ugc_err_true, psym = 5, symsize = ss, thick = ct
oploterror, ugc_clean_flux_2p, ugc_flux_true, ugc_clean_err_2p, ugc_err_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*ugc_clean_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*ugc_clean_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*ugc_clean_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 3, 40, /data, 'Coadded slope: '+string(ugc_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 36, /data, 'Nod 1 slope: '+string(ugc_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 32, /data, 'Nod 2 slope: '+string(ugc_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(ugc_clean_flux))+45, ugc_flux_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

; EW

plot, ugc_clean_ew, ugc_ew_true, $
	psym = 4, $
	title = 'UGC 5101', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 120], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, ugc_clean_ew_1p, ugc_ew_true, psym = 5, symsize = ss, thick = ct
oplot, ugc_clean_ew_2p, ugc_ew_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*ugc_clean_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*ugc_clean_ew_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*ugc_clean_ew_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 1, 65, /data, 'Coadded slope: '+string(ugc_ew_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 1, 60, /data, 'Nod 1 slope: '+string(ugc_ew_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 1, 55, /data, 'Nod 2 slope: '+string(ugc_ew_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(ugc_clean_ew))+100, ugc_ew_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

; Mrk 273

; Fluxes

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg065_clean.ps', /landscape
endif


ploterror, mrk_clean_flux, mrk_flux_true, mrk_clean_err, mrk_err_true, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 50], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oploterror, mrk_clean_flux_1p, mrk_flux_true, mrk_clean_err_1p, mrk_err_true, psym = 5, symsize = ss, thick = ct
oploterror, mrk_clean_flux_2p, mrk_flux_true, mrk_clean_err_2p, mrk_err_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*mrk_clean_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*mrk_clean_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*mrk_clean_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 3, 50, /data, 'Coadded slope: '+string(mrk_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 46, /data, 'Nod 1 slope: '+string(mrk_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 3, 42, /data, 'Nod 2 slope: '+string(mrk_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(mrk_clean_flux))+45, mrk_flux_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

; EW

plot, mrk_clean_ew, mrk_ew_true, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 120], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, mrk_clean_ew_1p, mrk_ew_true, psym = 5, symsize = ss, thick = ct
oplot, mrk_clean_ew_2p, mrk_ew_true, psym = 6, symsize = ss, thick = ct

oplot, xarr, xarr*mrk_clean_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr*mrk_clean_ew_result_1p(0), linestyle = 2, thick = ct
oplot, xarr, xarr*mrk_clean_ew_result_2p(0), linestyle = 3, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct
;xyouts, 5, 120, /data, 'Coadded slope: '+string(mrk_ew_result(0), format = '(f5.3)'), charsize = cs
;xyouts, 5, 110, /data, 'Nod 1 slope: '+string(mrk_ew_result_1p(0), format = '(f5.3)'), charsize = cs
;xyouts, 5, 100, /data, 'Nod 2 slope: '+string(mrk_ew_result_2p(0), format = '(f5.3)'), charsize = cs

xyouts, fltarr(n_elements(mrk_clean_ew))+110, mrk_ew_true, obj, charsize = cs, charthick = ct

legend, ['Coadded', 'Nod 1', 'Nod 2'], linestyle = [1,2,3], /top, /left, charsize = cs, thick = ct, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

print,''
print,'Cleaned data fits'
print,''
print, 'IRAS coadded best fit: ', iras_clean_result(0), ' +- ', iras_clean_perror(0)
print, 'IRAS chi^2: ', iras_clean_bestnorm / iras_clean_dof
print, 'UGC coadded best fit: ', ugc_clean_result(0), ' +- ', ugc_clean_perror(0)
print, 'UGC chi^2: ', ugc_clean_bestnorm / ugc_dof
print, 'Mrk coadded best fit: ', mrk_clean_result(0), ' +- ', mrk_clean_perror(0)
print, 'Mrk chi^2: ', mrk_clean_bestnorm / mrk_clean_dof
print,''

; Check to see if flux scaling has any wavelength dependence for coadded spectra

!p.multi = [0,1,1]
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg_ratio_wave_clean.ps', /landscape
endif
plot, iras_clean_flux / iras_flux_true, wave, /nodata,  $
	xtitle = 'Flux!Imeas!N / Flux!Ipub!N', $
	ytitle = 'Wavelength of line [!7l!3m]', $
	xrange = [0,2.5], $
	charsize = 1.5, charthick = ct, thick = ct
oplot, iras_clean_flux / iras_flux_true, wave, psym = -4, linestyle = 0, symsize = ss, thick = ct
oplot, ugc_clean_flux / ugc_flux_true,  wave, psym = -5, linestyle = 1, symsize = ss, thick = ct
oplot, mrk_clean_flux / mrk_flux_true,  wave, psym = -6, linestyle = 2, symsize = ss, thick = ct
legend, ['IRAS','UGC','Mrk'], linestyle = [0,1,2], /top, /right, thick = ct, charthick = ct
if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

; Lores data comparison for Mrk 273

!p.multi = [0,2,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg065_lores.ps', /landscape
endif

; Fluxes

ploterror, mrk_flux_lores, mrk_flux_true, mrk_err_lores, mrk_err_true, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Published flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 60], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, xarr, xarr*mrk_lores_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct

xyouts, fltarr(n_elements(mrk_flux_lores))+27, mrk_flux_true, obj, charsize = cs, charthick = ct

; EW

plot, mrk_ew_lores, mrk_ew_true, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured EW [10!E-3!N !7l!3m]', $
	ytitle = 'Published EW [10!E-3!N !7l!3m]', $
	xrange = [0, 140], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

;oplot, xarr, xarr*mrk_lores_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct

xyouts, fltarr(n_elements(mrk_ew_lores))+17, mrk_ew_true, obj, charsize = cs, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

print,''
print,'Lores data fits'
print,''
print, 'Mrk coadded lo vs. true best fit: ', mrk_lores_result(0), ' +- ', mrk_lores_perror(0)
print, 'Mrk chi^2: ', mrk_lores_bestnorm / mrk_lores_dof
print,''

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = fpath+'ulirg065_lohi.ps', /landscape
endif

; Fluxes

ploterror, mrk_flux_lores, mrk_clean_flux, mrk_err_lores, mrk_clean_err, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured lo-res flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	ytitle = 'Measured hi-res flux [10!E-14!N erg cm!E-2!N s!E-1!N]', $
	xrange = [0, 60], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

oplot, xarr, xarr*mrk_lohi_result(1) + mrk_lohi_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct

xyouts, fltarr(n_elements(mrk_flux_lores))+27, mrk_clean_flux, obj, charsize = cs, charthick = ct

; EW

plot, mrk_ew_lores, mrk_clean_ew, $
	psym = 4, $
	title = 'Mrk 273', $
	xtitle = 'Measured lo-res EW [10!E-3!N !7l!3m]', $
	ytitle = 'Measured hi-res EW [10!E-3!N !7l!3m]', $
	xrange = [0, 140], $
	charsize = cs, $
	symsize = ss, thick = ct, charthick = ct

;oplot, xarr, xarr*mrk_lores_ew_result(0), linestyle = 1, thick = ct
oplot, xarr, xarr, linestyle = 0, thick = ct

xyouts, fltarr(n_elements(mrk_ew_lores))+17, mrk_clean_ew, obj, charsize = cs, charthick = ct

if keyword_set(ps) then begin
	device, /close
	set_plot, 'x'
endif

print, 'Mrk coadded lo vs. hi best fit:   ', mrk_lohi_result, ' +- ', mrk_lohi_perror
print, 'Mrk chi^2: ', mrk_lohi_bestnorm / mrk_lohi_dof
print,''



end
