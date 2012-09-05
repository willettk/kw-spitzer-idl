pro sfr_cso, ps = ps, stop = stop
;+
; NAME:
;       
;	SFR_CSO
;
; PURPOSE:
;
;	Print SFR diagrams for Ne and PAH emission in CSO paper
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
;	IDL> sfr_cso, /ps
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Sep 08
;	Removed NGC 5793, 1245+676 - Nov 09
;-

mpc2cm = 3.086d24 
lsun   = 3.862d33 

; Restore neon data

restore,'~/Astronomy/Research/Spitzer/cso/linedata/neII.sav'
neII_cso = line
restore,'~/Astronomy/Research/Spitzer/cso/linedata/neIII.sav'
neIII_cso = line

neII_tag     = neII_cso.tag
neII_flux    = neII_cso.flux
neII_fluxerr = neII_cso.fluxerr

neIII_tag     = neIII_cso.tag
neIII_flux    = neIII_cso.flux
neIII_fluxerr = neIII_cso.fluxerr

taglist = csodat('tag')
lir = float(transpose(csodat('lir')))
lir = lir[where(finite(lir) eq 1)]
lir_tag = transpose(taglist)
lir_tag = lir_tag[where(finite(lir) eq 1)]

; Find data with both neon measurements and L_IR

match, neII_tag, neIII_tag, neII_ind, neIII_ind

neboth_tag = neII_tag[neII_ind]

neII_flux_temp     = neII_flux[neII_ind]
neII_fluxerr_temp  = neII_fluxerr[neII_ind]

neIII_flux_temp    = neIII_flux[neIII_ind]
neIII_fluxerr_temp = neIII_fluxerr[neIII_ind]

match, neboth_tag, lir_tag, neboth_ind, lirboth_ind

tag_both = neboth_tag[neboth_ind]

neII_flux_final     = neII_flux_temp[neboth_ind]
neII_fluxerr_final  = neII_fluxerr_temp[neboth_ind]

neIII_flux_final    = neIII_flux_temp[neboth_ind]
neIII_fluxerr_final = neIII_fluxerr_temp[neboth_ind]

lir_both = lir[lirboth_ind]

; Remove NGC 5793, 1245+676

goodboth = where(tag_both ne 'cso003' and tag_both ne 'cso010')

tag_both = neboth_tag[goodboth]
neII_flux_final     = neII_flux_temp[goodboth]
neII_fluxerr_final  = neII_fluxerr_temp[goodboth]
neIII_flux_final    = neIII_flux_temp[goodboth]
neIII_fluxerr_final = neIII_fluxerr_temp[goodboth]
lir_both = lir[goodboth]

; Find luminosity distances

dl = float(transpose(csodat('dl')))
dlerr = float(transpose(csodat('dlerr')))
match, transpose(taglist), tag_both, taglist_ind, tagboth_ind
dl_both = dl[taglist_ind]
dlerr_both = dlerr[taglist_ind]

; Compute luminosities and errors

netot_flux    = neII_flux_final + neIII_flux_final
netot_fluxerr = sqrt(neII_fluxerr_final^2 + neIII_fluxerr_final^2)

nelum = 4d * !dpi * (dl_both * mpc2cm)^2 * netot_flux * 1d7 / lsun
nelum_err = sqrt($
	dlerr_both^2    * (4d * !dpi * mpc2cm^2 * netot_flux * 2d * dl_both * 1d7 / lsun)^2 + $
	netot_fluxerr^2 * ((dl_both * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)^2 $
	)

nelum_log     = alog10(nelum)
nelum_log_err = alog10(nelum_err / nelum)

; Linear fit to Ne data

xarr = fillarr(1d-2,0,20)

expr='p[0]+x*p[1]'
start = [0,1]

fitboth = mpfitexpr(expr,lir_both, nelum_log, nelum_log_err,start,/quiet, perror = errboth)

; Upper limit on the neon flux for PKS 1413+135

lim = (linelim('cso005','neII',/noplot) + linelim('cso005','neIII',/noplot)) * 1d-21

nelumlim_cso005 = alog10(lim * 1d7 * 4d * !dpi * (dl[4] * mpc2cm)^2 / lsun)
lir_cso005 = lir[4]

; Plot results

red = fsc_color("Red")
green = fsc_color("Forest Green")

plotname='~/Astronomy/Research/Spitzer/cso/papers/sfr_cso.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color, /portrait
	cthick = 5
	lthick = 5
	cs = 2
	defcolor=fsc_color("Black")
endif else begin
	cs = 2
	cthick = 1
	lthick = 1
	defcolor=fsc_color("White")
endelse

!p.multi = [0,1,1]

plot, lir_both, $
	/nodata, $
	xtitle = 'log [L!IIR!N/L'+sunsymbol()+']', $
	ytitle = 'log [L!INe!N/L'+sunsymbol()+']', $
	charsize = cs, $
	xr = [9.0, 12.5], /xstyle, $
	yr = [6.0, 9.5], /ystyle, $
	color=defcolor, $
	thick = lthick, $
	xthick = lthick, $
	ythick = lthick, $
	charthick = cthick

;oplot, lir2, neII_lum, $
;	psym = symcat(16), $
;	color = red

oplot, lir_both, nelum_log, $
	psym = symcat(16), $
	color = defcolor, thick = lthick

oploterror, lir_both[0:2], nelum_log[0:2], nelum_err[0:2] / nelum[0:2], replicate(0.16,3), $
	/nohat, $
	psym = symcat(3), $
	color = defcolor, thick = lthick

oplot, [lir_cso005], [nelumlim_cso005], $
	psym = symcat(16), $
	color = defcolor, thick = lthick

arrow, lir_cso005, nelumlim_cso005, lir_cso005, nelumlim_cso005 - 0.5, $
	thick = lthick, color = defcolor, /data

; Fit to the CSOs (solid)

oplot, xarr, fitboth[0] + xarr*fitboth[1], color=defcolor, thick = lthick

; Ho and Keto (dotted)

oplot, xarr, xarr * 0.98 - 2.78, color = defcolor, linestyle = 1, thick = lthick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Plot the mean star formation rate using both Ne and PAH methods

; Neon

;nelumcgs = lsun * (10d)^neboth_lum
;sfr_neon = nelumcgs * 4.73d-41
;
;; PAH (spline)
;
;p6 = csodat('pah62lum')
;p11 = csodat('pah11lum')
;
;match, tag_list, ['cso005','cso010'], ca, cb
;ind = setdifference(indgen(n_elements(tag_list)),ca)
;p6 = p6[ind] & p11 = p11[ind]
;
;p6lum = lsun * (10d)^p6
;p11lum = lsun * (10d)^p11
;sfr_pah = (p6lum + p11lum) * 1.18d-41
;
;; PAH (PAHFIT)
;
;pf6 = csodat('pahfit62lum')
;pf11 = csodat('pahfit11lum')
;
;pf6 = pf6[ind] & pf11 = pf11[ind]
;
;pf6lum = lsun * (10d)^pf6
;pf11lum = lsun * (10d)^pf11
;sfr_pahfit = (pf6lum + pf11lum) * 1.18d-41
;
;print,'CSO SFR PAH spline = ',mean(sfr_pah), ' +_',stddev(sfr_pah)
;print,'CSO SFR neon =       ',mean(sfr_neon), ' +_',stddev(sfr_neon)
;print,'CSO SFR PAH PAHFIT = ',mean(sfr_pahfit), ' +_',stddev(sfr_pahfit)
;;print,'CSO NeII fit ',fit2, err2
;;print,'CSO Ne fit ',fitboth, errboth

if keyword_set(stop) then stop
end
