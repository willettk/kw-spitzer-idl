;+
; NAME:
;       
;	CSO_NEPAH
;
; PURPOSE:
;
;	
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
;
; REVISION HISTORY
;       Written by K. Willett                Dec 08
;	Removed NGC 5793, 1245+676 - Nov 09
;	Fixed incorrect indices, added LABEL keyword - Nov 09
;-

pro cso_nepah, stop=stop, ps = ps, label=label

both_det_ind = [0,1,3,5,8]
;both_det_ind = [0,1,2,3,5,8]

; PAH

pah62 = float(csodat('pah62lum'))
pah11 = float(csodat('pah11lum'))

pah_all = alog10(10d^(pah62) + 10d^(pah11))
pah_both = pah_all[both_det_ind]

; Neon

mpc2cm = 3.086d24 
lsun   = 3.862d33 

; Restore neon data

restore,'~/Astronomy/Research/Spitzer/cso/linedata/neII.sav'
neII_cso = line
restore,'~/Astronomy/Research/Spitzer/cso/linedata/neIII.sav'
neIII_cso = line

;taglist_ind = [0,1,2,3,5,6,7,8]
taglist_ind = [0,1,3,4,7]

neII_tag     = neII_cso.tag
neII_flux    = neII_cso.flux
neII_fluxerr = neII_cso.fluxerr

neIII_tag     = neIII_cso.tag
neIII_flux    = neIII_cso.flux
neIII_fluxerr = neIII_cso.fluxerr

ne_tag = neII_tag[taglist_ind]

neII_flux_both     = neII_flux[taglist_ind]
neII_fluxerr_both  = neII_fluxerr[taglist_ind]

neIII_flux_both    = neIII_flux[taglist_ind]
neIII_fluxerr_both = neIII_fluxerr[taglist_ind]

; Find luminosity distances for objects with both neon and PAH

dl = float(transpose(csodat('dl')))
dlerr = float(transpose(csodat('dlerr')))
dlind = [1,2,4,6,9] - 1
dl_both = dl[dlind]
dlerr_both = dlerr[dlind]

; Compute luminosities and errors

netot_flux    = neII_flux_both + neIII_flux_both
netot_fluxerr = sqrt(neII_fluxerr_both^2 + neIII_fluxerr_both^2)

neon = 4d * !dpi * (dl_both * mpc2cm)^2 * netot_flux * 1d7 / lsun
neon_err = sqrt($
	dlerr_both^2    * (4d * !dpi * mpc2cm^2 * netot_flux * 2d * dl_both * 1d7 / lsun)^2 + $
	netot_fluxerr^2 * ((dl_both * mpc2cm)^2 * 4d * !dpi * 1d7 / lsun)^2 $
	)

nelum_both = alog10(neon)
nelum_both_err = alog10(neon_err)

; PKS 1413+135

neii_cso005  = linelim('cso005','neii',/noplot)
neiii_cso005 = linelim('cso005','neiii',/noplot)

cso005_ne = alog10((neii_cso005 + neiii_cso005) * 4d * !dpi * (dl[4] * 3.086d24)^2 * 1d-21 * 1d7 / 3.862d33)
cso005_pah = alog10(10^9.05 + 10^8.76)

; 1946+70 and 4C 37.11

cso007_ind = where(neII_tag eq 'cso007')
cso007_ne = alog10(4d * !dpi * (dl[6] * mpc2cm)^2 * (neII_flux[cso007_ind] + neIII_flux[cso007_ind]) * 1d7 / lsun)
cso007_pah = pah11[6]

cso008_ind = where(neII_tag eq 'cso008')
cso008_ne = alog10(4d * !dpi * (dl[7] * mpc2cm)^2 * (neII_flux[cso008_ind] + neIII_flux[cso008_ind]) * 1d7 / lsun)
cso008_pah = pah11[7]

; Linear fit to Ne data

xarr = fillarr(0.1,0,30)
expr='p[0]+x*p[1]'
start = [0,1]

fit = mpfitexpr(expr,pah_both,nelum_both,nelum_both * 0.1,start,/quiet, perror = err)

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/cso/papers/cso_nepah.ps', /portrait
	th = 5
	cs = 2.0
	hsize = 300
endif else begin
	th = 1
	cs = 1.0
	hsize = 10
endelse

plot, pah_both, nelum_both, $
	ytitle = 'log (L!INe!N/L'+sunsymbol()+')', $
	xtitle = 'log (L!IPAH!N/L'+sunsymbol()+')', $
	xr = [6.5,9.5], /xstyle, $
	yr = [6,9.5], /ystyle, $
	thick = th, $
	xthick = th, $
	ythick = th, $
	charthick = th, $
	charsize = cs, $
	psym = symcat(16)

oploterror, pah_both, nelum_both, neon_err / neon, replicate(0.16,n_elements(pah_both)), $
	/nohat, $
	psym = symcat(3), $
	thick = th


; Limits on PKS 1413+135 and VII Zw 485

oplot, [cso005_pah], [cso005_ne], psym=symcat(16)
arrow, [cso005_pah], [cso005_ne], [cso005_pah],               [cso005_ne] - 0.3, /data, thick = th, hsize = hsize
arrow, [cso005_pah], [cso005_ne], [cso005_pah]- 0.3 * 8.5/11, [cso005_ne],       /data, thick = th, hsize = hsize
;oplot, [cso005_pah,cso010_pah], [cso005_ne,cso010_ne], psym=symcat(16)
;arrow, [cso005_pah,cso010_pah], [cso005_ne,cso010_ne], [cso005_pah,cso010_pah],               [cso005_ne,cso010_ne] - 0.5, /data, thick = th, hsize = 300
;arrow, [cso005_pah,cso010_pah], [cso005_ne,cso010_ne], [cso005_pah,cso010_pah]- 0.5 * 8.5/11, [cso005_ne,cso010_ne],       /data, thick = th, hsize = 300

; Limits on 1946+70 and 4C 37.11

oplot, [cso007_pah,cso008_pah], [cso007_ne,cso008_ne], psym=symcat(16)
arrow, [cso007_pah,cso008_pah], [cso007_ne,cso008_ne], [cso007_pah,cso008_pah]- 0.5 * 8.5/11, [cso007_ne,cso008_ne], $
	/data, thick = th, hsize = hsize

; Label objects

if keyword_set(label) then begin
	obj = csodat('obj',/ver)
	match,transpose(strtrim(obj[0,*],2)), strtrim(ne_tag,2), i1, i2
	xyouts, pah_both + 0.1, nelum_both, obj[1,i1]
	xyouts, cso007_pah + 0.1, cso007_ne, '1946+70'
	xyouts, cso008_pah + 0.1, cso008_ne, '4C +37.11'
	xyouts, cso005_pah + 0.1, cso005_ne, 'PKS 1413+135'
endif

; Overplot fit from Farrah ULIRGs

oplot, xarr, alog10(0.17) + 1.02 * xarr, linestyle = 2, thick = th

; Overplot my best fit

oplot, xarr, fit[0] + fit[1] * xarr, linestyle = 0, thick = th

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

print,''
print,'Correlation factor ',correlate(pah_both, nelum_both)
print,' Fit: ',fit
print,' Error: ',err
print,''

if keyword_set(stop) then stop

end
