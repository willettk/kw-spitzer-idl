pro sfr_compare, ps=ps, stop=stop
;+
; NAME:
;       
;	SFR_COMPARE
;
; PURPOSE:
;
;	Compare star formation rate indicators for all CSOs
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
;       Written by K. Willett                Sep 09
;-

mpc2cm = 3.086d24 
lsun   = 3.862d33 

obj = csodat('obj')

; IR data

cso_iras_string = float(csodat('lir'))
cso_iras = transpose(10.^cso_iras_string * lsun)		; L_IR in erg/s

; PAH

pah62 = float(csodat('pah62lum'))
pah11 = float(csodat('pah11lum'))

pah_all = transpose(10d^(pah62) + 10d^(pah11)) * lsun	; L_PAH in erg/s

; Neon

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

; Pad for non-detections in PKS 1413, VII ZW in neon

neII_tag = [neII_tag[0:3],'cso005',neII_tag[4:7],'cso010']
neII_flux = [neII_flux[0:3],0,neII_flux[4:7],0]
neII_fluxerr = [neII_fluxerr[0:3],0,neII_fluxerr[4:7],0]

neIII_tag = [neIII_tag[0:3],'cso005',neIII_tag[4:7],'cso010']
neIII_flux = [neIII_flux[0:3],0,neIII_flux[4:7],0]
neIII_fluxerr = [neIII_fluxerr[0:3],0,neIII_fluxerr[4:7],0]

; Find luminosity distances

dl = float(transpose(csodat('dl')))
dlerr = float(transpose(csodat('dlerr')))

; Compute luminosities and errors

netot_flux    = neII_flux + neIII_flux
netot_fluxerr = sqrt(neII_fluxerr^2 + neIII_fluxerr^2)

neon = 4d * !dpi * (dl * mpc2cm)^2 * netot_flux * 1d7
neon_err = sqrt($
	dlerr^2    * (4d * !dpi * mpc2cm^2 * netot_flux * 2d * dl * 1d7)^2 + $
	netot_fluxerr^2 * ((dl * mpc2cm)^2 * 4d * !dpi * 1d7)^2 $
	)

sfr_neon = 4.73d-41 * neon

sfr_pah = 1.18d-41 * pah_all

sfr_iras = 4.5d-44 * cso_iras

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/cso/papers/sfr_compare.ps',/color,/portrait
	defcolor = fsc_color("Black")
	cs = 1
	lthick = 5
	cthick = 3
	hsize = 200
endif else begin
	defcolor=fsc_color("White")
	cs = 2
	lthick = 1
	cthick = 1
	hsize = 10
endelse

ne_pah_ind = [0,1,3,5,8]
ne_iras_ind = [0,1,3,5]
pah_iras_ind = [0,1,3,5]

plot, alog10(sfr_iras[ne_iras_ind]), alog10(sfr_neon[ne_iras_ind]), $
	/nodata, $
;	/xlog, /ylog, $
	xr = [-1,3.5], /xstyle, $
	yr = [-1,3], /ystyle, $
	charsize = cs, $
	charthick = cthick, $
	thick=lthick, $
	xthick=lthick, $
	ythick=lthick, $
	xtitle = 'log (SFR!IIRAS!N/M'+sunsymbol()+'/yr)', $
	ytitle = 'log (SFR/M'+sunsymbol()+'/yr)'

oplot, alog10(sfr_iras[ne_iras_ind]), alog10(sfr_neon[ne_iras_ind]), psym=symcat(16)
oplot, findgen(10) - 2, findgen(10) - 2, linestyle=1, thick = lthick
oplot, alog10(sfr_iras[pah_iras_ind]), alog10(sfr_pah[pah_iras_ind]), psym=symcat(46), color=fsc_color("Red")

xyouts, alog10(sfr_iras[ne_iras_ind]) , alog10(sfr_neon[ne_iras_ind]) + 0.1, obj[ne_pah_ind], charthick = cthick

; Upper limit on neon, PAH flux for PKS 1413+135

lim = (linelim('cso005','neII',/noplot) + linelim('cso005','neIII',/noplot)) * 1d-21
nelumlim_cso005 = lim * 1d7 * 4d * !dpi * (dl[4] * mpc2cm)^2 * 4.73d-41
cso005_pah = (10^9.05 + 10^8.76) * lsun * 1.18d-41

oplot, [alog10(sfr_iras[4])], [alog10(nelumlim_cso005)], psym=symcat(16)
arrow, [alog10(sfr_iras[4])], [alog10(nelumlim_cso005)], [alog10(sfr_iras[4])], [alog10(nelumlim_cso005)] - 0.2, /data
oplot, [alog10(sfr_iras[4])], [alog10(cso005_pah)], psym=symcat(46), color=fsc_color("Red")
arrow, [alog10(sfr_iras[4])], [alog10(cso005_pah)], [alog10(sfr_iras[4])], [alog10(cso005_pah)] - 0.2, /data, color=fsc_color("Red")
xyouts, [alog10(sfr_iras[4])] + 0.1, [alog10(nelumlim_cso005)] - 0.2, [obj[4]], charthick = cthick

legend, /top, /left, color=[defcolor,fsc_color("Red")], psym=[16,46], ['Neon','PAH'], thick=lthick, charthick = cthick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

;plot, sfr_neon[ne_iras_ind], sfr_iras[ne_iras_ind], $
;	/ylog, $
;	/xlog, $
;	charsize = cs, $
;	psym = symcat(16), $
;	xtitle = 'SFR!INe!N', $
;	ytitle = 'SFR!IIRAS!N'
;
;oplot, findgen(10000)/10, findgen(10000)/10, linestyle=1
;xyouts, sfr_neon[ne_iras_ind], sfr_iras[ne_iras_ind], obj[ne_iras_ind]
;
;plot, sfr_pah[pah_iras_ind], sfr_iras[pah_iras_ind], $
;	/ylog, $
;	/xlog, $
;	charsize = cs, $
;	psym = symcat(16), $
;	xtitle = 'SFR!IPAH!N', $
;	ytitle = 'SFR!IIRAS!N'
;
;oplot, findgen(10000)/10, findgen(10000)/10, linestyle=1
;xyouts, sfr_pah[pah_iras_ind], sfr_iras[pah_iras_ind], obj[pah_iras_ind]

stop

;!p.multi=[0,2,1]
;cs = 2
;
;ne_pah_ind = [0,1,3,5,8]
;ne_iras_ind = [0,1,3,5]
;pah_iras_ind = [0,1,3,5]
;
;plot, sfr_iras[ne_iras_ind], sfr_neon[ne_iras_ind], $
;	/ylog, $
;	/xlog, $
;	charsize = cs, $
;	psym = symcat(16), $
;	ytitle = 'SFR!INe!N', $
;	xtitle = 'SFR!IIRAS!N'
;
;oplot, findgen(10000)/10, findgen(10000)/10, linestyle=1
;xyouts, sfr_iras[ne_iras_ind], sfr_neon[ne_iras_ind], obj[ne_iras_ind]
;
;plot, cso_iras_string[ne_iras_ind], alog10(neon[ne_iras_ind]/lsun), $
;	xr=[9,12.5], $
;	yr=[6,9.5], $
;	/xstyle, /ystyle, $
;	charsize=cs, $
;	psym=symcat(16), $
;	ytitle = 'log L!INe!N', $
;	xtitle = 'log L!IIR!N'
;
;oplot, [cso_iras_string[4]], [alog10(10^8.44+10^8.76)], $
;	psym = symcat(16)
;arrow, [cso_iras_string[4]], [alog10(10^8.44+10^8.76)], [cso_iras_string[4]], [alog10(10^8.44+10^8.76)] - 1, /data
;
;!p.multi=[0,1,1]

if keyword_set(stop) then stop

end
