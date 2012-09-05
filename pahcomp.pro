;pro pahcomp, NO_EXTINCTION=next,FP=fp,BLACK=black, SCALE_FAC=fac,WAVE_DIM=wd,_REF_EXTRA=e, stop = stop, ps = ps
;+
; NAME:
;       
;	PAHCOMP
;
; PURPOSE:
;
;	Plot an example of the PAHFIT fitting process vs. spline-fitting
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
;	PAHFIT must be compiled before compiling PAHCOMP
;
; REVISION HISTORY
;       Written by K. Willett                Dec 08
;-

fname = 'arch029'
tag, fname, dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/'+fname+'_fit.sav'
decoded_params = fit

restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
lambda = sed.wave_lr
intensity = sed.flux_lr
errors = sed.err_lr

units = 'Jy'

;if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/pahcomp.ps',/color,/portrait,xs=18,ys=8
	defcolor = fsc_color("Black")
	ls = 3
;endif else begin
;	defcolor = fsc_color("White")
;	ls = 1
;endelse

  !P.PSYM=0
;  device,get_decomposed=gd
;  device,decomposed=0
  tvlct,[255b,255b,0b,0b,160b,100b], $
        [0b,  0b,200b,0b,32b,100b], $
        [255b,0b,0b,255b,240b,100b],1
  if n_elements(units) eq 0 then units="MJy/sr"
  
  if !P.FONT eq 0 then begin 
     um='!Mm!Xm'
     nu='!Mn!X'
     lam='!Ml!X'
  endif else begin 
     um='!7l!Xm'
     nu='!7m!X'
     lam='!7k!X'
  endelse 
  
  rat=intensity/lambda
  rat = intensity
  if n_elements(fac) ne 0 then rat*=fac
  
  !p.multi = [0,1,1]
  !p.multi = [0,2,1]

  plot,lambda,rat,XTITLE='Wavelength (rest frame) ['+um+']', $
       YTITLE='Flux density [Jy]', $
       title = sed.obj, $
       xr=[5.5,7.0], $
       yr=[-0.0,0.5], $
       PSYM=10, XSTYLE=1, $
       YSTYLE=1, $
       _EXTRA=e, $
       charthick = 3, $
       thick = 3, $
       xthick = 3, $
       ythick = 3, $
       charsize = 1.0

;xyouts, 0.2, 0.8, 'PAHFIT', /normal, charsize = 1.5

  if n_elements(errors) ne 0 then begin 
     low=(intensity-errors)/lambda
     high=(intensity+errors)/lambda
     if n_elements(fac) ne 0 then begin 
        low*=fac
        high*=fac
     endif 
;     errplot,lambda,low,high
  endif 
  
  
  if keyword_set(fp) then lam=lambda else begin 
     if n_elements(wd) eq 0 then wd=1500
     mn=min(lambda,MAX=mx)
     lam=mn+findgen(wd)/(wd-1.)*(mx-mn)
  endelse 
  
  yfit=pahfit_components(lam,decoded_params,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot, $
                         EXTINCTION_FAC=ext,_EXTRA=e)
  
;  if ~keyword_set(next) then begin 
;     oplot,lam,!Y.CRANGE[0]+(!Y.CRANGE[1]-!Y.CRANGE[0])/1.05*ext, $
;           LINESTYLE=1,THICK=2,_EXTRA=e
;     axis,YAXIS=1,YRANGE=[0,1.05],/YSTYLE,YTITLE='Relative Extinction', $
;          CHARSIZE=1.3,_EXTRA=e
;  endif 
  
  rat=stars/lam
  rat = stars
  if n_elements(fac) ne 0 then rat*=fac
;  oplot,lam,rat,COLOR=1,THICK=2

  for i=0,(size(dusts,/DIMENSIONS))[1]-1 do begin 
     rat=dusts[*,i]/lam
     rat=dusts[*,i]
     if n_elements(fac) ne 0 then rat*=fac
;     oplot,lam,rat,COLOR=2,THICK=2,LINESTYLE=0,_EXTRA=e
  endfor 
  
  cont=dust_tot+stars
  
;  for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
;     rat=(cont+features[*,i])/lam
;     rat=(cont+features[*,i])
;     if n_elements(fac) ne 0 then rat*=fac
;     oplot,lam,rat,COLOR=4,LINESTYLE=0,_EXTRA=e
;  endfor 

   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])/lam
     rat=(cont+features[*,i])
     if n_elements(fac) ne 0 then rat*=fac
     oplot,lam,rat,COLOR=4,LINESTYLE=3,_EXTRA=e
  endfor 
  
     rat=(cont+features[*,2])/lam
     rat=(cont+features[*,2])
     if n_elements(fac) ne 0 then rat*=fac
     oplot,lam,rat,COLOR=4,LINESTYLE=3,_EXTRA=e, thick = 4

  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])/lam
     rat=(cont+lines[*,i])
     if n_elements(fac) ne 0 then rat*=fac
;     oplot,lam,rat,COLOR=5,LINESTYLE=0,_EXTRA=e
  endfor 
  
  rat=cont/lam
  rat=cont
  if n_elements(fac) ne 0 then rat*=fac  
  oplot,lam,rat,COLOR=fsc_color("Red"),THICK=4,_EXTRA=e
  
  rat=yfit/lam
  rat=yfit
  if n_elements(fac) ne 0 then rat*=fac  
  oplot,lam,rat,COLOR=fsc_color("Dark Grey"),_EXTRA=e, THICK=3, linestyle=2

;  device,decomposed=gd

limfac = 3.			; Factor by which to multiply positive PAH flux (for a conservative upper limit)
c = 299792.458d * 1d9		; um / s

; Read in the full LR spectrum

tag,fname,dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'

wave = lambda
flux = intensity
err = errors

; Spline fits for the 9.7 um continuum

		wcont = wave([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]) 
		fcont = flux([closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]) 
		yspl = spl_init(wcont, fcont)
		splresult = spl_interp(wcont, fcont, yspl, wave)
		pivotpoints = [closetomed(wave,flux,5.2),closetomed(wave,flux,5.6),$
			closetomed(wave,flux,14.0),closetomed(wave,flux,26.0)]

; Add feature to measure the PAH 6.2 um EW for placement in Spoon's classification scheme

; Designate continuum regions on either side

pivots_62 = [5.15, 5.95, 6.5, 7.1]

	npiv = n_elements(pivots_62)
	pahpiv = fltarr(npiv)
	for j = 0,npiv-1 do pahpiv[j] = closetomed(wave,flux,pivots_62[j])


bpah = closetomed(wave,flux,5.95)
epah = closetomed(wave,flux,6.55)

; Spline fit to PAH continuum

wp = wave(pahpiv)
fp = flux(pahpiv)
pspl1 = spl_init(wp,fp)
pahspl = spl_interp(wp,fp,pspl1,wave)

; Add the flux over designated area

newpah = flux - pahspl
addpah = fltarr(2,epah - bpah + 1)
for i = bpah, epah do begin
	addpah(0,i-bpah) = newpah(i)
	addpah(1,i-bpah) = (wave(i+1) - wave(i)) * c / (wave(i))^2
endfor

pah62fluxarr = addpah(0,*) * addpah(1,*) 					; Flux density in Jy, d_nu in Hz
pah62flux = total(pah62fluxarr) * 1d-30						; Flux in W/cm^2
cont62 = pahspl(closetomed(wave,flux,6.2))
pahew62 = pah62flux / cont62 * (6.2)^2 / c * 1d30				; EW in um
baseerr62 = stddev(flux(closetomed(wave,flux,6.1):closetomed(wave,flux,6.3)))	; Estimate error in baseline flux density [Jy]
pahew62err = baseerr62 * pah62flux / cont62^2 * (6.2)^2/ c*1d30		; Error in equivalent width

; 6.2 um FWHM

peak62 = max(newpah[bpah:epah])
peakwave62 = where(newpah[bpah:epah] eq peak62)
peakind62 = peakwave62 + bpah

i = 0 & nplevel = peak62
while nplevel gt 0.5 * peak62 do begin
	i = i+1
	nplevel = newpah[peakind62 - i]
endwhile

if i gt 0 then pah62fwhm = wave[peakind62 + i] - wave[peakind62 - i] else pah62fwhm = 0.

; Also define a PAH EW using the spline-based silicate continuum (correcting for water ice)

ice_62cont = splresult(closetomed(wave,splresult,6.2)) 
pahew62_ice = pah62flux / ice_62cont * (6.2)^2 / c * 1d30
baseerr_spl = stddev(splresult(closetomed(wave,splresult,6.1):closetomed(wave,splresult,6.3)))
pahew62_iceerr = baseerr_spl * pah62flux / ice_62cont^2 * (6.2)^2 / c*1d30

; Upper limits for PAH 6.2

if keyword_set(limit) then begin

	limrms = stddev(newpah[bpah:epah])
	fwhm62 = 0.14				; um - taken from average of all detected 6.2 um features in arch sample

	flux_cgs = limfac * limrms / (6.2)^2 * c * 1d-30
	sigwidth = fwhm62 / (2d * sqrt(2d * alog(2d)))
	pah62flux_lim = flux_cgs * sigwidth * sqrt(2d * !dpi) 			; Units of W/cm^2

	pahew62_lim = pah62flux_lim / cont62 * (6.2)^2 / c * 1d30
endif

; Plot the results

	; 6.2 um PAH pivot points, spline continuum, integrated area 
	
	oplot, wp, fp, psym = symcat(14), symsize = 1.5, color=fsc_color("Dark Green"), thick = ls
	oplot, wave, pahspl, color=fsc_color("Dark Green"), thick = ls, linestyle = 1
;	oplot, wave[bpah:epah], flux[bpah:epah], color=fsc_color("Dark Green"), thick = ls, linestyle = 0, psym=10
;	oplot, wave, splresult, color=fsc_color("Yellow"), thick = ls
;	ver, wave(bpah), linestyle = 1
;	ver, wave(epah), linestyle = 1
	
	;xyouts, 5.2, 0.07, 'PAH 6.2 um EW  : '+string(pahew62,format='(f6.3)')+' um', /data, charsize = cs
	;xyouts, 5.2, 0.05, 'PAH 6.2 um flux : '+string(pah62flux,format='(e9.3)')+' W/cm^2', /data, charsize = cs
	
plot, wave, flux, /nodata, $
        xr=[5.5,7.0], $
	yr=[-0.0,0.5], $
	charthick = 3, $
	thick = 3, $
	xthick = 3, $
	ythick = 3, $
	xtitle = 'Wavelength (rest frame) ['+um+']', $
	title = sed.obj

oplot, wave[bpah:epah-1], flux[bpah:epah-1] - pahspl[bpah:epah-1], color=fsc_color("Dark Green"),psym=10, thick = 3
;oplot, wp, intarr(n_elements(wp)), psym = symcat(14), symsize = 1, color=fsc_color("Dark Green"), thick = ls

     rat_pah=(features[*,2])
     oplot,lam,rat_pah,COLOR=4,LINESTYLE=3,_EXTRA=e, thick = 4

;if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
;endif

if keyword_set(stop) then stop
end


