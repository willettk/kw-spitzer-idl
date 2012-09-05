pro mypahfitplot,fname, $
                NO_EXTINCTION=next,FP=fp,BLACK=black, $
                SCALE_FAC=fac,WAVE_DIM=wd,_REF_EXTRA=e, $
		stop = stop, ps = ps
;+
; NAME:
;
;    PAHFIT_PLOT
;
; PURPOSE:
;
;    Plot the results of a PAHFIT decomposition.
;
; CALLING SEQUENCE:
;
;    pahfit_plot, pahfit_fit, lambda, flux, [flux_uncertainty,
;       UNITS=, /NO_EXTINCTION, /FAST_PLOT, /BLACK, SCALE_FAC=, _EXTRA=]
;
; INPUTS:
;
;    pahfit_fit: The decoded fit parameter structure returned by
;       PAHFIT.
;    
;    lambda: The wavelength in microns, in the rest frame.
;
;    flux: The flux (or flux intensity), in whatever units is specified
;       by the UNITS keyword (default: MJy/sr).
;
; OPTIONAL INPUTS:
;
;    flux_uncertainty: The flux intensity uncertainty, in the same
;       units as FLUX.
;
; KEYWORD PARAMETERS:
;
;    UNITS: The units to use in the Y axis title (default MJy/sr).
;
;    NO_EXTINCTION: Do not plot the extinction effect and "Relative
;       Extinction" right axis.
;
;    FAST_PLOT: Plot the model at the same sampling as the input
;       LAMBDA vector, rather than on a smoother grid.
;
;    BLACK: Plot the model as black instead of green.
;
;    SCALE_FAC: A factor by which to scale the model (for instance to
;       change the units).
;    
;    _EXTRA: Any extra plotting keyword parameters.
;      
; RESTRICTIONS:
;
;   The passed spectrum must be in the rest frame.
;   
; EXAMPLE:
;
;   IDL> pahfit_plot,lam/(1+cz/3.e5),flux,flux_unc
;
; MODIFICATION HISTORY:
;
; Adapted from PAHFIT_PLOT			- KW, Sep 08
;	Properly plots the difference spectrum - KW, Jan 10
;
;-

tag, fname, dirtag

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/'+fname+'_fit.sav'
decoded_params = fit

restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
lambda = sed.wave_lr
intensity = sed.flux_lr
errors = sed.err_lr

units = 'Jy'

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/plots/'+fname+'_decomp.ps',/color,/landscape
	defcolor = fsc_color("Black")
endif else defcolor = fsc_color("White")

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
  
  !p.multi = [0,1,2]

  plot,lambda,rat,XTITLE='Wavelength ('+um+')', $
       YTITLE='Flux density [Jy]', $
       PSYM=10, XSTYLE=1, $
       YSTYLE=1, $
       _EXTRA=e, $
       charthick = 2, $
       thick = 2, $
       charsize = 1.5

xyouts, 0.2, 0.8, 'PAHFIT', /normal, charsize = 1.5

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
  oplot,lam,rat,COLOR=1,THICK=2
  for i=0,(size(dusts,/DIMENSIONS))[1]-1 do begin 
     rat=dusts[*,i]/lam
     rat=dusts[*,i]
     if n_elements(fac) ne 0 then rat*=fac
     oplot,lam,rat,COLOR=2,THICK=2,LINESTYLE=0,_EXTRA=e
  endfor 
  
  cont=dust_tot+stars
  
  for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])/lam
     rat=(cont+features[*,i])
     if n_elements(fac) ne 0 then rat*=fac
     oplot,lam,rat,COLOR=4,LINESTYLE=0,_EXTRA=e
  endfor 
  
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])/lam
     rat=(cont+lines[*,i])
     if n_elements(fac) ne 0 then rat*=fac
     oplot,lam,rat,COLOR=5,LINESTYLE=0,_EXTRA=e
  endfor 
  
  rat=cont/lam
  rat=cont
  if n_elements(fac) ne 0 then rat*=fac  
  oplot,lam,rat,COLOR=6,THICK=4,_EXTRA=e			; Grey
  
  rat=yfit/lam
  rat=yfit
  if n_elements(fac) ne 0 then rat*=fac  
  oplot,lam,rat,THICK=keyword_set(black)?2:3,COLOR=keyword_set(black)?0:3, $
        _EXTRA=e
;  device,decomposed=gd

; Plot difference spectrum

!p.thick = 1
!p.charthick = 1

plot,lambda,intensity,XTITLE='Wavelength ('+um+')', $
       YTITLE='Flux density [Jy]', $
       PSYM=10, XSTYLE=1, $
       YSTYLE=1, $
       _EXTRA=e, $
       charthick = 2, $
       thick = 2, $
       charsize = 1.5, $
       yr=[-0.1,1.5], $
       /nodata

oplot, lambda, intensity, color=defcolor, psym = 10, thick = 2
oplot, lam, yfit, color=fsc_color("Dark Green")

;oplot, lambda, yfit - features_tot, color=fsc_color("Orange")		; everything but PAH
;oplot, lambda, features_tot, color=fsc_color("Red")			; Total PAH
oplot, lambda, yfit - cont, color=fsc_color("Orange Red")		; Total PAH
oplot, lambda, features[*,2]+cont, color=fsc_color("Purple"),thick=3	; Continuum
;oplot, lambda, features[*,2] + max(intensity[45:70]) - max(features[45:70,2]), color=fsc_color("Pink"),thick=1	; Shifted PAH feature

yfitdiff = fltarr(n_elements(lambda))
for i=0,n_elements(lambda)-1 do begin
	yfitdiff[i] = yfit[closeto(lam,lambda[i])]
endfor

oplot, lambda, (intensity - yfitdiff) , color=fsc_color("Blue"), psym = 10	; Differenced spectrum of data and total fit

ver, 5.95,linestyle=1
ver, 6.55,linestyle=1

legend,/top,/left, $
	['Data','Global fit','6.2 PAH + cont.','Fit - cont.','Diff. spectrum'], $
	color = [defcolor,fsc_color("Dark Green"),fsc_color("Purple"),fsc_color("Orange Red"),fsc_color("Blue")], $
	linestyle = intarr(5), thick = 3, charthick = 1, charsize = 1.6

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

!p.multi=[0,1,1]

if keyword_set(stop) then stop
end
