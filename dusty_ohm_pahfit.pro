;+
; NAME:
;       
;	DUSTY_OHM_PAHFIT
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
;       Written by K. Willett                Jan 10
;-

pro dusty_ohm_pahfit, fname, lamfit, fluxfit, errors, noplot = noplot

tag, fname, dirtag

restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fname+'.sav'
lambda = sed.wave_lr
intensity = sed.flux_lr
errors = sed.err_lr

restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/pahfit/'+fname+'_fit.sav'
decoded_params = fit

     wd=1500
     mn=min(lambda,MAX=mx)
     lam=mn+findgen(wd)/(wd-1.)*(mx-mn)
  
  yfit=pahfit_components(lam,decoded_params,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot, $
                         EXTINCTION_FAC=ext,_EXTRA=e)
  
  lamfit = lam
  fluxfit = dust_tot + stars

if ~keyword_set(noplot) then begin

	plot, lambda, intensity, $
		charsize = 2, $
		/xlog, /ylog, $
		xr=[4,35], /xstyle, $$
		psym = 10, $
		xtitle='Wavelength', $
		ytitle='Flux density [Jy]'
	
	red = fsc_color("Red")
	
	oplot, lamfit, fluxfit, $
		thick = 2, $
		color=red

endif

end
