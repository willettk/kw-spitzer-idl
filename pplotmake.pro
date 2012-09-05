device,window_state = state,decomposed=0

!p.multi=[0,1,1]

fpath = '~/Astronomy/Comps2/images/'

	restore,'~/Astronomy/Research/Spitzer/OHM/data/structures/mega027.sav'

	; Add code to create report file if none exists (check)

	pahfile = fpath+sed.tag+'_pahfit.txt'
	findpahfile = findfile(pahfile)
	if findpahfile(0) eq '' then spawn,'touch '+pahfile

	fit = pahfit(sed.wave_lr, sed.flux_lr, sed.err_lr, redshift = 0.0, $
	xsize = 800, ysize = 500, /no_megajansky_sr)

	set_plot,'ps'
	device,filename=fpath+'mega027.ps',/color,/landscape
	!p.thick = 2
	pahfit_plot,fit,sed.wave_lr,sed.flux_lr,sed.err_lr,units='Jy',title=sed.obj,symsize=0.4
	device,/close
	set_plot,'x'



end
