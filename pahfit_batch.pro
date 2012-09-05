; Run PAHFIT on the OHM and control samples

device,window_state = state,decomposed=0

fnames = ohmdat('tag')
cnames = condat('tag')

fobj = ohmdat('obj')
cobj = condat('obj')

nf = n_elements(fnames)
nc = n_elements(cnames)

!p.multi=[0,1,1]
if state(2) eq 0 then window, 2, xsize = 800, ysize = 500 else wset,2

fpath = '~/Astronomy/Research/Spitzer/OHM/pahfit/'
for i = 0, nf - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/OHM/data/structures/'+fnames(i)+'.sav'

	; Add code to create report file if none exists (check)

	pahfile = fpath+sed.tag+'_pahfit.txt'
	findpahfile = findfile(pahfile)
	if findpahfile(0) eq '' then spawn,'touch '+pahfile

	fit = pahfit(sed.wave_lr, sed.flux_lr, sed.err_lr, redshift = 0.0, /plot_progress, $
	report = fpath+sed.tag+'_pahfit.txt', xsize = 800, ysize = 500, /no_megajansky_sr)

	set_plot,'ps'
	device,filename=fpath+'plots/'+fnames(i)+'.ps',/color,/landscape
	!p.thick = 2
	pahfit_plot,fit,sed.wave_lr,sed.flux_lr,sed.err_lr,units='Jy',title=fobj(i),symsize=0.4
	device,/close
	set_plot,'x'

	print,'Plotted SED for '+fobj(i)
endfor

cpath = '~/Astronomy/Research/Spitzer/control/pahfit/'
for i = 0, nc - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/control/data/structures/'+cnames(i)+'.sav'

	; Add code to create report file if none exists (check)

	pahfile = cpath+sed.tag+'_pahfit.txt'
	findpahfile = findfile(pahfile)
	if findpahfile(0) eq '' then spawn,'touch '+pahfile

	fit = pahfit(sed.wave_lr, sed.flux_lr, sed.err_lr, redshift = 0.0, $ ;/plot_progress, $
	report = cpath+sed.tag+'_pahfit.txt', xsize = 800, ysize = 500, /no_megajansky_sr)

	set_plot,'ps'
	device,filename=cpath+'plots/'+cnames(i)+'.ps',/color,/landscape
	!p.thick = 2
	pahfit_plot,fit,sed.wave_lr,sed.flux_lr,sed.err_lr,units='Jy',title=cobj(i),symsize=0.4
	device,/close
	set_plot,'x'

	print,'Plotted SED for '+cobj(i)
endfor






end
