;pro zspitz_linestats, stop=stop
;+
; NAME:
;       
;	ZSPITZ_LINESTATS
;
; PURPOSE:
;
;	Compute stats on each mid-IR hi-res line for the measured IR redshifts from ZSPITZER	
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
;	IDL> .r zspitz_linestats
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Apr 08
;	Added archived OHMs w/measured z_OH - Oct 08
;-

c = 299792.458d ; km/s

; Darling OHMs

ohmlist = transpose(ohmdat('tag'))
ohmlist = ohmlist[where(ohmlist ne 'mega016')]

; Archived OHMs

archlist = transpose(archdat('tag'))

; Combine OHM samples

objlist = [ohmlist,archlist]
nobj = n_elements(objlist)

; Directories w/saved information

savdir_ohm = '~/Astronomy/Research/Spitzer/ohm/lines/hrsky/hires/saved/'
savdir_arch = '~/Astronomy/Research/Spitzer/archived/lines/hrsky/hires/saved/'
savdir_arch_nosky = '~/Astronomy/Research/Spitzer/archived/lines/nosky/hires/saved/'

; Archived OHMs with HR sky backgrounds

hrarch = 'arch'+string([5,10,20,29,31,34,35],format='(i03)')

; Ions to plot

ionlist = ['SIII', $
	'SIV', $
	'OIV', $
	'ArIII', $
	'H2s0', $
	'H2s1', $
	'H2s2', $
	'H2s3', $
	'NeII', $
	'NeIII', $
	'NeV', $
	'HI76']

!p.multi = [0,3,4]

; Retrieve individual IR line velocities for each ion

for k = 0, n_elements(ionlist)-1 do begin				

	; Empty arrays

	ionzarr = intarr(nobj)
	czoh_arr = intarr(nobj)
	dvel_ir_avgir = dblarr(nobj)
	dvel_ir_opt = dblarr(nobj)
	dvel_ir_ohm = dblarr(nobj)
	
	; Loop over objects

	for i = 0, nobj-1 do begin	

		; Find data directory

		listtag = strmid(objlist[i],0,4)
		if listtag eq 'mega' then savdir = savdir_ohm $
			else if listtag eq 'arch' and where(objlist[i] eq hrarch) ne -1 then savdir = savdir_arch $
			else savdir = savdir_arch_nosky

		; Restore data and test if ion is present

		restore,savdir+objlist[i]+'_zlines.sav'
		ion_now = ionlist[k]
		isionthere = where(ion eq ion_now)

		; Retrieve difference of measured ion velocity with IR average, optical, and OHM velocities 

		if isionthere ge 0 then begin

			; Increment for detection of ion in this object

			ionzarr[i] = 1

			dvel_ir_avgir[i] = velir_arr[isionthere] - wavg_vel
			dvel_ir_opt[i]   = velir_arr[isionthere] - cz_opt
			dvel_ir_ohm[i]   = velir_arr[isionthere] - czoh

		endif

		; Increment for detection of OHM velocity in this object

		if czoh gt 0. then czoh_arr[i] = 1

	endfor
	
	; Remove data where ion was not detected or OHM velocity is unknown

	remove_zeroes = where(ionzarr eq 1 and czoh_arr eq 1)
	dvel_ir_avgir = dvel_ir_avgir[remove_zeroes]
	dvel_ir_opt   = dvel_ir_opt[remove_zeroes]
	dvel_ir_ohm   = dvel_ir_ohm[remove_zeroes]

	; Print results to screen
	
	print,''
	print,'Average vel. difference of '+ionlist[k]+' from avg. IR: ',$
		string(mean(dvel_ir_avgir),format='(f7.1)'),' +- ',string(stddev(dvel_ir_avgir),format='(f7.1)'), ' km/s'
	print,'Average vel. difference of '+ionlist[k]+' from optical: ',$
		string(mean(dvel_ir_opt),format='(f7.1)'),' +- ',string(stddev(dvel_ir_opt),format='(f7.1)'), ' km/s'
	print,'Average vel. difference of '+ionlist[k]+' from OHM:     ',$
		string(mean(dvel_ir_ohm),format='(f7.1)'),' +- ',string(stddev(dvel_ir_ohm),format='(f7.1)'), ' km/s'
	
	; Scatter plot velocity trend for each species
	
	plot, dvel_ir_avgir, $
		/nodata, $
		psym = symcat(16), $
		xtitle = '!7D!3v!Iopt!N  [km/s]', $
		ytitle = '!7D!3v!IOHM!N [km/s]', $
		title = ionlist[k], $
		yrange=[-500,500], /ystyle, $
		charsize = 1.5
;		xrange = [-1000,1000], /xstyle, $
;		yrange = [-1000,1000], /ystyle
	
	oplot,dvel_ir_avgir, psym=symcat(16), symsize=1		; Optical vs. OH

;	oplot, fillarr(1,-1000,1000), fillarr(1,-1000,1000), linestyle=0

	hor, mean(dvel_ir_avgir), linestyle=2
	hor, mean(dvel_ir_avgir) - stddev(dvel_ir_avgir), linestyle=1
	hor, mean(dvel_ir_avgir) + stddev(dvel_ir_avgir), linestyle=1
	hor, 0, linestyle=0
	
;	ver, 0, linestyle=1
;	hor, 0, linestyle=1

endfor

; Obtain the bulk data for each object

vel_ir = dblarr(nobj)  & vel_ir_err = dblarr(nobj)
vel_opt = dblarr(nobj) & vel_opt_err = dblarr(nobj)
vel_oh = dblarr(nobj)  & vel_oh_err = dblarr(nobj)

for i = 0, nobj - 1 do begin

	listtag = strmid(objlist[i],0,4)
	if listtag eq 'mega' then savdir = savdir_ohm $
		else if listtag eq 'arch' and where(objlist[i] eq hrarch) ne -1 then savdir = savdir_arch $
		else savdir = savdir_arch_nosky

	restore,savdir+objlist[i]+'_zlines.sav'

	vel_ir[i]  = wavg_vel & vel_ir_err[i]  = wstddev_vel
	vel_opt[i] = cz_opt   & vel_opt_err[i] = cz_opt_err
	vel_oh[i]  = czoh     & vel_oh_err[i]  = czoh_err

endfor

; Remove velocities w/out measurements in all sets of velocity differences

clipped_indices = where(vel_ir ne 0 and vel_opt ne 0 and vel_oh ne 0 and vel_ir_err lt 1d4)
vel_iravg_clipped = vel_ir[clipped_indices]  & vel_iravg_clipped_err = vel_ir_err[clipped_indices]
vel_opt_clipped   = vel_opt[clipped_indices] & vel_opt_clipped_err   = vel_opt_err[clipped_indices]
vel_oh_clipped    = vel_oh[clipped_indices]  & vel_oh_clipped_err    = vel_oh_err[clipped_indices]

; Add errors for velocity differences in quadrature

vel_optir_err = sqrt(vel_iravg_clipped_err^2 + vel_opt_clipped_err^2)
vel_ohir_err  = sqrt(vel_iravg_clipped_err^2 + vel_oh_clipped_err^2)
vel_ohopt_err = sqrt(vel_opt_clipped_err^2 + vel_opt_clipped_err^2)

; Compute the weighted mean and uncertainty for all sets of velocity differences

wtdmean, vel_opt_clipped - vel_iravg_clipped, vel_optir_err, mean_optir_veldiff, sig_optir_veldiff
wtdmean, vel_oh_clipped - vel_iravg_clipped, vel_ohir_err, mean_ohir_veldiff, sig_ohir_veldiff
wtdmean, vel_oh_clipped - vel_opt_clipped, vel_ohopt_err, mean_ohopt_veldiff, sig_ohopt_veldiff

; Whisker plot of the difference between two measurements (ie, mean(IR) - optical) vs. another (ie, mean(IR))

!p.multi = [0,1,1]
cs = 2
set_plot,'ps'

device,filename='~/Astronomy/Research/Spitzer/plots/zspitzer_ir_opt.ps', /landscape,/color
ploterror, vel_opt_clipped - vel_iravg_clipped, vel_iravg_clipped, $
	sqrt(vel_iravg_clipped_err^2 + vel_opt_clipped_err^2), vel_iravg_clipped_err, $
	xtitle='v!Iopt!N - v!IIR!N [km/s]', $
	ytitle='v!IIR!N [km/s]', $
	charsize=cs, $
	psym = symcat(16)

ver,mean_optir_veldiff,linestyle=0,thick=4
;ver,mean_optir_veldiff - sig_optir_veldiff,linestyle=2, thick=3
;ver,mean_optir_veldiff + sig_optir_veldiff,linestyle=2, thick=3
device,/close

device,filename='~/Astronomy/Research/Spitzer/papers/zspitzer_ir_oh.ps', /portrait
ploterror, vel_oh_clipped - vel_iravg_clipped, vel_iravg_clipped, $
	sqrt(vel_iravg_clipped_err^2 + vel_oh_clipped_err^2), vel_iravg_clipped_err, $
	xtitle='v!IOH!N - v!IIR!N [km/s]', $
	ytitle='v!IIR!N [km/s]', $
	xr = [-1d3,1d3], $
	yr = [0, 9d4], $
	xthick=5, $
	ythick=5, $
	thick=5, $
	charthick=5, $
	errthick=3, $
	charsize=1.4, $
	psym = symcat(16)

ver,mean_ohir_veldiff,linestyle=2,thick=4
;ver,mean_ohir_veldiff - sig_ohir_veldiff,linestyle=2, thick=3
;ver,mean_ohir_veldiff + sig_ohir_veldiff,linestyle=2, thick=3
device,/close

device,filename='~/Astronomy/Research/Spitzer/plots/zspitzer_oh_opt.ps', /landscape,/color
ploterror, vel_oh_clipped - vel_opt_clipped, vel_oh_clipped, $
	sqrt(vel_oh_clipped_err^2 + vel_opt_clipped_err^2), vel_oh_clipped_err, $
	xtitle='v!IOH!N - v!Iopt!N [km/s]', $
	ytitle='v!IOH!N [km/s]', $
	charsize=cs, $
	psym = symcat(16)

ver,mean_ohopt_veldiff,linestyle=0,thick=4
;ver,mean_ohopt_veldiff - sig_ohopt_veldiff,linestyle=2, thick=3
;ver,mean_ohopt_veldiff + sig_ohopt_veldiff,linestyle=2, thick=3
device,/close

set_plot,'x'

; Print average value and uncertainties in various velocity differences

print,''
print,'v_opt - v_IR: ',string(mean_optir_veldiff,format='(f7.1)'),' +- ',string(sig_optir_veldiff,format='(f7.1)'),' km/s'
print,''
print,'v_OH - v_IR:  ',string(mean_ohir_veldiff,format='(f7.1)'),' +- ',string(sig_ohir_veldiff,format='(f7.1)'),' km/s'
print,''
print,'v_OH - v_opt: ',string(mean_ohopt_veldiff,format='(f7.1)'),' +- ',string(sig_ohopt_veldiff,format='(f7.1)'),' km/s'
print,''

!p.multi=[0,3,1]

; Scatter plots of velocity for each system - no systematic deviations shown for
; 	objects at higher velocities

plot, vel_iravg_clipped, vel_opt_clipped, $
	psym = symcat(15), $
	xtitle='IR avg. velocity', $
	ytitle='Optical velocity', $
	charsize=2

oplot,fillarr(1,0,1d5), fillarr(1,0,1d5)

plot, vel_iravg_clipped, vel_oh_clipped, $
	psym = symcat(15), $
	xtitle='IR avg. velocity', $
	ytitle='OHM velocity', $
	charsize=2

oplot,fillarr(1,0,1d5), fillarr(1,0,1d5)

plot, vel_oh_clipped, vel_opt_clipped, $
	psym = symcat(15), $
	xtitle='OHM velocity', $
	ytitle='Optical velocity', $
	charsize=2

oplot,fillarr(1,0,1d5), fillarr(1,0,1d5)

; Histograms of the difference between IR, OH, optical velocity measurements

plothist,vel_iravg_clipped - vel_opt_clipped,$
	bin=100, $
	xrange=[-600,600], /xstyle, $
	yrange = [0,20], /ystyle, $
	charsize = 3, $
	ytitle = 'Frequency', $
	xtitle = 'v!IIR!N - v!Iopt!N', $
	title='IR vs. optical'

ver, mean_optir_veldiff, linestyle=2

plothist,vel_iravg_clipped - vel_oh_clipped,$
	bin=1d2, $
	xrange=[-600,600], /xstyle, $
	yrange = [0,20], /ystyle, $
	charsize = 3, $
	ytitle = 'Frequency', $
	xtitle = 'v!IIR!N - v!IOH!N', $
	title='IR vs. OHM'

ver, mean_ohir_veldiff, linestyle=2

plothist,vel_opt_clipped - vel_oh_clipped, $
	bin=1d2, $
	xrange=[-600,600], /xstyle, $
	yrange = [0,20], /ystyle, $
	charsize = 3, $
	ytitle = 'Frequency', $
	xtitle = 'v!Iopt!N - v!IOH!N', $
	title='OHM vs. optical'

ver, mean_ohopt_veldiff, linestyle=2

; Determine if there is any significant difference between OHM and control sample IR-opt offset

	; Retrieve control sample optical and IR redshifts
	
	conlist = condat('tag')
	ncon = n_elements(conlist)
	hrcon = 'control'+string([23,24,25,26,40],format='(i03)')
	
	vel_con_ir = dblarr(ncon)  & vel_con_ir_err = dblarr(ncon)
	vel_con_opt = dblarr(ncon) & vel_con_opt_err = dblarr(ncon)
	
	for i = 0, ncon - 1 do begin
	
		listtag = strmid(conlist[i],0,4)
		savdir_con       = '~/Astronomy/Research/Spitzer/control/lines/hrsky/hires/saved/'
		savdir_con_nosky = '~/Astronomy/Research/Spitzer/control/lines/nosky/hires/saved/'
		if where(conlist[i] eq hrcon) ne -1 then savdir = savdir_con $
			else savdir = savdir_con_nosky
	
		restore,savdir+conlist[i]+'_zlines.sav'
	
		vel_con_ir[i] = wavg_vel & vel_con_ir_err[i] = wstddev_vel
		vel_con_opt[i] = cz_opt  & vel_con_opt_err[i] = cz_opt_err
	
	endfor
	
	; Remove missing data
	
	clipped_indices_con = where(vel_con_ir ne 0 and vel_con_opt ne 0)
	vel_con_iravg_clipped = vel_con_ir[clipped_indices_con]  & vel_con_iravg_clipped_err = vel_con_ir_err[clipped_indices_con]
	vel_con_opt_clipped   = vel_con_opt[clipped_indices_con] & vel_con_opt_clipped_err   = vel_con_opt_err[clipped_indices_con]
	
	; Weighted mean and uncertainty of velocity difference
	
	wtdmean, vel_con_opt_clipped - vel_con_iravg_clipped, sqrt(vel_con_iravg_clipped_err^2 + vel_con_opt_clipped_err^2), mean_optir_con, sig_optir_con
	
	; Run two-sided K-S test on IR-optical velocity difference in both sets
	
	var_ohm = vel_iravg_clipped - vel_opt_clipped 
	var_con = vel_con_iravg_clipped - vel_con_opt_clipped
	
	kstwo, var_ohm, var_con, D_nir, prob_nir
	gauss_nir = sqrt(2d) * inverf(1d - prob_nir)

	; Print results to screen
	
	print,''
	print,'D_KS    for IR-opt vel: '+string(D_nir,format='(f7.3)')
	print,'KS-prob for IR-opt vel: '+string(prob_nir,format='(f7.3)')
	print,'Gaussian probability:   '+string(gauss_nir,format='(f7.3)')+' sig'
	print,'Average value (OHM) :   '+string(mean_optir_veldiff,format='(f7.2)')+' +- '+string(sig_optir_veldiff,format='(f7.2)')
	print,'Average value (con) :   '+string(mean_optir_con,format='(f7.2)')+' +- '+string(sig_optir_con,format='(f7.2)')
	print,''
	
; Darling (2006) claim that blueshifted OHMs show a moderate correlation with OHM strength and width - test for OH-IR velocity offset

; Retrieve log L_OH and W_1667 data

w1667 = ohmpar('w1667',/v)
logloh = ohmpar('logl_oh',/v)

taglist_clipped = objlist[clipped_indices]
objlist_clipped = strarr(n_elements(taglist_clipped))
for i = 0, n_elements(taglist_clipped) - 1 do begin
	targets, taglist_clipped[i], r, o, z
	objlist_clipped[i] = strtrim(o,2)
endfor

match, objlist_clipped, transpose(strtrim(w1667[0,*],2)),   obj_w1667_ind,  w1667_obj_ind
match, objlist_clipped, transpose(strtrim(logloh[0,*],2)), obj_logloh_ind, logloh_obj_ind

!p.multi=[0,2,1]

plot, vel_oh_clipped[obj_w1667_ind] - vel_iravg_clipped[obj_w1667_ind], w1667[1,w1667_obj_ind], $
	charsize=2, $
	psym=symcat(16), $
	xtitle='v!IOH!N - v!IIR!N [km/s]', $
	ytitle='W!I1667!N [MHz]'

plot, vel_oh_clipped[obj_logloh_ind] - vel_iravg_clipped[obj_logloh_ind], logloh[1,logloh_obj_ind], $
	charsize=2, $
	yrange=[2,4], /ystyle, $
	psym=symcat(16), $
	xtitle='v!IOH!N - v!IIR!N [km/s]', $
	ytitle='log L!IOH!N [L!Isun!N]'

w1667_corr = r_correlate(vel_oh_clipped[obj_w1667_ind] - vel_iravg_clipped[obj_w1667_ind], w1667[1,w1667_obj_ind],d=wd,probd=wprobd,zd=wzd)
logloh_corr = r_correlate(vel_oh_clipped[obj_logloh_ind] - vel_iravg_clipped[obj_logloh_ind], logloh[1,logloh_obj_ind],d=ld,probd=lprobd,zd=lzd)

;print, 'W1667: ', w1667_corr,wd,wprobd,wzd
;print, 'log L_OH: ', logloh_corr,ld,lprobd,lzd
;print,''

w1667_pcorr = correlate(float(vel_oh_clipped[obj_w1667_ind] - vel_iravg_clipped[obj_w1667_ind]), float(w1667[1,w1667_obj_ind]))
logloh_pcorr = correlate(float(vel_oh_clipped[obj_logloh_ind] - vel_iravg_clipped[obj_logloh_ind]), float(logloh[1,logloh_obj_ind]))

;print, 'W1667: ', w1667_pcorr
;print, 'log L_OH: ', logloh_pcorr
;print,''

stop

if keyword_set(stop) then stop
end
