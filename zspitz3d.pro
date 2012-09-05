;pro zspitz3d
;+
; NAME:
;       
;	ZSPITZ3D
;
; PURPOSE:
;
;	Plot results for each line in 3D space
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
;       Written by K. Willett                Apr 08
;-

ohmlist = ohmdat('tag')
nolines = where(ohmlist eq 'mega016')
noind = setdifference(indgen(n_elements(ohmlist)),[nolines])
ohmlist = ohmlist(noind)
nohm = n_elements(ohmlist)

c = 299792.458d
oharr = double(ohmdat('czoh'))
veloharr = oharr(noind)
zoharr = oharr(noind) / c

savdir = '~/Astronomy/Research/Spitzer/ohm/lines/hrsky/hires/saved/'

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

!p.multi = [0,2,2]
wset,0

wavgarr = dblarr(nohm)
roptarr = dblarr(nohm)

ps = 0

for k = 0,2 do begin

	if keyword_set(ps) then begin
		set_plot,'ps'
		device,filename='~/Astronomy/Research/Spitzer/plots/zspitzer'+strtrim(k,2)+'.ps', /landscape, /color
		defcolor = fsc_color("Black")
		cs = 1
		ss = 0.5
	endif else begin
		defcolor = fsc_color("White")
		cs = 3
		ss = 1
	endelse
	
	if k eq 2 then j_end = 3 else j_end = 3
	for j = 0, j_end do begin
	
		ionzarr = dblarr(nohm)
		dvelarr_wavg = dblarr(nohm)
		dvelarr_ropt = dblarr(nohm)
		dvelarr_czoh = dblarr(nohm)
		delvarr = dblarr(nohm)
		
		for i = 0, nohm-1 do begin
			restore,savdir+ohmlist(i)+'_hrlines.sav'
			ion_now = ionlist(4*k + j)
			czoh = zoharr(i)
			isionthere = where(ion eq ion_now)
			if isionthere ge 0 then begin
				ionzarr(i) = zarr(isionthere)
				dvelarr_wavg(i) = veldiff(wavg,zarr(isionthere))
				dvelarr_ropt(i) = veldiff(ropt,zarr(isionthere))
				dvelarr_czoh(i) = veldiff(czoh,zarr(isionthere))
			endif
			wavgarr(i) = wavg
			roptarr(i) = ropt
			delvarr(i) = delv
		endfor
		
		remove_zeroes = where(ionzarr ne 0)
		ionzarr = ionzarr(remove_zeroes)
		dvelarr_wavg = dvelarr_wavg(remove_zeroes)
		dvelarr_ropt = dvelarr_ropt(remove_zeroes)
		dvelarr_czoh = dvelarr_czoh(remove_zeroes)
	
		print,''
		print,'Average vel. difference of '+ionlist(4*k + j)+' from wavg: ',$
			mean(dvelarr_wavg),' +- ',stddev(dvelarr_wavg), ' km/s'
		print,'Average vel. difference of '+ionlist(4*k + j)+' from ropt: ',$
			mean(dvelarr_ropt),' +- ',stddev(dvelarr_ropt), ' km/s'
		print,'Average vel. difference of '+ionlist(4*k + j)+' from czoh: ',$
			mean(dvelarr_czoh),' +- ',stddev(dvelarr_czoh), ' km/s'
		
		
		!p.multi=[0,1,1]
		surface, dist(5), /nodata, /save, $
;			psym = symcat(16), $
			xtitle = '!7D!3v from IR data [km/s]', $
			ytitle = '!7D!3v from optical data [km/s]', $
			ztitle = '!7D!3v from OH data [km/s]', $
			title = ionlist(4*k + j), $
			charsize = cs, $
			xrange = [-500,800], /xstyle, $
			yrange = [-500,800], /ystyle, $
			zrange = [-500,800], /zstyle, $
			xgridstyle = 1, $
			ygridstyle = 1, $
			xticklen = 1, $
			yticklen = 1, $
			ax = 40

		axis, xaxis = 1, /t3d, charsize = cs, xr=[-500,800], /xstyle, xticklen = 1
		axis, yaxis = 1, /t3d, charsize = cs, yr=[-500,800], /ystyle, yticklen = 1
		
;		axcolor=fsc_color("Yellow")
;		plots, [-500,800],[0,0],[0,0],linestyle=0,/t3d,color=axcolor
;		plots, [0,0],[-500,800],[0,0],linestyle=0,/t3d,color=axcolor
;		plots, [0,0],[0,0],[-500,800],linestyle=0,/t3d,color=axcolor
		
;		oplot,dvelarr_wavg, dvelarr_ropt, psym=symcat(16), symsize=ss		; IR vs. optical
;		oplot,dvelarr_wavg, dvelarr_czoh, psym=symcat(16), symsize=ss		; IR vs. OH
;		oplot,dvelarr_ropt, dvelarr_czoh, psym=symcat(16), symsize=ss		; Optical vs. OH
	
		plots,dvelarr_wavg,dvelarr_ropt,dvelarr_czoh,$
			psym=symcat(16),symsize=1.5,color=fsc_color("Red"), /t3d

		for m=0,n_elements(dvelarr_wavg) - 1 do $
			plots,[dvelarr_wavg(m),dvelarr_wavg(m)], [dvelarr_ropt(m),dvelarr_ropt(m)], [-500,dvelarr_czoh(m)], $
			color=fsc_color("Red"), /t3d
;		ver, 0, linestyle=1
;		hor, 0, linestyle=1
		
		wait,5
	endfor
endfor

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

stop
!p.multi=[0,3,1]
wset,1

plot, wavgarr, roptarr, $
	psym = symcat(15), $
	xrange=[0.1,0.2], /xstyle, $
	yrange=[0.1,0.2], /ystyle, $
	xtitle='Redshift from weighted IR average', $
	ytitle='Redshift from optical data', $
	charsize=2

oplot,fillarr(1,-500,500),fillarr(1,-500,500)
;xyouts,charsize=3,0.11,0.18,'IR emission is blueshifted',color=fsc_color("Blue")
;xyouts,charsize=3,0.14,0.12,'IR emission is redshifted',color=fsc_color("Red")

plot, wavgarr, zoharr, $
	psym = symcat(15), $
	xrange=[0.1,0.2], /xstyle, $
	yrange=[0.1,0.2], /ystyle, $
	xtitle='Redshift from weighted IR average', $
	ytitle='Redshift of OHM (1667 MHz)', $
	charsize=2
oplot,fillarr(1,-500,500),fillarr(1,-500,500)

plot, roptarr, zoharr, $
	psym = symcat(15), $
	xrange=[0.1,0.2], /xstyle, $
	yrange=[0.1,0.2], /ystyle, $
	xtitle='Redshift from optical data', $
	ytitle='Redshift of OHM (1667 MHz)', $
	charsize=2
oplot,fillarr(1,-500,500),fillarr(1,-500,500)

bsize = 5d-4
xrhist = [-3d-3,3d-3]
yrhist = [0,12]

plothist,wavgarr - roptarr,$
	bin=bsize, $
	xrange=xrhist, /xstyle, $
	yrange = yrhist, /ystyle, $
	charsize = 3, $
	ytitle = 'Frequency', $
	xtitle = 'z!IIR!N - z!Iopt!N', $
	title='IR vs. optical'

ver, 0, linestyle=2

plothist,wavgarr - zoharr,$
	bin=bsize, $
	xrange=xrhist, /xstyle, $
	yrange = yrhist, /ystyle, $
	charsize = 3, $
	ytitle = 'Frequency', $
	xtitle = 'z!IIR!N - z!IOH!N', $
	title='IR vs. OHM'

ver, 0, linestyle=2

plothist,zoharr - roptarr,$
	bin=bsize, $
	xrange=xrhist, /xstyle, $
	yrange = yrhist, /ystyle, $
	charsize = 3, $
	ytitle = 'Frequency', $
	xtitle = 'z!IOH!N - z!Iopt!N', $
	title='OHM vs. optical'

ver, 0, linestyle=2

stop
end
