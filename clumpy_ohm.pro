
;+
; NAME:
;       
;	CLUMPY_OHM
;
; PURPOSE:
;
;	Analyze results for OHMs from fitting CLUMPY models
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

red = fsc_color("Red")
blue = fsc_color("Blue")

; CLUMPY - torus model
;
;; OHMs 
;
;ohmtags = [transpose(archdat('tag')), transpose(ohmdat('tag'))]
;n = n_elements(ohmtags)
;
;ohm_sigarr = fltarr(n)
;ohm_Yarr = fltarr(n)
;ohm_N0arr = fltarr(n)
;ohm_qarr = fltarr(n)
;ohm_tauvarr = fltarr(n)
;ohm_incarr = fltarr(n)
;
;for i=0, n-1 do begin
;
;	tag, ohmtags[i], dirtag
;
;	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/torus_sav/'+ohmtags[i]+'_torus_grid.sav'
;
;	minind = where(errorgrid eq min(errorgrid))
;	minvalues = errorgrid_parse(minind)
;	
;	ohm_sigarr[i] = minvalues[0]
;	ohm_Yarr[i] = minvalues[1]
;	ohm_N0arr[i] = minvalues[2]
;	ohm_qarr[i] = minvalues[3]
;	ohm_tauvarr[i] = minvalues[4]
;	ohm_incarr[i] = minvalues[5]
;
;endfor
;
;; Non-masing galaxies
;
;contags = transpose(condat('tag'))
;ncon = n_elements(contags)
;
;con_sigarr = fltarr(ncon)
;con_Yarr = fltarr(ncon)
;con_N0arr = fltarr(ncon)
;con_qarr = fltarr(ncon)
;con_tauvarr = fltarr(ncon)
;con_incarr = fltarr(ncon)
;
;for i=0, ncon-1 do begin
;
;	tag, contags[i], dirtag
;
;	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/torus_sav/'+contags[i]+'_torus_grid.sav'
;
;	minind = where(errorgrid eq min(errorgrid))
;	minvalues = errorgrid_parse(minind)
;	
;	con_sigarr[i] = minvalues[0]
;	con_Yarr[i] = minvalues[1]
;	con_N0arr[i] = minvalues[2]
;	con_qarr[i] = minvalues[3]
;	con_tauvarr[i] = minvalues[4]
;	con_incarr[i] = minvalues[5]
;
;endfor
;
;!p.multi=[0,2,3]
;
;cs = 3.0
;
;histoplot, ohm_sigarr,  charsize = cs, xtitle='!7r!3', xr=[-5,95], bin=10.
;histoplot, con_sigarr,  charsize = cs, bin=10., color="Blue", /oplot
;ver, mean(ohm_sigarr),  color=fsc_color("Red"), linestyle = 2 & ver, mean(con_sigarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_Yarr,    charsize = cs, xtitle='Y', bin=10., xr=[0, 210]
;histoplot, con_Yarr,    charsize = cs, bin=10., color="Blue", /oplot
;ver, mean(ohm_Yarr),    color=fsc_color("Red"), linestyle = 2 & ver, mean(con_Yarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_N0arr,   charsize = cs, xtitle='N!I0!N', xr=[0,25], bin=1.
;histoplot, con_N0arr,   charsize = cs, bin=1., color="Blue", /oplot
;ver, mean(ohm_N0arr),   color=fsc_color("Red"), linestyle = 2 & ver, mean(con_N0arr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_qarr,    charsize = cs, xtitle='q', xr=[-0.5,3], /xstyle, bin=0.5, /half
;histoplot, con_qarr,    charsize = cs, bin=0.5, /half, color="Blue", /oplot
;ver, mean(ohm_qarr),    color=fsc_color("Red"), linestyle = 2 & ver, mean(con_qarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_tauvarr, charsize = cs, xtitle='!7s!3!IV!N', bin=10., xr=[0, 220]
;histoplot, con_tauvarr, charsize = cs, bin=10., color="Blue", /oplot
;ver, mean(ohm_tauvarr), color=fsc_color("Red"), linestyle = 2 & ver, mean(con_tauvarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_incarr,  charsize = cs, xtitle='i', xr=[-5, 95], bin=10., /half
;histoplot, con_incarr,  charsize = cs, bin=10., /half, color="Blue", /oplot
;ver, mean(ohm_incarr),  color=fsc_color("Red"), linestyle = 2 & ver, mean(con_incarr), color=fsc_color("Blue"), linestyle = 2
;
;!p.multi=[0,1,1]
;
;; Print results to screen
;
;print,''
;print,'OHMs'
;print,''
;
;print, mean(ohm_sigarr), ' +- ', stddev(ohm_sigarr)
;print, mean(ohm_Yarr), ' +- ', stddev(ohm_Yarr)
;print, mean(ohm_N0arr), ' +- ', stddev(ohm_N0arr)
;print, mean(ohm_qarr), ' +- ', stddev(ohm_qarr)
;print, mean(ohm_tauvarr), ' +- ', stddev(ohm_tauvarr)
;print, mean(ohm_incarr), ' +- ', stddev(ohm_incarr)
;
;
;print,''
;print,'Non-masing galaxies'
;print,''
;
;print, mean(con_sigarr), ' +- ', stddev(con_sigarr)
;print, mean(con_Yarr), ' +- ', stddev(con_Yarr)
;print, mean(con_N0arr), ' +- ', stddev(con_N0arr)
;print, mean(con_qarr), ' +- ', stddev(con_qarr)
;print, mean(con_tauvarr), ' +- ', stddev(con_tauvarr)
;print, mean(con_incarr), ' +- ', stddev(con_incarr)

;; CLUMPY - spherical model
;
;; OHMs 
;
;ohmtags = [transpose(archdat('tag')), transpose(ohmdat('tag'))]
;n = n_elements(ohmtags)
;
;ohm_Yarr = fltarr(n)
;ohm_N0arr = fltarr(n)
;ohm_qarr = fltarr(n)
;ohm_tauvarr = fltarr(n)
;
;for i=0, n-1 do begin
;
;	tag, ohmtags[i], dirtag
;
;	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/sphere_sav/'+ohmtags[i]+'_sphere_grid.sav'
;
;	minind = where(errorgrid_sphere eq min(errorgrid_sphere))
;	minvalues = errorgrid_parse(minind,/sphere)
;	
;	ohm_Yarr[i] = minvalues[0]
;	ohm_N0arr[i] = minvalues[1]
;	ohm_qarr[i] = minvalues[2]
;	ohm_tauvarr[i] = minvalues[3]
;
;endfor
;
;; Non-masing galaxies
;
;contags = transpose(condat('tag'))
;ncon = n_elements(contags)
;
;con_Yarr = fltarr(ncon)
;con_N0arr = fltarr(ncon)
;con_qarr = fltarr(ncon)
;con_tauvarr = fltarr(ncon)
;
;for i=0, ncon-1 do begin
;
;	tag, contags[i], dirtag
;
;	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/CLUMPY/sphere_sav/'+contags[i]+'_sphere_grid.sav'
;
;	minind = where(errorgrid_sphere eq min(errorgrid_sphere))
;	minvalues = errorgrid_parse(minind,/sphere)
;	
;	con_Yarr[i] = minvalues[0]
;	con_N0arr[i] = minvalues[1]
;	con_qarr[i] = minvalues[2]
;	con_tauvarr[i] = minvalues[3]
;
;endfor
;
;!p.multi=[0,2,2]
;
;cs = 1.5
;
;histoplot, ohm_Yarr,    charsize = cs, xtitle='Y', bin=10., xr=[0, 40]
;histoplot, con_Yarr,    charsize = cs, bin=10., color="Blue", /oplot, /half
;ver, mean(ohm_Yarr),    color=fsc_color("Red"), linestyle = 2 & ver, mean(con_Yarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_N0arr,   charsize = cs, xtitle='N!I0!N', xr=[0,25], bin=2.
;histoplot, con_N0arr,   charsize = cs, bin=2., color="Blue", /oplot
;ver, mean(ohm_N0arr),   color=fsc_color("Red"), linestyle = 2 & ver, mean(con_N0arr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_qarr,    charsize = cs, xtitle='q', xr=[-0.5,3], /xstyle, bin=0.5, /half
;histoplot, con_qarr,    charsize = cs, bin=0.5, /half, color="Blue", /oplot
;ver, mean(ohm_qarr),    color=fsc_color("Red"), linestyle = 2 & ver, mean(con_qarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_tauvarr, charsize = cs, xtitle='!7s!3!IV!N', bin=20., xr=[0, 220]
;histoplot, con_tauvarr, charsize = cs, bin=20., color="Blue", /oplot
;ver, mean(ohm_tauvarr), color=fsc_color("Red"), linestyle = 2 & ver, mean(con_tauvarr), color=fsc_color("Blue"), linestyle = 2
;
;!p.multi=[0,1,1]
;
;; Print results to screen
;
;print,''
;print,'CLUMPY, spherical geometry'
;print,''
;print,'OHMs'
;print,''
;
;print, mean(ohm_Yarr), ' +- ', stddev(ohm_Yarr)
;print, mean(ohm_N0arr), ' +- ', stddev(ohm_N0arr)
;print, mean(ohm_qarr), ' +- ', stddev(ohm_qarr)
;print, mean(ohm_tauvarr), ' +- ', stddev(ohm_tauvarr)
;
;
;print,''
;print,'Non-masing galaxies'
;print,''
;
;print, mean(con_Yarr), ' +- ', stddev(con_Yarr)
;print, mean(con_N0arr), ' +- ', stddev(con_N0arr)
;print, mean(con_qarr), ' +- ', stddev(con_qarr)
;print, mean(con_tauvarr), ' +- ', stddev(con_tauvarr)

;; DUSTY
;
;; OHMs 
;
;ohmtags = [transpose(archdat('tag')), transpose(ohmdat('tag'))]
;n = n_elements(ohmtags)
;
;ohm_Yarr = fltarr(n)
;ohm_qarr = fltarr(n)
;ohm_tauvarr = fltarr(n)
;ohm_tdustarr = fltarr(n)
;
;for i=0, n-1 do begin
;
;	tag, ohmtags[i], dirtag
;
;	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+ohmtags[i]+'_sphere_grid.sav'
;
;	minind = where(errorgrid eq min(errorgrid))
;	minvalues = errorgrid_parse(minind[0],/dusty)
;	
;	sppfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/kw_sphere_Y'+$
;		string(minvalues[0],format='(i04)')+'_q'+string(minvalues[1],format='(f3.1)')+'.spp'
;	
;	readcol, sppfile, model, tau0, /silent
;	
;	ohm_Yarr[i] = minvalues[0]
;	ohm_qarr[i] = minvalues[1]
;	ohm_tauvarr[i] = tau0[minvalues[2]]
;	ohm_tdustarr[i] = tdust
;
;	print,minind
;
;endfor
;
;; Non-masing galaxies
;
;contags = transpose(condat('tag'))
;ncon = n_elements(contags)
;
;con_Yarr = fltarr(ncon)
;con_qarr = fltarr(ncon)
;con_tauvarr = fltarr(ncon)
;con_tdustarr = fltarr(ncon)
;
;for i=0, ncon-1 do begin
;
;	tag, contags[i], dirtag
;
;	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+contags[i]+'_sphere_grid.sav'
;
;	minind = where(errorgrid eq min(errorgrid))
;	minvalues = errorgrid_parse(minind[0],/dusty)
;	
;	sppfile='~/Astronomy/Research/Spitzer/CSO/DUSTY/kw_sphere_Y'+$
;		string(minvalues[0],format='(i04)')+'_q'+string(minvalues[1],format='(f3.1)')+'.spp'
;	
;	readcol, sppfile, model, tau0, /silent
;	
;	con_Yarr[i] = minvalues[0]
;	con_qarr[i] = minvalues[1]
;	con_tauvarr[i] = tau0[minvalues[2]]
;	con_tdustarr[i] = tdust
;
;endfor
;
;!p.multi=[0,2,2]
;
;cs = 1.5
;
;histoplot, ohm_Yarr,     charsize = cs, xtitle='Y', bin=10., xr=[0, 210]
;histoplot, con_Yarr,     charsize = cs, bin=10., color="Blue", /oplot
;ver, mean(ohm_Yarr),     color=fsc_color("Red"), linestyle = 2 & ver, mean(con_Yarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_qarr,     charsize = cs, xtitle='q', xr=[-0.5,3], /xstyle, bin=0.5
;histoplot, con_qarr,     charsize = cs, bin=0.5, color="Blue", /oplot
;ver, mean(ohm_qarr),     color=fsc_color("Red"), linestyle = 2 & ver, mean(con_qarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_tauvarr,  charsize = cs, xtitle='!7s!3!IV!N', bin=10., xr=[0, 220]
;histoplot, con_tauvarr,  charsize = cs, bin=10., color="Blue", /oplot
;ver, mean(ohm_tauvarr),  color=fsc_color("Red"), linestyle = 2 & ver, mean(con_tauvarr), color=fsc_color("Blue"), linestyle = 2
;histoplot, ohm_tdustarr, charsize = cs, xtitle='T!Idust!N', xr=[0,200], bin=10.
;histoplot, con_tdustarr, charsize = cs, bin=10., color="Blue", /oplot
;ver, mean(ohm_tdustarr), color=fsc_color("Red"), linestyle = 2 & ver, mean(con_tdustarr), color=fsc_color("Blue"), linestyle = 2
;
;!p.multi=[0,1,1]
;
;; Print results to screen
;
;print,''
;print,'DUSTY - spherical shell'
;print,''
;print,'OHMs'
;print,''
;
;print, 'Y:      ',mean(ohm_Yarr), ' +- ', stddev(ohm_Yarr)
;print, 'q:      ',mean(ohm_qarr), ' +- ', stddev(ohm_qarr)
;print, 'tau_V:  ',mean(ohm_tauvarr), ' +- ', stddev(ohm_tauvarr)
;print, 'T_dust: ',mean(ohm_tdustarr), ' +- ', stddev(ohm_tdustarr)
;
;
;print,''
;print,'Non-masing galaxies'
;print,''
;
;print, 'Y:      ',mean(con_Yarr), ' +- ', stddev(con_Yarr)
;print, 'q:      ',mean(con_qarr), ' +- ', stddev(con_qarr)
;print, 'tau_V:  ',mean(con_tauvarr), ' +- ', stddev(con_tauvarr)
;print, 'T_dust: ',mean(con_tdustarr), ' +- ', stddev(con_tdustarr)

; DUSTY - PAHFIT

; OHMs 

ohmtags = [transpose(archdat('tag')), transpose(ohmdat('tag'))]
n = n_elements(ohmtags)

ohm_Yarr = fltarr(n)
ohm_qarr = fltarr(n)
ohm_tauvarr = fltarr(n)
ohm_tdustarr = fltarr(n)

for i=0, n-1 do begin

	tag, ohmtags[i], dirtag

	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+ohmtags[i]+'_sphere_pahfitgrid.sav'

	ohm_Yarr[i] = y_min
	ohm_qarr[i] = q_min
	ohm_tauvarr[i] = tauv_min
	ohm_tdustarr[i] = tdust

endfor

; Non-masing galaxies

contags = transpose(condat('tag'))
ncon = n_elements(contags)

con_Yarr = fltarr(ncon)
con_qarr = fltarr(ncon)
con_tauvarr = fltarr(ncon)
con_tdustarr = fltarr(ncon)

for i=0, ncon-1 do begin

	tag, contags[i], dirtag

	restore, filename='~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+contags[i]+'_sphere_pahfitgrid.sav'

	con_Yarr[i] = y_min
	con_qarr[i] = q_min
	con_tauvarr[i] = tauv_min
	con_tdustarr[i] = tdust

endfor

save, file='~/Astronomy/Research/Spitzer/dusty_pahfit_results.sav', $
	ohmtags, ohm_yarr, ohm_qarr, ohm_tauvarr, ohm_tdustarr, $
	contags, con_yarr, con_qarr, con_tauvarr, con_tdustarr

!p.multi=[0,2,2]

ps = 1

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/dusty_pahfit_hist.ps', /color, /portrait, decomposed=1
	cs = 1.0
	defcolor=fsc_color("Black")
	thick=4
endif else begin
	cs=2
	defcolor=fsc_color("White")
	thick=1
endelse

plothist, ohm_Yarr, charsize = cs, xtitle='Y', bin=200., xr=[0, 1100], /half, $
	thick=thick, xthick=thick, ythick=thick, charthick=thick, color=defcolor, /nodata
plothist, ohm_Yarr, bin=200., /half, /overplot, thick=thick, color=red
plothist, con_Yarr, bin=200., /half, /overplot, thick=thick, color=blue
ver, mean(ohm_Yarr), color=red, linestyle = 2, thick=thick 
ver, mean(con_Yarr), color=blue, linestyle = 2, thick=thick

plothist, ohm_qarr, charsize = cs, xtitle='q', bin=0.5, xr=[0,3.0], /xstyle, /half, $
	thick=thick, xthick=thick, ythick=thick, charthick=thick, color=defcolor, /nodata
plothist, ohm_qarr, bin=0.5, /half, /overplot, thick=thick, color=red
plothist, con_qarr, bin=0.5, /half, /overplot, thick=thick, color=blue
ver, mean(ohm_qarr), color=red, linestyle = 2, thick=thick 
ver, mean(con_qarr), color=blue, linestyle = 2, thick=thick

plothist, ohm_tauvarr, charsize = cs, xtitle='!7s!3!IV!N', bin=50., xr=[0,450], yr=[0,35], /ystyle, /half, $
	thick=thick, xthick=thick, ythick=thick, charthick=thick, color=defcolor, /nodata
plothist, ohm_tauvarr, bin=50., /half, /overplot, thick=thick, color=red
plothist, con_tauvarr, bin=50., /half, /overplot, thick=thick, color=blue
ver, mean(ohm_tauvarr), color=red, linestyle = 2, thick=thick 
ver, mean(con_tauvarr), color=blue, linestyle = 2, thick=thick

plothist, ohm_tdustarr, charsize = cs, xtitle='Dust temperature [K]', bin=15., xr=[0,200], /half, $
	thick=thick, xthick=thick, ythick=thick, charthick=thick, color=defcolor, /nodata
plothist, ohm_tdustarr, bin=15., /half, /overplot, thick=thick, color=red
plothist, con_tdustarr, bin=15., /half, /overplot, thick=thick, color=blue
ver, mean(ohm_tdustarr), color=red, linestyle = 2, thick=thick 
ver, mean(con_tdustarr), color=blue, linestyle = 2, thick=thick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/papers/tdust_logoh.ps', /color, /portrait
	cs = 1.5
	defcolor=fsc_color("Black")
	thick=5
	cthick=4
	hsize=250
endif else begin
	cs=2
	defcolor=fsc_color("White")
	thick=1
	cthick=1
	hsize=30
endelse

logoh = float([transpose(archdat('logoh')), transpose(ohmdat('logoh'))])
logoh_lim = float([transpose(condat('logoh'))])
plot, ohm_tdustarr, logoh, /nodata, $
	xr=[30,115], $
	yr = [0,4], $
	charsize = cs, $
	charthick = cthick, $
	thick=thick, $
	xthick=thick, $
	ythick=thick, $
	xtitle='Dust temperature [K]', $
	ytitle='log (L!IOH!N [L!I'+sunsymbol()+'])'
oplot, ohm_tdustarr, logoh, psym=symcat(14), color=defcolor, symsize=1.2
oplot, con_tdustarr, logoh_lim, psym=symcat(7), color=defcolor, thick=thick, symsize=1.2
arrow, con_tdustarr, logoh_lim, con_tdustarr, logoh_lim - 0.4, color=defcolor, /data, thick=thick, hsize=hsize
print,''
print,'Tdust - log (LOH) correlation: ',correlate(ohm_tdustarr,logoh)

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif


kstwo, ohm_tdustarr,con_tdustarr, D_tdust, prob_tdust
gauss_tdust = sqrt(2d) * inverf(1d - prob_tdust)

print,''
print,'D_KS    for Tdust: '+string(D_tdust,format='(f7.3)')
print,'KS-prob for Tdust: '+string(prob_tdust,format='(e8.2)')
print,'Gaussian probability:  '+string(gauss_tdust,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_tdustarr),format='(f7.2)')+' +- '+string(stddev(ohm_tdustarr),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_tdustarr),format='(f7.2)')+' +- '+string(stddev(con_tdustarr),format='(f7.2)')

kstwo, ohm_tauvarr,con_tauvarr, D_tauv, prob_tauv
gauss_tauv = sqrt(2d) * inverf(1d - prob_tauv)

print,''
print,'D_KS    for tau_V: '+string(D_tauv,format='(f7.3)')
print,'KS-prob for tau_V: '+string(prob_tauv,format='(e8.2)')
print,'Gaussian probability:  '+string(gauss_tauv,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_tauvarr),format='(f7.2)')+' +- '+string(stddev(ohm_tauvarr),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_tauvarr),format='(f7.2)')+' +- '+string(stddev(con_tauvarr),format='(f7.2)')

kstwo, ohm_yarr,con_yarr, D_y, prob_y
gauss_y = sqrt(2d) * inverf(1d - prob_y)

print,''
print,'D_KS    for Y: '+string(D_y,format='(f7.3)')
print,'KS-prob for Y: '+string(prob_y,format='(e8.2)')
print,'Gaussian probability:  '+string(gauss_y,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_yarr),format='(f7.2)')+' +- '+string(stddev(ohm_yarr),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_yarr),format='(f7.2)')+' +- '+string(stddev(con_yarr),format='(f7.2)')

kstwo, ohm_qarr,con_qarr, D_q, prob_q
gauss_q = sqrt(2d) * inverf(1d - prob_q)

print,''
print,'D_KS    for q: '+string(D_q,format='(f7.3)')
print,'KS-prob for q: '+string(prob_q,format='(e8.2)')
print,'Gaussian probability:  '+string(gauss_q,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_qarr),format='(f7.2)')+' +- '+string(stddev(ohm_qarr),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_qarr),format='(f7.2)')+' +- '+string(stddev(con_qarr),format='(f7.2)')

; Print results to screen

print,''
print,'DUSTY (PAHFIT) - spherical shell'
print,''
print,'OHMs'
print,''

print, 'Y:      ',mean(ohm_Yarr), ' +- ', stddev(ohm_Yarr)
print, 'q:      ',mean(ohm_qarr), ' +- ', stddev(ohm_qarr)
print, 'tau_V:  ',mean(ohm_tauvarr), ' +- ', stddev(ohm_tauvarr)
print, 'T_dust: ',mean(ohm_tdustarr), ' +- ', stddev(ohm_tdustarr)


print,''
print,'Non-masing galaxies'
print,''

print, 'Y:      ',mean(con_Yarr), ' +- ', stddev(con_Yarr)
print, 'q:      ',mean(con_qarr), ' +- ', stddev(con_qarr)
print, 'tau_V:  ',mean(con_tauvarr), ' +- ', stddev(con_tauvarr)
print, 'T_dust: ',mean(con_tdustarr), ' +- ', stddev(con_tdustarr)

;stop

end
