pro hac_avg, ps = ps
;+
; NAME:
;       
;	
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
;       Written by K. Willett                
;-

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Dark Grey")
white = fsc_color("White")
black = fsc_color("Black")

;######################
;#######  SH  #########
;######################

onames = ohmdat('tag')
no = n_elements(onames)

cnames = condat('tag')
nc = n_elements(cnames)

anames = archdat('tag')
na = n_elements(anames)

grid = fillarr(0.01,8,18)		; Avg. separation in wavelength bins is 0.0105 um

allspec = fltarr(n_elements(grid),no+na)
conspec = fltarr(n_elements(grid),nc)

; OHMs

allohms_sh = [transpose(onames),transpose(anames)]

for i = 0, no + na - 1 do begin

	tag, allohms_sh(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms_sh(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err =  sed.err_sh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15

	allspec(*,i) = newy
endfor

meanohm_sh = median(allspec,dim=2)

for i = 0, nc - 1 do begin

	tag, cnames(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames(i)+'.sav'

	wave = sed.wave_sh
	flux = sed.flux_sh
	err =  sed.err_sh

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,15.0))

	newy = newy * 0.01 / newy15

	allspec(*,i) = newy
endfor

meancon_sh = median(allspec,dim=2)

; Plot SH mean spectra without individual objects

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/hac_avg.ps', /color, /portrait 
	cs = 1.5
	lthick = 5
	medthick = 4 
	defcolor = black
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
	defcolor = white
endelse

plot,indgen(50),indgen(50), /nodata, $
	;xr = [7,19], /xstyle, $
	;yr = [-0.005,0.03], /ystyle, $
	xr = [13.0,15.5], /xstyle, $
	yr = [7.5,11.5], /ystyle, $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Flux density [mJy]', $
;	title = 'Mean SH IRS spectra of OHMs vs. non-masing ULIRGs', $
	thick = lthick, $
	xthick = lthick, $
	ythick = lthick, $
	charthick = lthick, $
	charsize = cs

oplot,newx,meanohm_sh * 1d3,color=defcolor,thick=lthick, psym=10
oplot,newx,(meancon_sh * 1d3) + 1,color=grey,thick=medthick, psym=10

;legend, /top, /left, ['OHMs ('+strtrim(n_elements(allohms_sh),2)+')', 'Control sample ('+strtrim(nc,2)+')'], $
;	color=[red, blue], $
;	linestyle=[0,0], $
;	thick = [2,2], $
;	charsize = cs


;ver, 15.0, linestyle=2

ver, 13.7, linestyle=1, thick = lthick
ver, 14.02, linestyle=1, thick = lthick
ver, 15.0, linestyle=1, thick = lthick

xyouts, /data, 13.4, 11, textoidl('C_2H_2'), charsize = cs, charthick = lthick
xyouts, /data, 14.1, 11, textoidl('HCN'), charsize = cs, charthick = lthick
xyouts, /data, 15.1, 9, textoidl('CO_2'), charsize = cs, charthick = lthick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

end
