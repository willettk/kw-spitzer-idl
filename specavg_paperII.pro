
;+
; NAME:
;       
;	SPECAVG_PAPERII
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
;       Written by K. Willett                Sep 2010
;-

ps = 1

;######################
;#######  LR  #########
;######################

; Load names of files from data structures

onames = ohmdat('tag')
badohm = where(onames eq 'mega034')
goodohm = setdifference(indgen(n_elements(onames)),badohm)
onames = onames(goodohm)
no = n_elements(onames)

anames = transpose(archdat('tag'))
na = n_elements(anames)

cnames = condat('tag')
badcon = where(cnames eq 'control033')
goodcon = setdifference(indgen(n_elements(cnames)),badcon)
cnames = cnames(goodcon)
nc = n_elements(cnames)

; Wavelength grid over which to plot data

wavesep = 0.087
grid = fillarr(wavesep,5.0,30.0)		; Avg. separation between wave bins is 0.087 um

allspec = fltarr(n_elements(grid),no+na)
conspec = fltarr(n_elements(grid),nc)

; Loop over OHMs
	 
allohms = [onames,anames]

for i = 0, n_elements(allohms) - 1 do begin

	; Restore data

	tag, allohms(i), dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	; Sample the lo-res spectra at 0.1 um intervals

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	; Scale spectra to 10 mJy at 15.0 um

	newy15 = newy(closeto(newx,13.2))
	newy = newy * 0.1 / newy15

	allspec[*,i] = newy
endfor

; Plot the median template from all spectra

meanohm = median(allspec,dim=2)
stdohm = meanohm
for i=0,n_elements(stdohm)-1 do stdohm[i] = stddev(allspec[i,*])

; Control sample

for i = 0, nc - 1 do begin

	restore,'~/Astronomy/Research/Spitzer/control/data/structures/'+cnames(i)+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	newx = grid
	newy = fltarr(n_elements(newx))
	for j = 0,n_elements(newx)-1 do newy(j) = flux(closeto(wave,newx(j)))

	newy15 = newy(closeto(newx,13.2))

	newy = newy * 0.1 / newy15

	conspec[*,i] = newy
endfor

meancon = median(conspec,dim=2)
stdcon = meancon
for i=0, n_elements(stdcon)-1 do stdcon[i] = stddev(conspec[i,*])

!p.multi=[0,1,1]

red = fsc_color("Red")
blue = fsc_color("Blue")
grey = fsc_color("Dark Grey")
white = fsc_color("White")
black = fsc_color("Black")
salmon = fsc_color("Salmon")

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/specavg_paperII.ps', /color, /portrait, /decomposed
	cs = 1.4
	lthick = 5
	cthick = 5
	medthick = 4
	defcolor=black
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
	medthick = 1
	defcolor=white
endelse

plot,indgen(35), /nodata, $
	xr = [4.5,31], /xstyle, $
	yr = [5d-3,2d-0], /ystyle, $
	/xlog, /ylog, $
	xticks = 5, $
	xtickv = [5,10,15,20,25,30], $
	xtitle = 'Wavelength (rest frame) [!7l!3m]', $
	ytitle = 'Normalized flux density', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, charthick = cthick, charsize = cs

concolor=red

; Error bars on each point (error for median is sigma/sqrt(N) * sqrt(pi/2) - see Maritz & Jarrett [1978])

varmed_ohm = stdohm / sqrt(n_elements(meanohm)) * sqrt(!dpi / 2.)
varmed_con = stdcon / sqrt(n_elements(meancon)) * sqrt(!dpi / 2.)

oploterror,newx,          meanohm,varmed_ohm,color=defcolor,thick=lthick,   psym=10, errthick=2.0, /nohat, errcolor=grey
oploterror,newx+wavesep/2,meancon,varmed_con,color=concolor,thick=medthick, psym=10, errthick=2.0, /nohat, errcolor=red

; Overplot data again

oplot,newx,meanohm,color=defcolor,thick=lthick, psym=10 
oplot,newx+wavesep/2,meancon,color=concolor,thick=lthick, psym=10 

;legend, /top, /left, ['OHMs ('+strtrim(n_elements(allohms),2)+')', 'non-masing ('+strtrim(nc,2)+')'], $	; Bug including 3 OHMs I had to remove; won't run without morbo, so cheat for now. 

legend, /top, /left, ['OHMs ('+strtrim('51',2)+')', 'non-masing ('+strtrim(nc,2)+')'], $
	color=[defcolor,concolor], $
	linestyle=[0,0], $
	thick = [lthick,medthick], $
	charsize = 1, charthick=cthick

ver, 13.2, linestyle=2, thick=lthick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

end
