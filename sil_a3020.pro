pro sil_a3020, ps = ps, stop = stop

;+
; NAME:
;       
;	SIL_A3020
;
; PURPOSE:
;
;	Plots 9.7 um silicate depth vs. the 30-20 um spectral index
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
;       Written by K. Willett                Feb 10
;-

asil = archdat('sil')
osil = ohmdat('sil')
csil = condat('sil')

ohmsil = float([transpose(asil[0,*]),transpose(osil[0,*])])
consil = float([transpose(csil[0,*])])

aspindex = archdat('spindex')
ospindex = ohmdat('spindex')
cspindex = condat('spindex')

ohmspindex = float([transpose(aspindex[1,*]),transpose(ospindex[1,*])])
conspindex = float([transpose(cspindex[1,*])])

aloh = archdat('logoh')
oloh = ohmdat('logoh')
cloh = condat('logoh')

ohmloh = float([transpose(aloh[0,*]),transpose(oloh[0,*])])
conloh = float([transpose(cloh[0,*])])

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/sil_a3020.ps'
	thick = 5
	cthick = 5
	symthick = 4
	csize = 2
	legsize = 1.3
endif else begin
	thick = 1
	cthick = 1
	csize = 2
	legsize = 1.0
endelse

plot, ohmsil, ohmspindex, $
	/nodata, $
	xr = [-4,1], $
	yr = [0,8], $
	xtitle='S!I9.7!N', $
	ytitle='!7a!3!I30-20!N', $
	thick = thick, $
	xthick = thick, $
	ythick = thick, $
	charthick = cthick, $
	charsize = csize

for i = 0, n_elements(ohmsil)-1 do $
	oplot, [ohmsil[i]], [ohmspindex[i]], psym = symcat(9,thick=symthick), symsize = ohmloh[i] - 1
oplot, consil, conspindex, psym = symcat(7), thick=thick

legend, /top, /right, $
	psym=[9,7], ['OHMs','non-masing'], $
	charsize=legsize, thick=thick, charthick=cthick, symthick=symthick

; Draw rough locus of separation between masing and non-masing populations; somewhat arbitrary
;	regarding exact location

xlocus = [-2, -0.5]
ylocus = [2.9, 4.3]

plots, xlocus, ylocus, linestyle=2, thick=cthick

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
