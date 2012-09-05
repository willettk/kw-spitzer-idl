; pro sil_hist
;+
; NAME:
;       
;	SIL_HIST
;
; PURPOSE:
;
;	Make a histogram of the silicate differences between OHMs and control sample
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
;       Written by K. Willett                Nov 08
; 	Added T_{30-20} - Apr 09
;-

; Load data

o = ohmdat('sil') & o = float(o[0,*])
a = archdat('sil') & a = float(a[0,*])
c = condat('sil') & c = float(c[0,*])

ohms = [transpose(o),transpose(a)]
con = transpose(c)

ohms = ohms[where(ohms ne 0.)]
con = con[where(con ne 0.)]

ps = 1

!p.multi=[0,2,1]

if ps eq 1 then begin
	set_plot,'ps'
	device, filename='~/Astronomy/Research/Spitzer/plots/sil_hist.ps',/landscape
	lthick = 3
	cthick = 3
	cs = 2
endif else begin
	lthick = 1
	cthick = 1
	cs = 1
endelse

plothist, ohms, $
	bin = 0.5, $
	xtitle = 'S!Isil!N', $
	ytitle = 'Frequency', $
	xthick = lthick, $
	ythick = lthick, $
	thick = lthick, $
	charthick = cthick, $
	charsize = cs, $
	linestyle = 2, $
	xr = [-4, 1.0]

plothist, con, $
	/overplot, $
	bin = 0.5, $
	thick = lthick, $
	linestyle = 0

if ps eq 1 then begin
	device,/close
	set_plot,'x'
endif

end
