pro spindex_sil, ps = ps, label = label, cso = cso
;+
; NAME: 
;       SPINDEX_SIL 
;
;	Written by KWW, Mar 08
;-


ohmspindex = ohmdat('spindex')
archspindex = archdat('spindex')
conspindex = condat('spindex',sz=2)

ohmsil = ohmdat('sil')
archsil = archdat('sil')
consil = condat('sil',sz=2)

!p.multi=[0,2,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Desktop/jila/spindex_sil_jila.ps', /color
	cs = 1
	ls = 2
	defcolor = fsc_color("Black")
endif else begin
	cs = 2
	ls = 1
	defcolor = fsc_color("White")
endelse

red = fsc_color("Red")
blue = fsc_color("Blue")
yellow = fsc_color("Yellow")

plot,ohmsil(0,*),ohmspindex(0,*),color=defcolor,/nodata, $
	xr=[1,-4], /xstyle, $
	yr=[-1,7], /ystyle, $
	title='NIR slope vs. silicate', $
	xtitle='9.7 !7l!3m silicate strength ', $
	ytitle = 'NIR spectral index', $
	charsize = cs, $
	thick = ls, $
	charthick = ls
oplot, ohmsil(0,*), ohmspindex(0,*), color=red, psym=symcat(15)
oplot, archsil(0,*), archspindex(0,*), color=red, psym=symcat(15)
oplot, consil(*,0), conspindex(*,0), color=blue, psym=symcat(16)

plot,ohmsil(0,*),ohmspindex(1,*),color=defcolor,/nodata, $
	xr=[1,-4], /xstyle, $
	yr=[-1,7], /ystyle, $
	title='MIR slope vs. silicate', $
	xtitle='9.7 !7l!3m silicate strength ', $
	ytitle = 'MIR spectral index', $
	charsize = cs, $
	thick = ls, $
	charthick = ls
oplot, ohmsil(0,*), ohmspindex(1,*), color=red, psym=symcat(15)
oplot, archsil(0,*), archspindex(1,*), color=red, psym=symcat(15)
oplot, consil(*,0), conspindex(*,1), color=blue, psym=symcat(16)


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

stop
end
