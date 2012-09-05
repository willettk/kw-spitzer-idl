
;+
; NAME:
;       
;	SIROCKY_BATCH
;
; PURPOSE:
;
;	Run SIROCKY.pro on all OHMs and non-masing control targets
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
;	Requires SIROCKY.pro
;
; REVISION HISTORY
;       Written by K. Willett                Apr 09
;-

allobjects = [transpose(ohmdat('tag')),transpose(archdat('tag')),transpose(condat('tag'))]
allobjects = allobjects[where(allobjects ne 'mega034' and allobjects ne 'control033')]

continuum_ohm = ['mega'+string([7,32],format='(i03)'),'arch'+string([20,23,25],format='(i03)')]
continuum_con = ['control'+string([37,35,8,36],format='(i03)')]
continuum_both = [continuum_ohm, continuum_con]

absorption_ohm = ['arch'+string([10,13,14,17,26,48,30],format='(i03)')]

match, allobjects, [continuum_both,absorption_ohm], a, b
pah_both = allobjects[setdifference(indgen(n_elements(allobjects)-1),a)]

; PAH-dominated spectra

npah = n_elements(pah_both)
tau10_pah = fltarr(npah)
tau18_pah = fltarr(npah)

for i = 0, npah - 1 do begin

	sirocky, spl=3, pah_both[i], tau10, tau10err, tau18, tau18err, /quiet

	tau10_pah[i] = tau10
	tau18_pah[i] = tau18

	if tau18 gt 0.5 then print, pah_both[i], tau18

endfor

nabs = n_elements(absorption_ohm)
tau10_abs = fltarr(nabs)
tau18_abs = fltarr(nabs)

; Absorption-dominated spectra

for i = 0, nabs - 1 do begin

	sirocky, spl=2, absorption_ohm[i], tau10, tau10err, tau18, tau18err, /quiet

	tau10_abs[i] = tau10
	tau18_abs[i] = tau18

	if tau18 gt 0.5 then print, absorption_ohm[i], tau18
endfor

ncon = n_elements(continuum_both)
tau10_con = fltarr(ncon)
tau18_con = fltarr(ncon)

; Continuum-dominated spectra

for i = 0, ncon - 1 do begin

	sirocky, spl=1, continuum_both[i], tau10, tau10err, tau18, tau18err, /quiet

	tau10_con[i] = tau10
	tau18_con[i] = tau18

	if tau18 gt 0.5 then print, continuum_both[i], tau18
endfor

tau10_ohm_pah = tau10_pah[where(strmid(pah_both,0,3) eq 'meg' or strmid(pah_both,0,3) eq 'arc')]
tau10_control_pah = tau10_pah[where(strmid(pah_both,0,3) eq 'con')]

tau10_ohm_abs = tau10_abs

tau10_ohm_con = tau10_con[where(strmid(continuum_both,0,3) eq 'meg' or strmid(continuum_both,0,3) eq 'arc')]
tau10_control_con = tau10_con[where(strmid(continuum_both,0,3) eq 'con')]

tau18_ohm_pah = tau18_pah[where(strmid(pah_both,0,3) eq 'meg' or strmid(pah_both,0,3) eq 'arc')]
tau18_control_pah = tau18_pah[where(strmid(pah_both,0,3) eq 'con')]

tau18_ohm_abs = tau18_abs

tau18_ohm_con = tau18_con[where(strmid(continuum_both,0,3) eq 'meg' or strmid(continuum_both,0,3) eq 'arc')]
tau18_control_con = tau18_con[where(strmid(continuum_both,0,3) eq 'con')]

ohm10 = [tau10_ohm_pah, tau10_ohm_abs, tau10_ohm_con]
control10 = [tau10_control_pah, tau10_control_con]

ohm18 = [tau18_ohm_pah, tau18_ohm_abs, tau18_ohm_con]
control18 = [tau18_control_pah, tau18_control_con]

ps = 1

if ps eq 1 then begin
	set_plot,'ps'
	device,filename = '~/Astronomy/Research/Spitzer/plots/sirocky.ps', /color,/landscape
	defcolor=fsc_color("Black")
endif else defcolor=fsc_color("White")

plot, ohm10, ohm18, $
	/nodata, $
	xtitle='S!I10!N', $
	ytitle='S!I18!N', $
	xrange=[-6,2], /xstyle, $
	yrange=[-1.5,1.0], /ystyle, $
	thick=4, xthick=4, ythick=4, $
	charsize = 1.5, charthick=4

oplot, ohm10, ohm18, psym=symcat(14), color=fsc_color("Red")
oplot, control10, control18, psym=symcat(15), color=fsc_color("Blue")

; Optically thin point

oplot,[1.25],[0.65], psym=symcat(16), color=defcolor, symsize=2

legend, ['OHMs','non-masing'], color=[fsc_color("Red"),fsc_color("Blue")], psym=[14,15], /top, /left, $
	thick=4, charsize=1, charthick=3

if ps eq 1 then begin
	device,/close
	set_plot,'x'
endif

end
