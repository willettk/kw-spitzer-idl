
;+
; NAME:
;       
;	SFR_IR_OHM
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
;       Written by K. Willett              Mar 2010  
;-

mpc2cm = 3.086d24 
lsun   = 3.862d33 

aobj = archdat('obj')
oobj = ohmdat('obj')

ohmobj = [transpose(aobj),transpose(oobj)]
conobj = transpose(condat('obj'))

atag = archdat('tag')
otag = ohmdat('tag')

ohmtag = [transpose(atag),transpose(otag)]
contag = transpose(condat('tag'))

; IR data

airas_string = float(archdat('lfir'))
airas = transpose(10.^airas_string * lsun)		; L_FIR in erg/s

oiras_string = float(ohmdat('lfir'))
oiras = transpose(10.^oiras_string * lsun)		; L_FIR in erg/s

ciras_string = float(condat('lfir'))
ciras = transpose(10.^ciras_string * lsun)		; L_FIR in erg/s

ohm_iras = [airas,oiras]
con_iras = ciras

; PAH

apah62 = float(archdat('pah62lum'))
apah11 = float(archdat('pah11lum'))
opah62 = float(ohmdat('pah62lum'))
opah11 = float(ohmdat('pah11lum'))
cpah62 = float(condat('pah62lum'))
cpah11 = float(condat('pah11lum'))

apah = transpose(10d^(apah62) + 10d^(apah11)) * lsun	; L_PAH in erg/s
opah = transpose(10d^(opah62) + 10d^(opah11)) * lsun	; L_PAH in erg/s
cpah = transpose(10d^(cpah62) + 10d^(cpah11)) * lsun	; L_PAH in erg/s

ohm_pah = [apah, opah]
con_pah = cpah

; Neon

ohm_ne_flux = fltarr(n_elements(ohmtag))
con_ne_flux = fltarr(n_elements(contag))

for i=0, n_elements(ohmtag)-1 do begin
	neii = getlineflux(ohmtag[i], 'neII',/quiet)
	neiii = getlineflux(ohmtag[i], 'neIII',/quiet)
	neii = neii[0] & neiii = neiii[0]

	if neii gt 0 and neiii gt 0 then ohm_ne_flux[i] = neii + neiii $
		else if neii gt 0 and neiii eq -1 then ohm_ne_flux[i] = neii $
		else if neiii gt 0 and neii eq -1 then ohm_ne_flux[i] = neiii

endfor

for i=0, n_elements(contag)-1 do begin
	neii  = getlineflux(contag[i], 'neII',/quiet)
	neiii = getlineflux(contag[i], 'neIII',/quiet)
	neii = neii[0] & neiii = neiii[0]

	if neii gt 0 and neiii gt 0 then con_ne_flux[i] = neii + neiii $
		else if neii gt 0 and neiii eq -1 then con_ne_flux[i] = neii $
		else if neiii gt 0 and neii eq -1 then con_ne_flux[i] = neiii

endfor

; Find luminosity distances

adl = float(transpose(archdat('dl')))
odl = float(transpose(ohmdat('dl')))

ohm_dl = [adl, odl]
con_dl = float(transpose(condat('dl')))

; Compute luminosities and errors

ohm_neon_lum = 4d * !dpi * (ohm_dl * mpc2cm)^2 * ohm_ne_flux * 1d7
con_neon_lum = 4d * !dpi * (con_dl * mpc2cm)^2 * con_ne_flux * 1d7

ohm_sfr_neon = 4.73d-41 * ohm_neon_lum
con_sfr_neon = 4.73d-41 * con_neon_lum

ohm_sfr_pah = 1.18d-41 * ohm_pah
con_sfr_pah = 1.18d-41 * con_pah

ohm_sfr_lir = 4.5d-44 * ohm_iras
con_sfr_lir = 4.5d-44 * con_iras

print,[transpose(ohmobj[sort(ohmobj)]),transpose(string(ohm_sfr_neon[sort(ohmobj)],format='(i4)'))]
print,''
print,[transpose(conobj[sort(conobj)]),transpose(string(con_sfr_neon[sort(conobj)],format='(i4)'))]

print,[transpose(ohmobj[sort(ohmobj)]),transpose(string(ohm_sfr_pah[sort(ohmobj)],format='(i4)'))]
print,''
print,[transpose(conobj[sort(conobj)]),transpose(string(con_sfr_pah[sort(conobj)],format='(i4)'))]

print,[transpose(ohmobj[sort(ohmobj)]),transpose(string(ohm_sfr_lir[sort(ohmobj)],format='(i4)'))]
print,''
print,[transpose(conobj[sort(conobj)]),transpose(string(con_sfr_lir[sort(conobj)],format='(i4)'))]

ohm_good = where(finite(ohm_sfr_neon) eq 1 and finite(ohm_sfr_pah) eq 1 and finite(ohm_sfr_lir) eq 1)
con_good = where(finite(con_sfr_neon) eq 1 and finite(con_sfr_pah) eq 1 and finite(con_sfr_lir) eq 1)

; Plot

!p.multi = [0,3,1]

red = fsc_color("Red")
blue = fsc_color("Blue")

plot, ohm_sfr_neon, ohm_sfr_pah, /nodata, xtit='Neon SFR', ytit='PAH SFR', charsize=3
oplot, ohm_sfr_neon, ohm_sfr_pah, color=red, psym=symcat(9)
oplot, con_sfr_neon, con_sfr_pah, color=blue, psym=symcat(7)
oplot, indgen(1d4), linestyle=1
xyouts, 100, 350, string(correlate(ohm_sfr_neon[ohm_good], ohm_sfr_pah[ohm_good]),format='(f5.2)'), color=red, charsize=2
xyouts, 100, 300, string(correlate(con_sfr_neon[con_good], con_sfr_pah[con_good]),format='(f5.2)'), color=blue, charsize=2

plot, ohm_sfr_neon, ohm_sfr_lir, /nodata, xtit='Neon SFR', ytit='L_IR SFR', charsize=3
oplot, ohm_sfr_neon, ohm_sfr_lir, color=red, psym=symcat(9)
oplot, con_sfr_neon, con_sfr_lir, color=blue, psym=symcat(7)
oplot, indgen(1d4), linestyle=1
xyouts, 100, 1000, string(correlate(ohm_sfr_neon[ohm_good], ohm_sfr_lir[ohm_good]),format='(f5.2)'), color=red, charsize=2
xyouts, 100, 900, string(correlate(con_sfr_neon[con_good], con_sfr_lir[con_good]),format='(f5.2)'), color=blue, charsize=2

plot, ohm_sfr_lir, ohm_sfr_pah, /nodata, xtit='L_IR SFR', ytit='PAH SFR', charsize=3
oplot, ohm_sfr_lir, ohm_sfr_pah, color=red, psym=symcat(9)
oplot, con_sfr_lir, con_sfr_pah, color=blue, psym=symcat(7)
oplot, indgen(1d4), linestyle=1
xyouts, 100, 350, string(correlate(ohm_sfr_lir[ohm_good], ohm_sfr_pah[ohm_good]),format='(f5.2)'), color=red, charsize=2
xyouts, 100, 300, string(correlate(con_sfr_lir[con_good], con_sfr_pah[con_good]),format='(f5.2)'), color=blue, charsize=2

print,''
print,'Mean OHM SFR (neon): ', mean(ohm_sfr_neon[ohm_good]),' +- ', stddev(ohm_sfr_neon[ohm_good])
print,'Mean OHM SFR (PAH): ',  mean(ohm_sfr_pah[ohm_good]),' +- ',  stddev(ohm_sfr_pah[ohm_good])
print,'Mean OHM SFR (LIR): ',  mean(ohm_sfr_lir[ohm_good]),' +- ',  stddev(ohm_sfr_lir[ohm_good])
print,''
print,'Mean non-masing SFR (neon): ', mean(con_sfr_neon[con_good]),' +- ', stddev(con_sfr_neon[con_good])
print,'Mean non-masing SFR (PAH): ',  mean(con_sfr_pah[con_good]),' +- ',  stddev(con_sfr_pah[con_good])
print,'Mean non-masing SFR (LIR): ',  mean(con_sfr_lir[con_good]),' +- ',  stddev(con_sfr_lir[con_good])
print,''

end
