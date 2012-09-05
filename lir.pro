;+
; LIR
;
; Computes the total IR luminosity based on all IRAS fluxes (12, 25, 60, 100) um
;
; Takes fluxes in Jy, D_L in Mpc; outputs result in L_sun
;
; Format: lir(f12, f25, f60, f100, dl)
; Example: print,lir(0.22, 0.45, 2.56, 3.45, 104)
;
; Formula comes from Sanders & Mirabel (1996), ARAA
;-

function lir, f12, f25, f60, f100, dl

mpc2cm = 3.09d24 ; cm/Mpc
mks2cgs = 1d3	; 1 W m^-2 = 10^3 erg s^-1 cm^-2
lsun = 3.862d33 ; erg s^-1

ir_flux=1.8d-14 * (13.48 * f12 + 5.16 * f25 + 2.58 * f60 + f100)	; W m^-2

ir_lum = (4d * !dpi * (dl*mpc2cm)^2 * ir_flux * mks2cgs) / lsun

return, alog10(ir_lum)

end
