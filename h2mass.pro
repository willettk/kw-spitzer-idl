function h2mass, fj, tex, dl, line, msun = msun, quiet = quiet
;+
; NAME:
;       H2MASS
;
; PURPOSE:
; 	Compute mass of molecular hydrogen using luminosity of mid-IR line and T_ex
;
; INPUTS:
;
;	FJ - flux of H2 line [W/cm^2]
;
;	TEX - excitation temperature [K]
;
;	DL - luminosity distance of target [Mpc]
;
;	LINE - rotational level of line to use
;
;
; OUTPUTS:
;
;	MTOT - mass of molecular hydrogen [g]
;
; KEYWORDS:
;
;	MSUN - gives output in units of 10^7 M_sun
;
; EXAMPLE:
;	IDL> mass = h2mass(1d-21, 300, 530, 3, /msun)
;
; REQUIRES:
;
; REVISION HISTORY
;       Written by K. Willett                Sep 2007
;-

; Define physical constants (all cgs)

c = 3d10
h = 6.626d-27
mh2 = 3.345d-24
kb = 1.381d-16
solarmass = 1.99d40

; Constants for H2 lines

phi_orthotopara = 4d/3d
phi_paratoortho = 4d

lambda0 = 28.221d * 1d-4
lambda1 = 17.035d * 1d-4
lambda2 = 12.279d * 1d-4
lambda3 = 9.6649d * 1d-4
lambda4 = 8.0258d * 1d-4
lambda5 = 6.9091d * 1d-4
lambda6 = 6.1088d * 1d-4
lambda7 = 5.5115d * 1d-4

a0 = 2.94d-11
a1 = 4.76d-10
a2 = 2.75d-9
a3 = 9.83d-9
a4 = 2.64d-8
a5 = 5.88d-8
a6 = 1.14d-7
a7 = 2.00d-7

g0 = 5d
g1 = 21d
g2 = 9d
g3 = 33d
g4 = 13d
g5 = 45d
g6 = 17d
g7 = 57d

e0 = 510d 
e1 = 1015d 
e2 = 1682d 
e3 = 2504d 
e4 = 3474d 
e5 = 4586d 
e6 = 5879d 		; Roussel gives 5828
e7 = 7197d 


; Calculate ortho-H2 partition function

case line of
	0: begin
		main_e = e0
		main_g = g0
		main_a = a0
		main_dele = h * c / lambda0
	   end
	1: begin
		main_e = e1
		main_g = g1
		main_a = a1
		main_dele = h * c / lambda1
	   end
	2: begin
		main_e = e2
		main_g = g2
		main_a = a2
		main_dele = h * c / lambda2
	   end
	3: begin
		main_e = e3
		main_g = g3
		main_a = a3
		main_dele = h * c / lambda3
	   end
	4: begin
		main_e = e4
		main_g = g4
		main_a = a4
		main_dele = h * c / lambda4
	   end
	5: begin
		main_e = e5
		main_g = g5
		main_a = a5
		main_dele = h * c / lambda5
	   end
	6: begin
		main_e = e6
		main_g = g6
		main_a = a6
		main_dele = h * c / lambda6
	   end
	7: begin
		main_e = e7
		main_g = g7
		main_a = a7
		main_dele = h * c / lambda7
	   end
endcase

if line mod 2 eq 0 then begin
	phi = phi_paratoortho 	
	zbottom = g0 * exp(-1d * e0 / tex) + g2 * exp(-1d * e2 / tex) $
		+ g4 * exp(-1d * e4 / tex) + g6 * exp(-1d * e6 / tex)
	deltemp_term = exp(-1d * (e0 + e2 + e4 + e6 - main_e) / tex)  * $
		(exp((e2 + e4 + e6)/tex) * (e0 - main_e) * g0 + $
		 exp((e0 + e4 + e6)/tex) * (e2 - main_e) * g2 + $
		 exp((e0 + e2 + e6)/tex) * (e4 - main_e) * g4 + $
		 exp((e0 + e2 + e4)/tex) * (e6 - main_e) * g6)
	if not keyword_set(quiet) then print,'Para line'
endif else begin 
	phi = phi_orthotopara
	zbottom = g1 * exp(-1d * e1 / tex) + g3 * exp(-1d * e3 / tex) $
		+ g5 * exp(-1d * e5 / tex) + g7 * exp(-1d * e7 / tex)
	deltemp_term = exp(-1d * (e1 + e3 + e5 + e7 - main_e) / tex)  * $
		(exp((e3 + e5 + e7)/tex) * (e1 - main_e) * g1 + $
		 exp((e1 + e5 + e7)/tex) * (e3 - main_e) * g3 + $
		 exp((e1 + e3 + e7)/tex) * (e5 - main_e) * g5 + $
		 exp((e1 + e3 + e5)/tex) * (e7 - main_e) * g7)
	if not keyword_set(quiet) then print,'Ortho line'
endelse

ztop = main_g * exp(-1d * main_e / tex)
z = ztop / zbottom

; Convert luminosity to flux

fj_cgs = fj * 1d7			; Convert flux from W/cm^2 to erg/cm^2/s
dl_mpc = dl * 3.086d24		; Convert luminosity distance from Mpc to cm
lj = 4d * !dpi * dl_mpc^2 * fj_cgs	; Flux to luminosity

; Calculate gas mass

mtot = mh2 * phi * lj  / (main_a * main_dele * z)

df = 0.38d-21 * 1d7
dt = 100.

errterm1 = 4d * !dpi * phi * mh2 * df * dl_mpc^2 / (main_a * main_dele * z)

errterm2 = 4d * !dpi * phi * mh2 * fj_cgs * dl_mpc^2 * dt * deltemp_term / $
	(main_a * main_dele * main_g * tex^2)

mtot_err = sqrt(errterm1^2 + errterm2^2)

if keyword_set(msun) then mtot = mtot / solarmass ; Gives mass in units of 10^7 M_sun

; Return mass in grams

;if not keyword_set(quiet) then print, mtot_err / solarmass

return, mtot

end
