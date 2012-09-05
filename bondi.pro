;+
; NAME:
;       
;	BONDI.pro
;
; PURPOSE:
;
;	Compute Bondi power vs. radio jet power + luminosity for CSO sample
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
;	IDL> bondi
;
; NOTES:
;
;	Tests the model in Allen et al. (2006) for equipartition of the jet power and Bondi accretion luminosity
;
; REVISION HISTORY
;
;       Written by K. Willett                Oct 09
;-

pro bondi, ps=ps, stop=stop, label=label

; P_Bondi (slow, ADAF accretion onto a black hole that would not necessarily produce fine-structure emission lines)

lambda = 0.25		; numerical coefficient depending on adiabatic index of accreting gas
eta = 0.1		; efficiency of accretion
gamma_1 = 5./3.		; adiabatic index of accreting gas
G = 6.670d-8		; gravitaitional constant [dyne cm^2/g^2]
c = 299492758. * 1d2	; speed of light [cm/s]
mu = 0.62		; mean atomic weight of gas, assuming cosmic abundances
mp = 1.673d-24		; mass of a proton [g]
k = 1.381d-16		; Boltzmann's constant [erg/K]
M_sun = 1.989d33	; mass of the sun [g]
pc2cm = 3.086d18	; 1 parsec in cm
keV2kelvin = 1.16d7	; keV to Kelvin conversion

mbh_default = 1d8

;M_BH = [mbh_default, 1.5d8, 10^9.23, mbh_default, mbh_default, mbh_default, mbh_default, 1d7, 10^8.62, mbh_default] * M_sun	; mass of central black hole [g]; assume 10^8 M_sun if no other information is available

	; New method - use bulge luminosity-M_BH relation from Bentz et al. (2009) with SDSS photometry

	M_BH = 10.^[8.78, 8.74, 8.39, 8.59, 8.81, 8.81, 8.22, 8.05, 8.48, 8.54] * M_sun	; A_V corrected from NED Schlegel

	; V magnitude for NGC 5793 is from de Vaucouleurs catalog (m_V = 13.22)
	; B mag for PKS 1718-649 is m_B = 13.16; assume typical elliptical color of B - V = 0.9
	; 	Bulge luminosity BH mass of 8.48 is close to [O IV] dispersion mass of 8.62 for PKS 1718-649
	; H absolute mag for 1946+70 is -23.77 (Perlman+01)

rho = 1d0 * mp		; density of gas [g/cm^3]
T = 1d0 * kev2kelvin	; temperature of gas [K] (assume ~ 1 keV)

cs = sqrt(gamma_1 * k * T / (mu * mp))		; sound speed in the gas [cm/s]

pbondi = 4d * !dpi * lambda * G^2 * M_BH^2 * rho * eta * c^2 / cs^3	; Bondi power (using L_bulge M_BH)

	; Estimate systematic errors in P_Bondi

	dPb_drho = 4d * !dpi * lambda * G^2 * M_BH^2 * eta * c^2 / cs^3

	dPb_dT = 4d * !dpi * lambda * G^2 * M_BH^2 * rho * eta * c^2 * (-1.5) * $
		(gamma_1 * k * T / (mu * mp))^(-5./2) * gamma_1 * k / (mu * mp)

	dPb_dMBH = 4d * !dpi * lambda * G^2 * rho * eta * c^2 / cs^3 * 2d * M_BH

	sigma_rho = 0.92 * mp				; Errors are from scatter for ellipticals in Allen sample
	sigma_T = 0.21 * keV2kelvin
	sigma_MBH = mbh_default * 0.23 * M_sun

	sigma_pbondi = sqrt(sigma_rho^2 * dPb_drho^2 + sigma_T^2 * dPb_dT^2 + sigma_MBH^2 * dPb_dMBH^2)

; Power from the radio galaxy - mechanical work + radio luminosity

rw = [ 13.5,  4., 0.5,  2.,    5., 2.,  15.,  5.,  0.25,  10.] * pc2cm	; radius along prolate axis of jet
rl = [112.1, 42., 2.08, 9.6, 146., 8., 170., 14.4, 1.75,  58.] * pc2cm	; distance between jet hotspots

tage = [550, 502, 34, 188, 1700, 92, 1d4, 1d5, 1d5, 4000] * 3.16d7	; Kinematic age of CSO [sec]

p = k * T * rho / (mu * mp)			; pressure in the jet [dyne/cm^2]
V = rw^2 * !dpi * rl				; volume of cavity excavated by jet [cm^3] - assume cylindrical shape
gamma_2 = 4/3.					; adiabatic index of relativistic plasma

ejet = gamma_2 / (gamma_2 - 1) * p * V		; energy in the jets

pjet = ejet / tage

;l1420 = [25.34,25.05,23.10,25.02,26.30,25.08,26.31,23.55,24.25,25.39]	; radio power at 1.4 GHz [W/Hz]

pradio = [3.0d42, 4.2d42, 5.0d41, 3.4d42, 5.1d43, 1.2d43, 9.7d44, 5.9d40, 7.6d41, 7.5d42]	 ; From LRADIO.pro

	; Systematic errors in jet power

	dPj_drho = gamma_2 / (gamma_2 - 1) * V * k * T / (mu * mp) / tage

	dPj_dT = gamma_2 / (gamma_2 - 1) * V * k * rho / (mu * mp) / tage

	dPj_drw = gamma_2 / (gamma_2 - 1) * p * 2d * rw * !dpi * rl / tage

	dPj_drl = gamma_2 / (gamma_2 - 1) * p * rw^2 * !dpi / tage

	dPj_dtage = gamma_2 / (gamma_2 - 1) * V * p / (-1 * tage^2)

	sigma_rw = 0.1 * rw		; Assume 10%
	sigma_rl = 0.1 * rw 		; Assume 10%
	sigma_tage = 0.25 * tage	; Assume 25%

	sigma_pjet = sqrt(sigma_rho^2 * dPj_drho^2 + sigma_T^2 * dPj_dT^2 + $
		sigma_rw^2 * dPj_drw^2 + sigma_rl^2 * dPj_drl^2 + sigma_tage^2 * dPj_dtage^2)

	; Measured error in radio luminosity

	sigma_pradio = 0.1 * pradio		; Assume 10% error

; Plot

; Data from Allen et al. (2006)

allen_jet = alog10([6.60  +  3.62  , 0.868 +  0.658 , 0.445 +  0.362 , 1.26  +  2.18  , 0.0727+  0.0829, 0.0147+  0.0152, 0.313 +  0.478 , 0.0367+  0.0374, 0.741 +  0.838])

jet_err = alog10(sqrt([2.91^2 + 1.70^2, 0.371^2 + 0.268^2, 0.181^2+ 0.150^2, 0.84^2+ 1.44^2, 0.0274^2+ 0.0307^2, $
	0.0060^2+0.0060^2, 0.177^2+ 0.244^2, 0.0185^2+ 0.0188^2, 0.316^2+0.353^2]))

allen_bondi = [1.41, 0.69, 0.79, 1.16, 0.37, -0.71, 0.40, -0.15, 0.49] 
lower_err   = [0.09, 0.29, 0.23, 0.40, 0.21, 0.24, 0.55, 0.40, 0.26]
upper_err   = [0.09, 0.30, 0.25, 0.28, 0.22, 0.24, 0.56, 0.43, 0.34]

!p.multi=[0,1,1]

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename = '~/Astronomy/Research/Spitzer/cso/papers/bondi.ps', $
		/color, /encap, xoff=1, yoff=1, /portrait
	cs = 2.0
	ls = 5
	defcolor=fsc_color("Black")
	jetcolor = fsc_color("Blue")
;	!p.font=0
endif else begin
	cs = 2
	ls = 1
	defcolor=fsc_color("White")
	jetcolor = fsc_color("Cyan")
endelse
	
plot, indgen(10), $
	/nodata, $
;	/iso, $
	xr = [-2.5,3.0], /xstyle, $
	yr = [-2,2], /ystyle, $
	charsize = cs, $
	charthick = ls, $
	thick = ls, $
	xthick = ls, $
	ythick = ls, $
	xtitle=' log (P!Ijet+radio!N / 10!E43!N erg s!E-1!N)', $
	ytitle=' log (P!IBondi!N / 10!E43!N erg s!E-1!N)'

oploterror, $
	allen_jet, $
	allen_bondi, $
	allen_jet - alog10(10.^allen_jet - 10.^jet_err), $
	lower_err, $
	errcolor=fsc_color("Dark Grey"), $
	psym=symcat(9), /lobar, /nohat, thick = ls

oploterror, $
	allen_jet, $
	allen_bondi, $
	allen_jet - alog10(10.^allen_jet + 10.^jet_err), $
	upper_err, $
	color=fsc_color("Dark Grey"), $
	errcolor=fsc_color("Dark Grey"), $
	psym=symcat(9), /hibar, /nohat, thick = ls

xarr = fillarr(0.1,-10,10)
a = 0.65
b = 0.77

oplot, xarr, a + b * xarr, linestyle=2, thick = ls, color=fsc_color("Dark Grey")
;oplot, xarr, (xarr - a) / b, linestyle=2, thick = ls

; CSO data

csoind = [0,1,2,4,5,6,8,9]

oploterror, $
	alog10((pjet[csoind] + pradio[csoind]) / 1d43), $
	alog10(pbondi[csoind] / 1d43), $
	alog10((pjet[csoind] + pradio[csoind]) / 1d43) - alog10((pjet[csoind]+pradio[csoind] - sqrt(sigma_pjet[csoind]^2 + sigma_pradio[csoind]^2))/1d43), $
	alog10(pbondi[csoind]/1d43) - alog10((pbondi[csoind] + sigma_pbondi[csoind]) / 1d43), $
	color=fsc_color("Black"), psym=symcat(16), /nohat, errcolor=fsc_color("Black"), $
	thick = ls

; Draw line indicating fraction of power coming from jet vs. radio luminosity

;for i=0, n_elements(csoind)-1 do $
;	plots, $
;		[alog10((pjet[csoind[i]] + pradio[csoind[i]]) / 1d43), $
;		alog10((pradio[csoind[i]]) / 1d43)], $
;		[alog10(pbondi[csoind[i]] / 1d43), $
;		alog10(pbondi[csoind[i]] / 1d43)], $
;		color=jetcolor, thick=6, linestyle=0

;cso_agelimits = [4,6,7,8]
cso_agelimits = [4,6,8]

arrow, $
	alog10((pjet[cso_agelimits] + pradio[cso_agelimits]) / 1d43), $
	alog10(pbondi[cso_agelimits] / 1d43), $
	alog10((pjet[cso_agelimits] + pradio[cso_agelimits]) / 1d43) + 0.5, $
	alog10(pbondi[cso_agelimits] / 1d43), $
	color=fsc_color("Black"), /data, $
	hthick = ls, thick = ls

;legend, /top, /left, $
;	['CSOs', 'ellipticals'], $
;	color=[fsc_color("Black"), defcolor], $
;	psym=[16,16], $
;	charthick = ls, $
;	charsize = 1.5

; Label CSOs that deviate from Allen relation

xyouts, 1.2, -1.5, $
	'PKS 1413+135', $
	color=fsc_color("Black"), $
	charsize = 1.0, $
	charthick = 4

xyouts, 1.0, 0.7, $
	'4C 12.50', $
	color=fsc_color("Black"), $
	charsize = 1.0, $
	charthick = 4

if keyword_set(label) then begin

	obj = csodat('obj',/rasort)
	obj = obj[0,*]
	
	xyouts, $
	alog10((pjet[csoind] + pradio[csoind]) / 1d43), $
	alog10(pbondi[csoind] / 1d43) -0.2, $
	obj[csoind], $
	color=fsc_color("Black"), $
	charsize = cs, charthick = ls

endif

if keyword_set(ps) then begin
;	!p.font=-1
	device, /close
	set_plot,'x'
endif

if keyword_set(stop) then stop

end
