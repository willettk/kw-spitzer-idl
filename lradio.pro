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

;function lradio, d, sc, nuc, nu0, nu1, nu2, alpha_low, alpha_high

mpc2cm = 3.086d24
jy2cgs = 1d-23

;l = 4d * !dpi * (d * mpc2cm)^2 * sc * jy2cgs / nuc * $
;	(1d / (alpha_low + 1) * (nu1^(alpha_low+1) - nu0^(alpha_low+1)) + $
;	1d / (alpha_high + 1) * (nu2^(alpha_high+1) - nu1^(alpha_high+1)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 4C +31.04
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Data from NED

d = 264
nu_low = 178d6	; 408 MHz
s_low = 2.4 	; Jy
nu_high = 4.8d9	; 4.8 GHz
s_high = 1.25	; Jy		Giroletti et al. 2003

alpha_low = 2.5
alpha_high = -1.1		; gir03 - steepest index in the lobe portions

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 300d6			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'4C +31.04'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 4C +37.11
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Data from NED

d = 242
nu_low = 178d6	; 408 MHz		NED
s_low = 3.2 	; Jy
nu_high = 4.8d9	; 4.8 GHz		NED
s_high = 1.04	; Jy		

alpha_high = -0.3		; NED

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 300d6			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'4C +37.11'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 1146+59
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 48.2
nu_low = 365d6	; 408 MHz
s_low = 0.327 	; Jy
nu_high = 8.4d9	; 4.8 GHz
s_high = 0.516	; Jy

alpha_high = -0.34

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 4d9			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'1146+59'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 1245+676
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 495
nu_low = 1.4d9	; 408 MHz
s_low = 0.252 	; Jy
nu_high = 4.8d9	; 4.8 GHz
s_high = 0.188	; Jy

alpha_high = -0.24

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 1.4d9			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'1245+676'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 4C +12.50
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 571
nu_low = 178d6	; 408 MHz
s_low = 4.6 	; Jy
nu_high = 4.8d9	; 4.8 GHz
s_high = 2.93	; Jy

alpha_high = -0.46

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 408d6			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'4C +12.50'
print,string(l, format='(e7.1)')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; OQ 208
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 349
nu_low = 1.4d9	; 408 MHz
s_low = 0.910 	; Jy
nu_high = 15d9	; 4.8 GHz
s_high = 1.139	; Jy

alpha_high = -0.92

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 5d9			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'OQ 208'
print,string(l, format='(e7.1)')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PKS 1413+135
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 1244
nu_low = 8.4d9	; 408 MHz
s_low = 0.61 	; Jy
nu_high = 230d9	; 4.8 GHz
s_high = 1.50	; Jy

alpha_high = -0.79

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 37d9			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'PKS 1413+135'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NGC 5793
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 53
nu_low = 1.4d9	; 408 MHz
s_low = 1.047 	; Jy
nu_high = 4.8d9	; 4.8 GHz
s_high = 0.508	; Jy

alpha_high = -0.59

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 1.4d9			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'NGC 5793'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PKS 1718-649
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Data from Tingay et al. (1997)

d = 62.6
nu_low = 408d6	; 408 MHz
s_low = 2.47 	; Jy
nu_high = 4.8d9	; 4.8 GHz
s_high = 4.73	; Jy

alpha_low  = 0.36
alpha_high = -0.65

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 408d6			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'PKS 1718-649'
print,string(l, format='(e7.1)')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 1946+70
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

d = 460
nu_low = 1.4d9	; 408 MHz
s_low = 0.970 	; Jy
nu_high = 4.8d9	; 4.8 GHz
s_high = 0.643	; Jy

alpha_low = 2.5
alpha_high = -0.4

nu0 = 1d7			; Frequency at which synchrotron absorption effectively cuts off low frequency emission
nu1 = 1.4d9			; Peak frequency
nu2 = 1d11			; Break frequency?

l = 4d * !dpi * (d * mpc2cm)^2 * $
	(nu_low * s_low * jy2cgs / (alpha_low + 1) * ((nu1/nu_low)^(alpha_low + 1) - (nu0/nu_low)^(alpha_low + 1)) + $
	nu_high * s_high * jy2cgs / (alpha_high + 1) * ((nu2/nu_high)^(alpha_high + 1) - (nu1/nu_high)^(alpha_high + 1)))

print,''
print,'1946+70'
print,string(l, format='(e7.1)')



print,''

end
