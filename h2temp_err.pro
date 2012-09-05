function h2temp_err, f1, f2, del_f1, del_f2, dl, hot=hot

; Quick IDL function to compute the error in the calculated H2 temp
; using only two points in the ln(Nj / gj) vs E / kb space. 

; Fluxes and errors are in erg/s/cm^2, D_L in Mpc

; Assume all lines (as in lo-res) are H2 S(1) and S(3)

kb = 1.38d-16

if keyword_set(hot) then begin

	e1 = 2504*kb
	e2 = 7197*kb
	
	a1 = 9.84d-9
	a2 = 2.00d-7
	
	g1 = 33
	g2 = 57

endif else begin

	e1 = 1015*kb
	e2 = 2504*kb
	
	a1 = 4.8d-10
	a2 = 9.84d-9
	
	g1 = 21
	g2 = 33
	
endelse

dl = dl * 3.09d24

l1 = 4d * !dpi * dl^2 * f1 * 1d7
l2 = 4d * !dpi * dl^2 * f2 * 1d7
del_l1 = 4d * !dpi * dl^2 * del_f1 * 1d7
del_l2 = 4d * !dpi * dl^2 * del_f2 * 1d7

dt_dl1 = ((e2 - e1) / (l1)) / (alog((l2 * a1 * e1 * g1)/(l1 * a2 * e2 * g2)))^2  / kb
dt_dl2 = (-(e2 - e1) / (l2)) / (alog((l2 * a1 * e1* g1)/(l1 * a2 * e2 * g2)))^2 / kb

del_temp = sqrt(del_l1^2 * dt_dl1^2 + del_l2^2 * dt_dl2^2)


return, del_temp

end
