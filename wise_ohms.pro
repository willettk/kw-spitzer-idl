
;+
; NAME:
;       
;	WISE_OHMS
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

file='~/Desktop/wise_allsky.wise_allsky_4band_p3as_psd22202.tbl.tbl.tbl'
limfile = '~/Desktop/wise_allsky.wise_allsky_4band_p3as_psd10344.tbl.tbl.tbl'
o = read_ipac_table(file)
l = read_ipac_table(limfile)

!p.multi=[0,2,2]

w1 = o.w1mpro
w2 = o.w2mpro
w3 = o.w3mpro
w4 = o.w4mpro

w1lim = l.w1mpro
w2lim = l.w2mpro
w3lim = l.w3mpro
w4lim = l.w4mpro

cs = 1.5

cgplot, w2-w3, w1-w2, $
	charsize=cs, $
	psym=16, $
	color='red', $
	xtitle='W2 - W3', $
	ytitle='W1 - W2', $
	xr=[-1,5], $
	yr=[-1,2]

cgplot, w2lim-w3lim, w1lim-w2lim, $
	/overplot,$
	psym=9, $
	color='blue'

cgplot, w3-w4, w2-w3, $
	charsize=cs, $
	psym=16, $
	color='red', $
	xtitle='W3 - W4', $
	ytitle='W2 - W3'

cgplot, w3lim-w4lim, w2lim-w3lim, $
	/overplot,$
	psym=9, $
	color='blue'

cgplot, w2, w1-w2, $
	charsize=cs, $
	psym=16, $
	color='red', $
	yr=[-1,2] ,$
	xtitle='W2', $
	ytitle='W1 - W2'

;cgplot, w2lim, w1lim-w2lim, $
;	/overplot,$
;	psym=9, $
;	color='blue'

cgplot, w3, w2-w3, $
	charsize=cs, $
	psym=16, $
	color='red', $
	yr=[-1,5] ,$
	xtitle='W3', $
	ytitle='W2 - W3'

;cgplot, w3lim, w2lim-w3lim, $
;	/overplot,$
;	psym=9, $
;	color='blue'

kstwo, w1-w2, w1lim-w2lim, d12, p12
kstwo, w2-w3, w2lim-w3lim, d23, p23
kstwo, w3-w4, w3lim-w4lim, d34, p34

print,'W1-W2',sqrt(2d) * inverf(1.-p12)
print,'W2-W3',sqrt(2d) * inverf(1.-p23)
print,'W3-W4',sqrt(2d) * inverf(1.-p34)


end
