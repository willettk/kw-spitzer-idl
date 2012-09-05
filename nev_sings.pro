;+
; NAME:
;       
;	NEV_SINGS
;
; PURPOSE:
;
;	Plot data from SINGS
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
;       Written by K. Willett                Aug 08
;-

fname = '~/Astronomy/Research/Spitzer/ddale/dale_NeV_Xray_SINGS2.csv'

readcol,fname,$
	name,$
	nev14,nev14err,nev24,nev24err,spitzarea,$
	softx,softerr,softrad,softtype,$
	hardx,harderr,hardrad,hardtype,$
	xflux,xerr,xrad,xtype,$
	lco,lhcn,hcncoratio,$
	oiv,oiverr,$
	format='a,f,f,f,f,f,f,f,f,a,f,f,f,a,f,f,f,a,f,f,f,f,f', /silent

nevind = where(nev14err ne -9.99)
nevlim = where(nev14err eq -9.99)

!p.multi = [0,1,1]

plot, nev14 / softx, nev24 / softx, $
	psym = 4, $
	/nodata, $
	title = 'SINGS data', $
	xtitle='[Ne V] 14 / X-ray', $
	ytitle='[Ne V] 24 / X-ray'

oploterror, $
	nev14[nevind]/softx[nevind], $
	nev24[nevind]/softx[nevind], $
	nev14err[nevind]/softx[nevind], $
	nev24err[nevind]/softx[nevlim], $
	psym = symcat(14)

oplot, $
	nev14[nevlim]/softx[nevlim], $
	nev24[nevlim]/softx[nevlim], $
	psym = 4

end
