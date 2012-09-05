
;+
; NAME:
;       
;	IR8_OHMS
;
; PURPOSE:
;
;	Plot the IR8 diagnostic of Elbaz et al. (2011) for the IRS sample of OHMs and non-masing galaxies
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
;       Written by K. Willett                Aug 11
;-

mnames = transpose(ohmdat('tag'))
badohm = where(mnames eq 'mega034')
goodohm = setdifference(indgen(n_elements(mnames)),badohm)
mnames = mnames(goodohm)
nm = n_elements(mnames)

anames = transpose(archdat('tag'))
na = n_elements(anames)

cnames = condat('tag')
badcon = where(cnames eq 'control033')
goodcon = setdifference(indgen(n_elements(cnames)),badcon)
cnames = cnames(goodcon)
nc = n_elements(cnames)

; Retrieve L_IR for the Spitzer OHMs and non-masing galaxies

mega_lir = transpose(ohmdat('lir'))
arch_lir = transpose(archdat('lir'))

ohm_lir = [mega_lir[goodohm],arch_lir]

lir_ohm = 10.^float(ohm_lir)

con_lir = transpose(condat('lir'))
lir_con = 10.^float(con_lir[goodcon])

; Calculate L8 based on rest-frame 8 um continuum in IRS spectra

; Loop over OHMs
	 
allohms = [mnames,anames]
no = n_elements(allohms)

fnu8_ohm = fltarr(no)
fnu8_ohm_err = fltarr(no)
dl_ohm = fltarr(no)
fnu8_con = fltarr(nc)
fnu8_con_err = fltarr(nc)
dl_con = fltarr(nc)

c = 299792.458 * 1d5
nu8 = c / (8d * 1d-4)
mpc2cm = 3.086d24
lsun = 3.862d33
jy2cgs = 1d-23

for i = 0, no - 1 do begin

	; Restore data

	tag, allohms[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+allohms[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	fnu8_ohm[i] = flux[closeto(wave,8.0)]
	fnu8_ohm_err[i] = err[closeto(wave,8.0)]
	dl_ohm[i] = lumdist(sed.redshift,/wmap5,/silent) * mpc2cm

endfor

for i = 0, nc - 1 do begin

	; Restore data

	tag, cnames[i], dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames[i]+'.sav'

	wave = sed.wave_lr
	flux = sed.flux_lr
	err =  sed.err_lr

	fnu8_con[i] = flux[closeto(wave,8.0)]
	fnu8_con_err[i] = err[closeto(wave,8.0)]
	dl_con[i] = lumdist(sed.redshift,/wmap5,/silent) * mpc2cm

endfor

l8_ohm     = nu8 * fnu8_ohm * (double(dl_ohm))^2 * 4d * !dpi * jy2cgs / lsun
l8_ohm_err = nu8 * fnu8_ohm_err * (double(dl_ohm))^2 * 4d * !dpi * jy2cgs / lsun
l8_con     = nu8 * fnu8_con * (double(dl_con))^2 * 4d * !dpi * jy2cgs / lsun
l8_con_err = nu8 * fnu8_con_err * (double(dl_con))^2 * 4d * !dpi * jy2cgs / lsun

ir8_ohm = lir_ohm/l8_ohm
ir8_con = lir_con/l8_con

; Replicate Fig. 14 from Elbaz paper

ps_start, file='~/Astronomy/Research/Spitzer/plots/ir8_ohms.ps', /quiet
; Top

x0 = 0.20
x1 = 0.9
y0 = 0.15
y1 = 0.4
y2 = 0.9

cgplot, l8_ohm, lir_ohm, psym=16, color="Red", $
	position=[x0,y1,x1,y2], $
	/ylog, /xlog, $
	xrange=[1d8,2d12], /xstyle, $
	yrange=[1d9,1d13], /ystyle, $
	ytitle='L!IIR!N!Etot!N [L!I'+sunsymbol()+'!N]', $
	xtickname=replicate(' ',10)

cgplot, l8_con, lir_con, psym=15, color="Blue", /over

; Local relationship for IR8 follows a Gaussian distribution with mu = 3.9, sigma = 1.25

xarr = (findgen(1d5) + 1) * 1d8
cgplot, /over, xarr, xarr*3.9, color="Black", thick=2, linestyle=0
cgplot, /over, xarr, xarr*(3.9-1.25), color="Black", thick=2, linestyle=1
cgplot, /over, xarr, xarr*(3.9+1.25), color="Black", thick=2, linestyle=1

al_legend, /top, /left, psym=[16,15], color=['Red','Blue'], ['OHMs','Non-masing'], charsize=1.2

; Bottom

cgplot, l8_ohm, ir8_ohm, psym=16, color="Red", $
	/noerase, $
	position=[x0,y0,x1,y1], $
	/ylog, /xlog, $
	xrange=[1d8,2d12], /xstyle, $
	yrange=[5d-1,1d2], /ystyle, $
	xtitle='L!I8!N [L!I'+sunsymbol()+'!N]', $
	ytitle='IR8'
cgplot, l8_con, ir8_con, psym=15, color="Blue", /over
cgplots, 10.^!x.crange, replicate(mean(ir8_ohm),2), color="Red", thick=2, linestyle=2
cgplots, 10.^!x.crange, replicate(mean(ir8_con),2), color="Blue", thick=2, linestyle=2
;cgplots, 10.^!x.crange, replicate(mean(ir8_ohm),2)+stddev(ir8_ohm), color="Red", thick=2, linestyle=1
;cgplots, 10.^!x.crange, replicate(mean(ir8_ohm),2)-stddev(ir8_ohm), color="Red", thick=2, linestyle=1
;cgplots, 10.^!x.crange, replicate(mean(ir8_con),2)+stddev(ir8_con), color="Blue", thick=2, linestyle=1
;cgplots, 10.^!x.crange, replicate(mean(ir8_con),2)-stddev(ir8_con), color="Blue", thick=2, linestyle=1

; Local relationship for IR8 follows a Gaussian distribution with mu = 3.9, sigma = 1.25

cgplots, 10.^!x.crange, [3.9, 3.9], color="Black", thick=2
cgplots, 10.^!x.crange, [3.9, 3.9]-1.25, color="Black", thick=1, linestyle=1
cgplots, 10.^!x.crange, [3.9, 3.9]+1.25, color="Black", thick=1, linestyle=1

kstwo, ir8_ohm, ir8_con, d_ks, prob_ks
gauss_ir8 = sqrt(2d) * inverf(1d - prob_ks)

ps_end

print,''
print,'OHM: ',mean(ir8_ohm),stddev(ir8_ohm)
print,'Non-masing: ',mean(ir8_con),stddev(ir8_con)
print,'K-S test:', prob_ks, gauss_ir8
print,''

; Write data to an IDL file for James (Feb 8 2012)

l8_nonmasing  = l8_con
ir8_nonmasing = ir8_con
lir_nonmasing = lir_con
dl_nonmasing = dl_con

names_ohm = strarr(no)
for i=0,no-1 do begin
	targets, allohms[i], z, name, dl
	names_ohm[i] = name
endfor

names_nonmasing = strarr(nc)
for i=0,nc-1 do begin
	targets, cnames[i], z, name, dl
	names_nonmasing[i] = name
endfor

save, filename='~/Astronomy/Research/Spitzer/ir8_ohms.sav', $
	l8_ohm, ir8_ohm, lir_ohm, names_ohm, dl_ohm, $
	l8_nonmasing, ir8_nonmasing, lir_nonmasing, names_nonmasing, dl_nonmasing

stop

end
