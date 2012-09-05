pro texplot, ps = ps

; Make H2 T_ex for slide in Comps 2 defense

;j = 1,3,5,7
h = [1.24,0.75,0.75,0.66]*1d-21	; W/cm^2
herr = [0.34,0.09,0.07,0.24]*1d-21	; W/cm^2

hlines = ['H!I2!N S(1)', 'H!I2!N S(3)', 'H!I2!N S(5)', 'H!I2!N S(7)']
ej = [1015, 2504, 4586, 7197]				; Energy of upper state
aj = [4.8d-3, 9.84d-2, 0.588, 2.00]			; Einstein-A coefficient
lambda = [17.035,9.6649,6.9091,5.5115]			; Rest wavelength of transition
gj = [21,33,45,57]					; Statistical weight of upper level

plancks = 6.6261d-27					; Planck's constant
c = 3d10						; Speed of light [cm/s]
del_ej = plancks * c / lambda				; Energy of transition

; Linear fit for systems with a single excitation temperature

expr = 'p[0] + x*p[1]'
start = [10, 3d-3]
xarr = fillarr(1, 0, 8000)

nj     = h     / (aj * del_ej)
nj_err = herr / (aj * del_ej)

; PS output

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Comps2/images/texplot.ps', /landscape
	cs = 2
	lthick = 2.5
	cthick = 2.5
endif else begin
	cs = 2
	lthick = 1
	cthick = 1
endelse

!p.multi=[0,1,1]
plot, ej, alog(nj / gj), $
	xtitle = 'Upper energy level E!IJ!N/k!Ib!N [K]', $
	ytitle = 'ln(N!IJ!N / g!IJ!N) + constant', $
	psym = 4, $
	charsize = cs, $
	thick = lthick, $
	charthick = lthick, $
	title = 'IRAS 18588+3517', $
	xr = [0,9000], $
	yr = [-20, 2]

oploterror, ej, alog(nj / gj), alog(nj / gj) - alog((nj - nj_err) / gj), psym = 3, /lobar, thick = lthick
oploterror, ej, alog(nj / gj), -alog(nj / gj) + alog((nj + nj_err) / gj), psym = 3, /hibar, thick = lthick

plotsym2,1,3
oplot, ej([2]), alog(nj[2] / gj([2])), psym = 8

hotind = [1,3]
hotfit = mpfitexpr(expr, ej(hotind), alog(nj(hotind) / gj(hotind)), alog(nj_err / gj), $
	perror = perror, start, /quiet)

hotarr = fillarr(1,ej(1),ej(3))
warmarr = fillarr(1,ej(0),ej(1))
oplot, hotarr, hotfit(0) + hotfit(1) * hotarr, linestyle = 1, thick = lthick

fourarr = [0,1]
sub_yvalues = ej(fourarr) * hotfit(1) + hotfit(0)
sub_yvalues = sub_yvalues - sub_yvalues(0)
new_yvalues = alog(nj(fourarr) / gj(fourarr)) - sub_yvalues

oplot, ej(fourarr), new_yvalues, psym = 2

warmfit = mpfitexpr(expr, ej(fourarr), new_yvalues, replicate(stddev(0.1*new_yvalues),n_elements(new_yvalues)), $
	start, /quiet)
oplot, warmarr, warmfit(0) + warmarr*warmfit(1), linestyle = 2, thick = lthick

xyouts, ej+100, alog(nj / gj)+2, hlines, charsize = cs, charthick = lthick

xyouts, 2000, -2, 'T!Iwarm!N = '+string(-1d/warmfit(1), format = '(i4)')+' K', charsize = cs, charthick = lthick
xyouts, 5000, -8, 'T!Ihot!N  = '+string(-1d/hotfit(1), format = '(i4)')+' K', charsize = cs, charthick = lthick


; End PS output


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif


end
