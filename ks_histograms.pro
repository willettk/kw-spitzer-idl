pro ks_histograms, ps=ps, wait = wait, stop = stop
;+
; NAME:
;       
;	KS_HISTOGRAMS
;
; PURPOSE:
;
;	Compute two-sided Kolmogorov-Smirnov stat. for mid-IR data and plot histograms
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
;       Written by K. Willett                Mar 08
;	Updated for both samples - Nov 08
;-

otag = ohmdat('tag')
atag = archdat('tag')
ctag = condat('tag')

match,otag,['mega034'],oa,ob    & goodo = setdifference(indgen(n_elements(otag)),oa)
match,ctag,['control033'],ca,cb & goodc = setdifference(indgen(n_elements(ctag)),ca)

!p.multi=[0,2,2]

red = fsc_color("Red")
blue = fsc_color("Blue")
lightblue = fsc_color("Cyan")
yellow = fsc_color("Yellow")
orange = fsc_color("Orange")
green = fsc_color("Dark Green")
purple = fsc_color("Purple")
color6 = fsc_color("Rosy Brown")

;####################
; SET 1
;####################

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/ks_histograms.ps', /color, /landscape
	cs = 1
	ls = 3
	defcolor = fsc_color("Black")
endif else begin
	cs = 1.5
	ls = 1
	defcolor = fsc_color("White")
endelse

; NIR spectral index

o = ohmdat('spindex')
o15 = transpose(float(o(0,*)))
o15 = o15[goodo]

a = archdat('spindex')
a15 = transpose(float(a(0,*)))

c = condat('spindex')
c15 = transpose(float(c(0,*)))
c15 = c15[goodc]

bs = 0.5

plothist,a15,/nodata,$
	xr=[-1,6], $
	yr=[0,20], $
	axiscolor=defcolor, $
	title='NIR index', $
	xtitle='!7a!3!I15-6!N', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,[o15,a15], bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,c15, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean([o15,a15]),color=red, linestyle = 1, thick = ls
ver,mean(c15),color=blue,linestyle=1, thick = ls

kstwo, [o15,a15],c15, D_nir, prob_nir
gauss_nir = sqrt(2d) * inverf(1d - prob_nir)

print,''
print,'D_KS    for NIR index: '+string(D_nir,format='(f7.3)')
print,'KS-prob for NIR index: '+string(prob_nir,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_nir,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean([o15,a15]),format='(f7.2)')+' +- '+string(stddev([o15,a15]),format='(f7.2)')
print,'Average value (con) :  '+string(mean([c15]),format='(f7.2)')+' +- '+string(stddev([c15]),format='(f7.2)')

xyouts,0.2,0.85,/normal,string(gauss_nir,format='(f3.1)')+'!7r!3',charsize=cs

; MIR spectral index

o = ohmdat('spindex')
o30 = transpose(float(o(1,*)))
o30_slope_err = ohmdat('spindexerr')
o30_err = transpose(float(o30_slope_err[1,*]))

a = archdat('spindex')
a30 = transpose(float(a(1,*)))
a30_slope_err = archdat('spindexerr')
a30_err = transpose(float(a30_slope_err[1,*]))

c = condat('spindex')
c30 = transpose(float(c(1,*)))
c30_slope_err = condat('spindexerr')
c30_err = transpose(float(c30_slope_err[1,*]))

bs = 1

plothist,a30,/nodata,$
	xr=[-1,8], $
	yr=[0,20], $
	axiscolor=defcolor, $
	title='MIR index', $
	xtitle='!7a!3!I30-20!N', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,[o30,a30], bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,c30, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean([o30,a30]),color=red, linestyle = 1, thick = ls
ver,mean(c30),color=blue,linestyle=1, thick = ls

kstwo, [o30,a30],c30, D_mir, prob_mir
gauss_mir = sqrt(2d) * inverf(1d - prob_mir)

print,''
print,'D_KS    for MIR index: '+string(D_mir,format='(f7.3)')
print,'KS-prob for MIR index: '+string(prob_mir,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_mir,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean([o30,a30]),format='(f7.2)')+' +- '+string(stddev([o30,a30]),format='(f7.2)')
print,'Average value (con) :  '+string(mean([c30]),format='(f7.2)')+' +- '+string(stddev([c30]),format='(f7.2)')

xyouts,0.7,0.85,/normal,string(gauss_mir,format='(f3.1)')+'!7r!3',charsize=cs

; PAH 6.2 um EW

o = ohmdat('pahfit62ew')
o62ew = transpose(float(o[0,*]))
o62ew = o62ew[where(o62ew ne 0.)]

a = archdat('pahfit62ew')
a62ew = transpose(float(a[0,*]))
a62ew = a62ew[where(a62ew ne 0.)]

c = condat('pahfit62ew')
c62ew = transpose(float(c[0,*]))
c62ew = c62ew[where(c62ew ne 0.)]

ohm_pahew62 = [o62ew,a62ew]
con_pahew62 = [c62ew]

bs = 0.3

plothist,ohm_pahew62,/nodata,$
	xr=[-3,3], $
	yr=[0,20], $
	axiscolor=defcolor, $
	title='PAH 6.2 EW ', $
	xtitle='log EW [!7l!3m]', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,ohm_pahew62, bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,con_pahew62, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean(ohm_pahew62),color=red, linestyle = 1, thick = ls
ver,mean(con_pahew62),color=blue,linestyle=1, thick = ls

kstwo, ohm_pahew62,con_pahew62, D_ew62, prob_ew62
gauss_ew62 = sqrt(2d) * inverf(1d - prob_ew62)

print,''
print,'D_KS    for 6.2 EW  : '+string(D_ew62,format='(f7.3)')
print,'KS-prob for 6.2 EW  : '+string(prob_ew62,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_ew62,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_pahew62),format='(f7.2)')+' +- '+string(stddev([ohm_pahew62]),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_pahew62),format='(f7.2)')+' +- '+string(stddev([con_pahew62]),format='(f7.2)')

xyouts,0.2,0.35,/normal,string(gauss_ew62,format='(f3.1)')+'!7r!3',charsize=cs

; PAH 11.3 um EW

o = ohmdat('pahfit11ew')
o11ew = transpose(float(o[0,*]))
o11ew = o11ew[where(o11ew ne 0.)]

a = archdat('pahfit11ew')
a11ew = transpose(float(a[0,*]))
a11ew = a11ew[where(a11ew ne 0.)]

c = condat('pahfit11ew')
c11ew = transpose(float(c[0,*]))
c11ew = c11ew[where(c11ew ne 0.)]

ohm_pahew11 = alog10([o11ew,a11ew])
con_pahew11 = alog10([c11ew])

bs = 0.3

plothist,ohm_pahew11,/nodata,$
	xr=[-3,3], $
	yr=[0,20], $
	axiscolor=defcolor, $
	title='PAH 11.3 EW ', $
	xtitle='log EW [!7l!3m]', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,ohm_pahew11, bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,con_pahew11, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean(ohm_pahew11),color=red, linestyle = 1, thick = ls
ver,mean(con_pahew11),color=blue,linestyle=1, thick = ls

kstwo, ohm_pahew11,con_pahew11, D_ew11, prob_ew11
gauss_ew11 = sqrt(2d) * inverf(1d - prob_ew11)

print,''
print,'D_KS    for 11.3 EW : '+string(D_ew11,format='(f7.3)')
print,'KS-prob for 11.3 EW : '+string(prob_ew11,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_ew11,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_pahew11),format='(f7.2)')+' +- '+string(stddev([ohm_pahew11]),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_pahew11),format='(f7.2)')+' +- '+string(stddev([con_pahew11]),format='(f7.2)')

xyouts,0.7,0.35,/normal,string(gauss_ew11,format='(f3.1)')+'!7r!3',charsize=cs

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(wait) then junk1=get_kbrd()

;####################
; SET 2
;####################

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/ks_histograms2.ps', /color, /landscape
	cs = 1
	ls = 2
	defcolor = fsc_color("Black")
endif else begin
	cs = 1.5
	ls = 1
	defcolor = fsc_color("White")
endelse

; PAH 6.2 um luminosity

o = ohmdat('pahfit62lum')
o62lum = transpose(float(o))
o62lum = o62lum[where(finite(o62lum) eq 1)]

a = archdat('pahfit62lum')
a62lum = transpose(float(a))
a62lum = a62lum[where(finite(a62lum) eq 1)]

c = condat('pahfit62lum')
c62lum = transpose(float(c))
c62lum = c62lum[where(finite(c62lum) eq 1)]

opah62lum = [o62lum,a62lum]
cpah62lum = [c62lum]

bs = 0.3

plothist,opah62lum,/nodata,$
	xr=[8,11], /xstyle, $
	yr=[0,20], $
	axiscolor=defcolor, $
	title='PAH 6.2 !7l!3m luminosity', $
	xtitle='log L!I7.3!N [L'+sunsymbol()+']', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist, opah62lum, bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist, cpah62lum, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean(opah62lum),color=red, linestyle = 2, thick = ls
ver,mean(cpah62lum),color=blue,linestyle=2, thick = ls

kstwo, opah62lum,cpah62lum, D_pah62lum, prob_pah62lum
gauss_pah62lum = sqrt(2d) * inverf(1d - prob_pah62lum)

print,''
print,'D_KS    for 6.2 luminosity: '+string(D_pah62lum,format='(f7.3)')
print,'KS-prob for 6.2 luminosity: '+string(prob_pah62lum,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_pah62lum,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(opah62lum),format='(f7.2)')+' +- '+string(stddev([opah62lum]),format='(f7.2)')
print,'Average value (con) :  '+string(mean(cpah62lum),format='(f7.2)')+' +- '+string(stddev([cpah62lum]),format='(f7.2)')

xyouts,0.2,0.85,/normal,string(gauss_pah62lum,format='(f3.1)')+'!7r!3',charsize=cs

; PAH 11.3 um luminosity

o = ohmdat('pahfit11lum')
o11lum = transpose(float(o))
o11lum = o11lum[where(finite(o11lum) eq 1)]

a = archdat('pahfit11lum')
a11lum = transpose(float(a))
a11lum = a11lum[where(finite(a11lum) eq 1)]

c = condat('pahfit11lum')
c11lum = transpose(float(c))
c11lum = c11lum[where(finite(c11lum) eq 1)]

opah11lum = [o11lum,a11lum]
cpah11lum = [c11lum]

bs = 0.3

plothist,opah11lum,/nodata,$
	xr=[8,11], /xstyle, $
	yr=[0,20], $
	axiscolor=defcolor, $
	title='PAH 11.3 !7l!3m luminosity', $
	xtitle='log L!I11.3!N [L'+sunsymbol()+']', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist, opah11lum, bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist, cpah11lum, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean(opah11lum),color=red, linestyle = 2, thick = ls
ver,mean(cpah11lum),color=blue,linestyle=2, thick = ls

kstwo, opah11lum,cpah11lum, D_pah11lum, prob_pah11lum
gauss_pah11lum = sqrt(2d) * inverf(1d - prob_pah11lum)

print,''
print,'D_KS    for 11.3 luminosity: '+string(D_pah11lum,format='(f7.3)')
print,'KS-prob for 11.3 luminosity: '+string(prob_pah11lum,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_pah11lum,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(opah11lum),format='(f7.2)')+' +- '+string(stddev([opah11lum]),format='(f7.2)')
print,'Average value (con) :  '+string(mean(cpah11lum),format='(f7.2)')+' +- '+string(stddev([cpah11lum]),format='(f7.2)')

xyouts,0.7,0.85,/normal,string(gauss_pah11lum,format='(f3.1)')+'!7r!3',charsize=cs

; Silicate depth 

o = ohmdat('sil')
osil_raw     = transpose(float(o[0,*]))
osil         = osil_raw[where(osil_raw ne 0.)]
osil_err_raw = transpose(float(o[1,*]))
osil_err     = osil_err_raw[where(osil_raw ne 0.)]

a = archdat('sil')
asil_raw     = transpose(float(a[0,*]))
asil         = asil_raw[where(asil_raw ne 0.)]
asil_err_raw = transpose(float(a[1,*]))
asil_err     = asil_err_raw[where(asil_raw ne 0.)]

c = condat('sil')
csil_raw     = transpose(float(c[0,*]))
csil         = csil_raw[where(csil_raw ne 0.)]
csil_err_raw = transpose(float(c[1,*]))
csil_err     = csil_err_raw[where(csil_raw ne 0.)]

ohm_silicate = [osil,asil]
con_silicate = csil

ohm_silicate_err = [osil_err,asil_err]
con_silicate_err = csil_err

bs = 0.5

plothist,ohm_silicate, /nodata, $
	axiscolor=defcolor, $
	xr=[-5,1], $
	yr=[0,20], $
	title='Silicate strength', $
	xtitle='Silicate strength', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,ohm_silicate,bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,con_silicate,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean(ohm_silicate),color=red, linestyle = 2, thick = ls
ver,mean(con_silicate),color=blue, linestyle = 2, thick = ls

kstwo, ohm_silicate,con_silicate, D_sil, prob_sil
gauss_sil = sqrt(2d) * inverf(1d - prob_sil)

print,''
print,'D_KS    for silicate strength: '+string(D_sil,format='(f7.3)')
print,'KS-prob for silicate strength: '+string(prob_sil,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_sil,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean(ohm_silicate),format='(f7.2)')+' +- '+string(stddev([ohm_silicate]),format='(f7.2)')
print,'Average value (con) :  '+string(mean(con_silicate),format='(f7.2)')+' +- '+string(stddev([con_silicate]),format='(f7.2)')

xyouts,0.2,0.35,/normal,string(gauss_sil,format='(f3.1)')+'!7r!3',charsize=cs

; Dust temperature

o = ohmdat('dtemp')
odtemp = transpose(float(o[0,*]))

a = archdat('dtemp')
adtemp = transpose(float(a[0,*]))

c = condat('dtemp')
cdtemp = transpose(float(c[0,*]))

bs = 10

plothist,adtemp,/nodata,$
	xr=[30,100], $
	yr=[0,25], $
	axiscolor=defcolor, $
	title='Dust temperature', $
	xtitle='T!Idust!N [K]', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,[odtemp,adtemp], bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,cdtemp, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean([odtemp,adtemp]),color=red, linestyle = 2, thick = ls
ver,mean(cdtemp),color=blue,linestyle=2, thick = ls

kstwo, [odtemp,adtemp],cdtemp, D_dtemp, prob_dtemp
gauss_dtemp = sqrt(2d) * inverf(1d - prob_dtemp)

print,''
print,'D_KS    for dust temperature: '+string(D_dtemp,format='(f7.3)')
print,'KS-prob for dust temperature: '+string(prob_dtemp,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_dtemp,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean([odtemp,adtemp]),format='(i4)')+' +- '+string(stddev([odtemp,adtemp]),format='(i4)')
print,'Average value (con) :  '+string(mean([cdtemp]),format='(i4)')+' +- '+string(stddev([cdtemp]),format='(i4)')

xyouts,0.7,0.35,/normal,string(gauss_dtemp,format='(f3.1)')+'!7r!3',charsize=cs

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(wait) then junk2=get_kbrd()

;####################
; SET 3
;####################

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/plots/ks_histograms3.ps', /color, /landscape
	cs = 1
	ls = 2
	defcolor = fsc_color("Black")
endif else begin
	cs = 1.5
	ls = 1
	defcolor = fsc_color("White")
endelse

; L_FIR 

o = ohmdat('lfir')
o_lfir = transpose(float(o))

a = archdat('lfir')
a_lfir = transpose(float(a))

c = condat('lfir')
c_lfir = transpose(float(c))

bs = 0.3

plothist, o_lfir, /nodata, $
	axiscolor=defcolor, $
	xr=[10,14], $
	yr=[0,25], $
	title='L!IFIR!N', $
	xtitle='log L!IFIR!N [L'+sunsymbol()+']', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,[o_lfir,a_lfir],bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,c_lfir,         bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean([o_lfir,a_lfir]),color=red, linestyle = 1, thick = ls
ver,mean(c_lfir),color=blue, linestyle = 1, thick = ls

kstwo, [o_lfir,a_lfir],c_lfir, D_lfir, prob_lfir
gauss_lfir = sqrt(2d) * inverf(1d - prob_lfir)

print,''
print,'D_KS    for L_FIR: '+string(D_lfir,format='(f7.3)')
print,'KS-prob for L_FIR: '+string(prob_lfir,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_lfir,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean([o_lfir,a_lfir]),format='(f7.2)')+' +- '+string(stddev([o_lfir,a_lfir]),format='(f7.2)')
print,'Average value (con) :  '+string(mean([c_lfir]),format='(f7.2)')+' +- '+string(stddev([c_lfir]),format='(f7.2)')

xyouts,0.2,0.85,/normal,string(gauss_lfir,format='(f3.1)')+'!7r!3',charsize=cs

; IRAS f60/f100

o = ohmdat('iras')
o60 =  float(o[2,*])
o100 = float(o[3,*])

a = archdat('iras')
a60 =  float(a[2,*])
a100 = float(a[3,*])

c = condat('iras')
c60 =  float(c[2,*])
c100 = float(c[3,*])

goodoindices = where(o60 gt 0 and o100 gt 0)
goodaindices = where(a60 gt 0 and a100 gt 0)
goodcindices = where(c60 gt 0 and c100 gt 0)

oratio = alog10(o60[goodoindices]/o100[goodoindices])
aratio = alog10(a60[goodaindices]/a100[goodaindices])
cratio = alog10(c60[goodcindices]/c100[goodcindices])

bs = 1d-1

plothist,oratio,/nodata,$
	xr=[-0.5,0.3], $
	yr=[0,25], $
	axiscolor=defcolor, $
	title='FIR flux ratio', $
	xtitle='log (f!I60!N/f!I100!n)', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	charthick = ls

plothist,[oratio,aratio], bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,cratio,          bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls
ver,mean([oratio,aratio]),color=red, linestyle = 2, thick = ls
ver,mean(cratio),color=blue,linestyle=2, thick = ls

kstwo, [oratio,aratio],cratio, D_iras, prob_iras
gauss_iras = sqrt(2d) * inverf(1d - prob_iras)

print,''
print,'D_KS    for f60/f100: '+string(D_iras,format='(f7.3)')
print,'KS-prob for f60/f100: '+string(prob_iras,format='(f7.3)')
print,'Gaussian probability:  '+string(gauss_iras,format='(f7.3)')+' sig'
print,'Average value (OHM) :  '+string(mean([oratio,aratio]),format='(f7.2)')+' +- '+string(stddev([oratio,aratio]),format='(f7.2)')
print,'Average value (con) :  '+string(mean([cratio]),format='(f7.2)')+' +- '+string(stddev([cratio]),format='(f7.2)')

xyouts,0.7,0.85,/normal,string(gauss_iras,format='(f3.1)')+'!7r!3',charsize=cs

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

; Make figure for Paper II

if keyword_set(ps) then begin
	set_plot,'ps'
	device,filename='~/Astronomy/Research/Spitzer/papers/ks_hist_II.ps', /color, /portrait, xs=18, ys=8
	cs = 1.2
	ls = 4
	defcolor = fsc_color("Black")
endif else begin
	cs = 1.5
	ls = 1
	defcolor = fsc_color("White")
endelse

bs = 1.

!p.multi=[0,2,1]

red = fsc_color("Red")
blue = fsc_color("Blue")

plothist,[o30,a30],/nodata,$
	xr=[-1,8], $
	yr=[0,25], $
	/xstyle, /ystyle, $
	axiscolor=defcolor, $
	xtitle='!7a!3!I30-20!N', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	xthick = ls, $
	ythick = ls, $
	charthick = ls

plothist,[o30,a30], bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,c30,       bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls

meanerr2, [o30, a30], [o30_err, a30_err], ohm30_mean, ohm_30_sig
meanerr2, c30, c30_err, con30_mean, con_30_sig

ver, ohm30_mean,color=red, linestyle = 2, thick = ls
ver, con30_mean,color=blue,linestyle=2, thick = ls

bs = 0.5

plothist,ohm_silicate,/nodata,$
	xr=[-4,1], $
	yr=[0,16], $
	/xstyle, /ystyle, $
	axiscolor=defcolor, $
	xtitle='S!I9.7!N', $
	ytitle = 'Frequency', $
	charsize = cs, $
	thick = ls, $
	xthick = ls, $
	ythick = ls, $
	charthick = ls

plothist,ohm_silicate, bin=bs,/halfbin,datacolor=red, /overplot, thick = ls
plothist,con_silicate, bin=bs,/halfbin,datacolor=blue,/overplot, thick = ls

meanerr2, ohm_silicate, ohm_silicate_err, ohm_sil_weighted_mean, ohm_sil_weighted_sigma
meanerr2, con_silicate, con_silicate_err, con_sil_weighted_mean, con_sil_weighted_sigma

ver,ohm_sil_weighted_mean,color=red, linestyle = 2, thick = ls
ver,con_sil_weighted_mean,color=blue,linestyle=2, thick = ls

if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif

if keyword_set(stop) then stop
end
