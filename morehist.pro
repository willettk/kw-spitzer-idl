; More histograms of my various control samples to look for deviations
pro morehist, ps=ps

ohmdust = ohmdat('dtemp')
condust = condat('dtemp',sz=2)

ohmdust=float(ohmdust(0,*))
condust=float(condust(*,0))

ohmslopes = ohmdat('spindex')
conslopes = condat('spindex',sz=2)

os1 = float(ohmslopes(0,*))
os2 = float(ohmslopes(1,*))
cs1 = float(conslopes(*,0))
cs2 = float(conslopes(*,1))

!p.multi=[0,2,2]

plotname='~/Astronomy/Comps2/figures/morehist.ps'
if keyword_set(ps) then begin
	set_plot,'ps'
	device, filename = plotname, /color,/landscape
	cs = 1
	cthick = 2
	lthick = 2
endif else begin
	cs = 2
	cthick = 1
	lthick = 1
endelse


red = fsc_color("Red")
blue = fsc_color("Blue")

bs=13
plothist, ohmdust, /nodata, $
	xtitle='Dust temperature [K]', $
	ytitle='Frequency', $
	title='Dust temp. distribution', $
	xr=[20,110], $
	yr=[0,15], $
	charsize=cs, thick=lthick, charthick=cthick

plothist, ohmdust, color=red, bin = bs, /overplot, thick = lthick
plothist, condust, color=blue, bin = bs, /overplot, thick = lthick
ver,mean(ohmdust),color=red,linestyle=1, thick = lthick
ver,mean(condust),color=blue,linestyle=1, thick = lthick


bs = 1
plothist, cs1, /nodata, $
	xtitle='Spectral index 15-6!7l!3m', $
	ytitle='Frequency', $
	title='!7a!3!I15-6!N distribution', $
	xr=[-1,7], $
	yr=[0,15], $
	charsize=cs, thick=lthick, charthick=cthick

plothist, os1, color=red, bin = bs, /overplot, thick = lthick
plothist, cs1, color=blue, bin = bs, /overplot, thick = lthick
ver,mean(os1),color=red,linestyle=1, thick = lthick
ver,mean(cs1),color=blue,linestyle=1, thick = lthick

plothist, cs2, /nodata, $
	xtitle='Spectral index 30-20!7l!3m', $
	ytitle='Frequency', $
	title='!7a!3!I30-20!N distribution', $
	xr=[-1,8], $
	yr=[0,15], $
	charsize=cs, thick=lthick, charthick=cthick

plothist, os2, color=red, bin = bs, /overplot, thick = lthick
plothist, cs2, color=blue, bin = bs, /overplot, thick = lthick
ver,mean(os2),color=red,linestyle=1, thick = lthick
ver,mean(cs2),color=blue,linestyle=1, thick = lthick

fnames=ohmdat('tag')
cnames=condat('tag')
nf = n_elements(fnames)
nc = n_elements(cnames)
f60 = fltarr(nf)
f100 = fltarr(nf)
c60 = fltarr(nc)
c100 = fltarr(nc)

for i = 0,nf-1 do begin
	tag,fnames(i),dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+fnames(i)+'.sav'

	f60(i) = sed.iras(2)
	f100(i) = sed.iras(3)
endfor

for i = 0,nc-1 do begin
	tag,cnames(i),dirtag
	restore,'~/Astronomy/Research/Spitzer/'+dirtag+'/data/structures/'+cnames(i)+'.sav'

	c60(i) = sed.iras(2)
	c100(i) = sed.iras(3)
endfor

bothf = where(f60 ne 0 and f100 ne 0)
bothc = where(c60 ne 0 and c100 ne 0)

f60 = f60(bothf) & f100 = f100(bothf)
c60 = c60(bothc) & c100 = c100(bothc)

bs=0.15
plothist, f60/f100, /nodata, $
	xtitle='IRAS f!I60!N/f!I100!N', $
	ytitle='Frequency', $
	title='FIR color distribution', $
	xr=[-0.5,1.5], $
	yr=[0,15], $
	charsize=cs, thick=lthick, charthick=cthick

plothist, f60/f100, color=red, bin = bs, /overplot, thick = lthick
plothist, c60/c100, color=blue, bin = bs, /overplot, thick = lthick
ver,mean(f60/f100),color=red,linestyle=1, thick = lthick
ver,mean(c60/c100),color=blue,linestyle=1, thick = lthick


if keyword_set(ps) then begin
	device,/close
	set_plot,'x'
endif


stop
end
