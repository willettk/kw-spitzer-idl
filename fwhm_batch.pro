;+
; NAME:
;       
;	FWHM_BATCH
;
; PURPOSE:
;
;	Compute FWHMs of line detections and upper limits on flux for non-detections
;
; INPUTS:
;
;	ION - 		string with ion name to compute (eg, 'neII')
;
; OUTPUTS:
;
;	Prints list of FWHM and limits to screen
;
; KEYWORDS:
;
;	CONTROL - 	searches control sample data
;
;	ARCH - 		searches archived OHM data
;
;	CSO - 		searches CSO data
;
;	MEGA - 		searches Darling OHM data
;
;
; EXAMPLE:
;
;	IDL> fwhm_batch,'neII',/mega
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett       May 08         
;-

function linelim, file, ion, linecen

; Read in data

readcol, file, skipline = 1, order, wave, flux, err, bit, format = 'i,f,f,f,i', /silent

if strmid(file,29,3) eq 'con' then begin
	btag = 71
	etag = 10
endif else begin
	btag = 67
	etag = 7
endelse

tag = strmid(file,btag,etag)
targets,tag,r,o

waveind = sort(wave)
wavetemp = wave(waveind) & fluxtemp = flux(waveind)
wave = wavetemp / (1d + r) & flux = fluxtemp

plot, wave, flux, psym=10, xr=[linecen-2,linecen+2],yr=[0,0.5],/xstyle
ver, linecen

width = 0.2

bind = closeto(wave,linecen-width)
eind = closeto(wave,linecen+width)

if bind eq eind then sigflux = 0 else begin

	ver,wave(bind), linestyle = 1
	ver,wave(eind), linestyle = 1
	
	meanflux = mean(flux(bind:eind))
	sigflux  = stddev(flux(bind:eind))
	
	;print,'Mean: ',meanflux
	;print,'Sigma: ',sigflux
	
	hor, meanflux + 3*sigflux, color=fsc_color("Yellow")
	hor, meanflux - 3*sigflux, color=fsc_color("Yellow")
	
endelse

return, sigflux

end

;;;;;;;

pro fwhm_batch, ion, control = control, arch = arch, cso = cso, mega = mega, lores = lores

; Control sample vs. OHMs

if keyword_set(control) then begin
	objs = 'control' 
	btag = 70
	ltag = 10
endif else begin
	objs = 'OHM'
	btag = 66
	ltag = 7
endelse

; Finding the average linewidth

if keyword_set(lores) then res = 'lores' else res = 'hires'

flist = '~/Astronomy/Research/Spitzer/'+objs+'/lines/round4/'+res+'/*'+ion+'_lines.txt'
files = file_search(flist)
nfiles = n_elements(files)
fwhmarr = fltarr(nfiles)
tagarr = strarr(nfiles)


for i = 0, nfiles - 1 do begin
	fwhm_extract = fwhm(files(i))
	tagextract = strmid(files(i),btag,ltag)
	fwhmarr(i) = fwhm_extract
	tagarr(i) = tagextract
endfor

print,'FWHMs of line detections:'
print, transpose([[tagarr],[string(fwhmarr)]])
print,''
print,'Avg. linewidth [um]: ',mean(fwhmarr),' +- ',stddev(fwhmarr)
print,''

; Finding the S/N in the region

alltags = ohmdat('tag')

limtags = strarr(n_elements(alltags))
for i = 0, n_elements(alltags)-1 do begin
	hp = where(alltags(i) eq tagarr)
	if hp eq -1 then limtags(i) = alltags(i)
endfor

junk1 = where(limtags ne '',count)
if count ne 0 then begin

	limtags = limtags(where(limtags ne ''))
	nltags = n_elements(limtags)
	
	ionlist = ['neII','neIII','sIII','sIV','oIV','arIII','h2s0','h2s1_sh','h2s1_lh','h2s1','h2s2','h2s3','h2s5','h2s7']
	linelist = [12.814,15.555,18.713,10.511,25.890,8.993,28.221,17.035,17.035,17.035,12.279,9.665,6.9,5.51]
	linecen = linelist(where(ion eq ionlist))
	linecen = linecen(0)
	
	temp2 = replicate('~/Astronomy/Research/Spitzer/'+objs+'/data/idl_spectra/calibrated/coadd/',nltags)
	if keyword_set(lores) then temp3 = replicate('_lr_cal.tbl',nltags) else	if keyword_set(lh) then temp3 = replicate('_lh_cal.tbl',nltags) else temp3 = replicate('_sh_cal.tbl',nltags)
	flist2 = temp2+limtags+temp3
		
	avgfwhm = mean(fwhmarr)
	
	jyconv = 1d-23
	micron2hz = 3d14 / (linecen)^2		; dnu = c / lambda^2 dlambda
	
	limarr = fltarr(nltags)
	upperlimit = fltarr(nltags)
	for i = 0, nltags - 1 do begin
		sigflux = linelim(flist2(i), ion, linecen)
		limarr(i) = sigflux
		upperlimit(i) = 2.26 * (avgfwhm * micron2hz) * (sigflux * jyconv) / 1d7
	endfor
	
	print,'Upper limits on line flux [10^-21 W/cm^2]: '
	print, transpose([[limtags],[string(upperlimit/1d-21)]])

endif else print, 'Lines detected for all targets in sample'

end
