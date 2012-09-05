;+
; NAME:
;       
;	LINEBATCH
;
; PURPOSE:
;
;	Record line fluxes and fits from the SMART text outputs
;
; INPUTS:
;
;	ion - 	string designating the ion or molecule to measure
;
; OUTPUTS:
;
;	Prints out list of line fluxes for each object in which it has been measured
;
; KEYWORDS:
;
;
;
; EXAMPLE:
;
;	IDL> linebatch, 'neIII'
;
; NOTES:
;
;
;
; REVISION HISTORY
;       Written by K. Willett                Apr 08
;	Filepaths now obsolete		KW, Jul 08
;-

function line, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

line_loc = where(strmid(emptyarr, 0, 11) eq  ' Line Flux:')
line_string = strmid(emptyarr(line_loc),33,10)
line = float(line_string)

err_string = strmid(emptyarr(line_loc),47,11)
err = float(err_string)

linecen_loc = where(strmid(emptyarr, 0, 13) eq  ' Line Center:')
linecen_string = strmid(emptyarr(linecen_loc),33,8)
linecen = float(linecen_string)

lrest_loc = where(strmid(emptyarr, 0, 10) eq  'Wavelength')
lrest_string = strsplit(emptyarr(lrest_loc+2),/ex)
lrest = float(lrest_string(0))

product = [line,err,linecen,lrest]
return, product

end

;;;;;;;

function ionex, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

erest_loc = where(strmid(emptyarr, 0, 10) eq  'Wavelength')
erest_string = strsplit(emptyarr(erest_loc+2),/ex)

element = erest_string(2)
transition = strjoin(erest_string(3:n_elements(erest_string)-1))

return, [element,transition]

end

;;;;;;;

function tagex, file

junk1 = strsplit(file,'/',/ex)
junk1 = junk1(n_elements(junk1) - 1)
junk2 = strsplit(junk1,'.',/ex)
junk2 = junk2(0)
junk3 = strsplit(junk2,'_',/ex)

tag = junk3(0)
ion = junk3(1)

return, tag
end

;;;;;;;

pro linebatch, ion, linearr, errarr, linecenarr, lrestarr, lores = lores, stop = stop

if keyword_set(lores) then res = 'lores' else res = 'hires'

tag = 'ohm'

flist = '~/Astronomy/Research/Spitzer/'+tag+'/lines/round3/'+res+'/*'+ion+'_*txt'

; Finding the line strength

files = file_search(flist)
nfiles = n_elements(files)

linearr = fltarr(nfiles)
errarr = fltarr(nfiles)
linecenarr = fltarr(nfiles)
lrestarr = fltarr(nfiles)
zarr = fltarr(nfiles)
tagarr = strarr(nfiles)


for i = 0, nfiles - 1 do begin
	ptemp = line(files(i))
	linearr(i) = ptemp(0)
	errarr(i) = ptemp(1)
	linecenarr(i) = ptemp(2)
	lrestarr(i) = ptemp(3)
	tagarr(i) = tagex(files(i))
	targets,tagarr(i),r,o
	zarr(i) = r
endfor

linearr = linearr / 1d-21
errarr = errarr / 1d-21
pmarr = replicate(' +- ',nfiles)

ionarr = ionex(files(0))
element = ionarr(0)
transition = ionarr(1)

; Check to make sure fitted lines are all the same

testarr = where(lrestarr eq lrestarr(0)) - lindgen(n_elements(lrestarr)) 
if abs(max(testarr)) ne 0 then begin
	print, 'Not all lines are the same'
	junkind = where(abs(max(testarr)) ne 0)
	print,'Offending object: '+tagarr(junkind)
	stop
endif

; Find rest frequency of line

zshift = linecenarr / lrestarr - 1d

; Print out the tag, flux, flux error, measured line center, expected line center, and shift

print,''
print,'Ion:              '+element+' '+transition
print,'lambda_0 = '+string(lrestarr(0))+' um'
print,''
print,'Tag          Flux [10^-21 W/cm^2]      New z      Old z    delta z'
print, transpose([[tagarr],[string(linearr,format='(f10.2)')],[pmarr],[string(errarr,format='(f6.2)')], $
	[string(zshift,format='(f13.4)')],[string(zarr,format='(f10.4)')],$
	[string(zarr-zshift,format='(f10.4)')]])
print,''
if n_elements(linearr) gt 1 then begin
	print,'Avg. flux  [10^-21 W/cm^2]: ',mean(linearr),' +- ',stddev(linearr)
	print,'Avg. error [10^-21 W/cm^2]: ',mean(errarr),' +- ',stddev(errarr)
	print,''
endif

if keyword_set(stop) then stop

end
