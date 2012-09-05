function pahread, file

nlines = file_lines(file)
emptyarr = strarr(nlines)

openr, lun, file, /get_lun
readf, lun, emptyarr
close, lun & free_lun, lun

line_loc = where(strmid(emptyarr, 0, 22) eq  ' Integrated Line Flux:')
line_string = strmid(emptyarr(line_loc),33,10)
line = float(line_string)

err_loc = where(strmid(emptyarr, 0, 22) eq  ' Baseline Flux Density')
err_string = strmid(emptyarr(err_loc),59,11)
err = float(err_string)
base_string = strmid(emptyarr(err_loc),41,13)
baseline = float(base_string)

product = [line, err, baseline]
return, product

end

;;;;;;;

pro pahbatch

flist62 = '~/Astronomy/Research/Spitzer/OHM/lines/round2/lores/*pah62_lines.txt' 
flist11 = '~/Astronomy/Research/Spitzer/OHM/lines/round2/lores/*pah11_lines.txt' 

; Finding the line strength

files62 = file_search(flist62)
nfiles62 = n_elements(files62)
linearr62 = fltarr(nfiles62)
errarr62 = fltarr(nfiles62)
baseline62 = fltarr(nfiles62)
tagarr62 = strarr(nfiles62)

files11 = file_search(flist11)
nfiles11 = n_elements(files11)
linearr11 = fltarr(nfiles11)
errarr11 = fltarr(nfiles11)
baseline11 = fltarr(nfiles11)
tagarr11 = strarr(nfiles11)

for i = 0, nfiles62 - 1 do begin
	ptemp = pahread(files62(i))
	line_extract = ptemp(0)
	err_extract = ptemp(1)
	base_extract = ptemp(2)
	tagextract = strmid(files62(i),66,7)
	linearr62(i) = line_extract
	errarr62(i) = err_extract
	baseline62(i) = base_extract
	tagarr62(i) = tagextract
endfor

for i = 0, nfiles11 - 1 do begin
	ptemp = pahread(files11(i))
	line_extract = ptemp(0)
	err_extract = ptemp(1)
	base_extract = ptemp(2)
	tagextract = strmid(files11(i),66,7)
	linearr11(i) = line_extract
	errarr11(i) = err_extract
	baseline11(i) = base_extract
	tagarr11(i) = tagextract
endfor


linearr62 = linearr62 
errarr62 = errarr62 
ew62 = linearr62 / baseline62
pmarr62 = replicate(' +- ',nfiles62)
ewerr62 = ew62 * errarr62 / (baseline62)

linearr11 = linearr11 
errarr11 = errarr11 
ew11 = linearr11 / baseline11
pmarr11 = replicate(' +- ',nfiles11)
ewerr11 = ew11 * errarr11 / (baseline11)

ew62 = ew62 / 1d-3 & ew11 = ew11 / 1d-3
ewerr62 = ewerr62 / 1d-3 & ewerr11 = ewerr11 / 1d-3
print,''
print,'Fluxes of PAH detections [10^-3 um]:'
print,'           PAH 6.2 EW               PAH 11.3 EW'
print,''
print, transpose([[tagarr62],[string(ew62,format='(f10.3)')],[pmarr62],[string(ewerr62,format='(f6.3)')], $
	[string(ew11,format='(f10.3)')],[pmarr11],[string(ewerr11,format='(f6.3)')]])
print,''
print,'Avg. PAH 6.2 um flux [10^-3 um]: ',mean(ew62),' +- ',stddev(ew62)
print,'Avg. PAH 6.2 um error [10^-3 um]: ',mean(ewerr62),' +- ',stddev(ewerr62)
print,''
print,'Avg. PAH 11.3 um flux [10^-3 um]: ',mean(ew11),' +- ',stddev(ew11)
print,'Avg. PAH 11.3 um error [10^-3 um]: ',mean(ewerr11),' +- ',stddev(ewerr11)
print,''


stop
end


