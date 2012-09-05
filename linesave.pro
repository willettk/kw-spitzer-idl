;+
;
; ********* OBSOLETE ***************
;
; NAME:
;       LINESAVE 
;
; PURPOSE:
;	Read output files of line fluxes in Spitzer IRS spectra (as measured in SMART) and write them to an IDL structure.
;
; INPUTS:
;
;	ION - string detailing the line transition to measure (ex, 'neIII')
;
; OUTPUTS:
;
;	- IDL save file with line fluxes and errors for all objects. If the line was not detected, then flux and error are 0.
;
; KEYWORDS:
;
;	LORES - measure the fluxes from the (calibrated) lo-res spectra; default will measure the deblended hi-res lines
;
;	CON - save the line results for the non-masing control sample. Default operates on OHMs.
;
;	CSO - save the line results for the CSOs. Default operates on OHMs.
;
; EXAMPLE:
;	
;	IDL> linesave, 'neIII', /cso
;
; NOTES:
;
;	Line fluxes must be saved from SMART in specific directories for the files to be found. See location of OHM
;		files for examples. 
;
;	This routine is obsolete - use LINEHR.pro now. - KW, Aug 09
;
; REVISION HISTORY
;       Written by K. Willett                	Sep 2007
;	Added header, CSO keyword		Oct 2007
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

product = [line, err]
return, product

end

;;;;;;;

pro linesave, ion, lores = lores, con = con, cso = cso

; Set values for number of characters in directory paths

if keyword_set(con) then begin
	objs = 'control' 
	btag = 70
	ltag = 10
endif else if keyword_set(cso) then begin
	objs = 'cso'
	btag = 66
	ltag = 6
endif else begin
	objs = 'OHM'
	btag = 66
	ltag = 7
endelse
if keyword_set(lores) then dir = 'lores' else dir = 'hires'

; Directory in which line fluxes are stored as text files

if not keyword_set(cso) then flist = '~/Astronomy/Research/Spitzer/'+objs+'/lines/round2/'+dir+'/*'+ion+'_lines.txt' $
	else flist = '~/Astronomy/Research/Spitzer/'+objs+'/lines/round1/'+dir+'/*'+ion+'_lines.txt'

; Finding the line strength

files = file_search(flist)
nfiles = n_elements(files)
linearr = fltarr(nfiles)
errarr = fltarr(nfiles)
tagarr = strarr(nfiles)


for i = 0, nfiles - 1 do begin
	ptemp = line(files(i))
	line_extract = ptemp(0)
	err_extract = ptemp(1)
	tagextract = strmid(files(i),btag,ltag)
	linearr(i) = line_extract
	errarr(i) = err_extract
	tagarr(i) = tagextract
endfor

linearr = linearr / 1d-21
errarr = errarr / 1d-21
pmarr = replicate(' +- ',nfiles)

print,''
print,'Fluxes of line detections [10^-21 W/cm^2]:'
print, transpose([[tagarr],[string(linearr,format='(f10.2)')],[pmarr],[string(errarr,format='(f6.2)')]])
print,''
if n_elements(linearr) gt 1 then print,'Avg. flux [W/cm^2]: ',mean(linearr),' +- ',stddev(linearr)
if n_elements(linearr) gt 1 then print,'Avg. error [W/cm^2]: ',mean(errarr),' +- ',stddev(errarr)
print,''

; Name and path for IDL structure

if keyword_set(con) then begin
	conname = condat('tag')
	nf = n_elements(conname)
	vals = fltarr(nf,2)
	for i = 0,n_elements(linearr)-1 do begin
		ind = where(tagarr(i) eq conname)
		vals(ind,0) = linearr(i)
		vals(ind,1) = errarr(i)
	endfor

	if keyword_set(lores) then savename = '~/Astronomy/Research/Spitzer/'+objs+'/lines/structures/'+ion+'_con_lr.sav' else $
		savename = '~/Astronomy/Research/Spitzer/'+objs+'/lines/structures/'+ion+'_con.sav'
endif else if keyword_set(cso) then begin
	fnames = csodat('tag')
	nf = n_elements(fnames)
	vals = fltarr(nf,2)
	for i = 0,n_elements(linearr)-1 do begin
		ind = where(tagarr(i) eq fnames)
		vals(ind,0) = linearr(i)
		vals(ind,1) = errarr(i)
	endfor
	if keyword_set(lores) then savename = '~/Astronomy/Research/Spitzer/'+objs+'/lines/structures/'+ion+'_lr.sav' else $
		savename = '~/Astronomy/Research/Spitzer/'+objs+'/lines/structures/'+ion+'.sav'
endif else begin
	fnames = ohmdat('tag')
	nf = n_elements(fnames)
	vals = fltarr(nf,2)
	for i = 0,n_elements(linearr)-1 do begin
		ind = where(tagarr(i) eq fnames)
		vals(ind,0) = linearr(i)
		vals(ind,1) = errarr(i)
	endfor
	if keyword_set(lores) then savename = '~/Astronomy/Research/Spitzer/'+objs+'/lines/structures/'+ion+'_lr.sav' else $
		savename = '~/Astronomy/Research/Spitzer/'+objs+'/lines/structures/'+ion+'.sav'
endelse

; Save data in IDL structure with line fluxes and errors; units are 10^-21 W/cm^2

save, vals, filename=savename

end


