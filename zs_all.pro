pro zs_all, stop=stop
;+
; NAME:
;       
;	ZS_ALL
;
; PURPOSE:
;
;	Compute average delta V for the IR-optical and IR-OHM from the ensemble of lines detected for the entire sample (ie, look for a systematic bias)
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
;	IDL> .r zs_all
;
; NOTES:
;
;	Adapt for the entire OHM and control samples
;
; REVISION HISTORY
;       Written by K. Willett                
;	Add archived OHM data - Oct 08	
;-

c = 299792.458d

ohmlist = ohmdat('tag')
nolines = where(ohmlist eq 'mega016')
noind = setdifference(indgen(n_elements(ohmlist)),[nolines])
ohmlist = ohmlist(noind)
nohm = n_elements(ohmlist)

oharr_ohm = double(ohmdat('czoh'))
oherrarr_ohm = double(ohmdat('czoh_err'))
veloharr_ohm = oharr_ohm(noind)
zoharr_ohm = oharr_ohm(noind) / c
veloherrarr_ohm = oherrarr_ohm(noind)
zoherrarr_ohm = oherrarr_ohm(noind) / c

; Archived OHMs

archlist = transpose(archdat('tag'))
narch = n_elements(archlist)

oharr_arch = double(archdat('czoh'))
oherrarr_arch = double(archdat('czoh_err'))
veloharr_arch = oharr_arch
zoharr_arch = oharr_arch / c
veloherrarr_arch = oherrarr_arch
zoherrarr_arch = oherrarr_arch / c

; Combine

objlist = [ohmlist,archlist]
nobj = n_elements(objlist)

veloharr = [transpose(oharr_ohm), transpose(oharr_arch)]
zoharr = [zoharr_ohm, transpose(zoharr_arch)]
veloherrarr = [veloherrarr_ohm, transpose(veloherrarr_arch)]
zoherrarr = [zoherrarr_ohm,transpose(zoherrarr_arch)]

; Directories w/saved information

savdir_ohm = '~/Astronomy/Research/Spitzer/ohm/lines/hrsky/hires/saved/'
savdir_arch = '~/Astronomy/Research/Spitzer/archived/lines/hrsky/hires/saved/'
savdir_arch_nosky = '~/Astronomy/Research/Spitzer/archived/lines/nosky/hires/saved/'

; Archived objects with HR sky backgrounds

hrarch = 'arch'+string([5,10,20,29,31,34,35],format='(i03)')

savdir = '~/Astronomy/Research/Spitzer/ohm/lines/hrsky/hires/saved/'

dvoptarr = dblarr(nobj)
dvohmarr = dblarr(nobj)

print,''
print,'No v_OH for the following targets:'
print,''

for i=0, nobj - 1 do begin

	listtag = strmid(objlist[i],0,4)
	if listtag eq 'mega' then savdir = savdir_ohm $
		else if listtag eq 'arch' and where(objlist[i] eq hrarch) ne -1 then savdir = savdir_arch $
		else savdir = savdir_arch_nosky

	file = savdir+objlist[i]+'_zlines.sav'
	restore,file

	dvoptarr[i] = delv_opt
	if czoh ne 0 then begin
		dvohmarr[i] = delv_ohm		; Do not include objects without a measured OHM velocity
	endif else print,objlist[i]

endfor

remove_zeroes = where(dvohmarr ne 0)
dvohmarr = dvohmarr[remove_zeroes]

print,''
print,'mean del v IR-opt = ',mean(dvoptarr),' +- ',stddev(dvoptarr)
print,'mean del v IR-OHM = ',mean(dvohmarr),' +- ',stddev(dvohmarr)
print,''

if keyword_set(stop) then stop

end
