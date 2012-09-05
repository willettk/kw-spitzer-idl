pro stitch_avg, cso = cso, control = control, arch = arch, verbose = verbose, stop = stop, both = both
;+
; NAME:
;       
;	STITCH_AVG
;
; PURPOSE:
;
;	Compute average scaling fractions for LR module stitching and photometric scaling
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
;	CSO - 		find averages for CSOs	
;
;	CONTROL -	find averages for OHM control sample
;
;	ARCH - 		find averages for archived OHMs	
;
; EXAMPLE:
;
;	IDL> stitch_avg,/cso
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Aug 08
;	Added VERBOSE keyword - Dec 08
;-

dirtag = 'ohm'
if keyword_set(cso) then dirtag = 'cso'
if keyword_set(arch) then dirtag = 'archived'
if keyword_set(control) then dirtag = 'control' 

dir='~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/scaling_frac/'
photfile = '~/Astronomy/Research/Spitzer/'+dirtag+'/data/idl_sav/lr_frac.sav'

filelist = file_search(dir+'*sav')

if keyword_set(both) then begin
	dir1='~/Astronomy/Research/Spitzer/archived/data/idl_sav/scaling_frac/'
	photfile1 = '~/Astronomy/Research/Spitzer/archived/data/idl_sav/lr_frac.sav'
	dir2='~/Astronomy/Research/Spitzer/ohm/data/idl_sav/scaling_frac/'
	photfile2 = '~/Astronomy/Research/Spitzer/ohm/data/idl_sav/lr_frac.sav'

	filelist1 = file_search(dir1+'*sav')
	filelist2 = file_search(dir2+'*sav')

	filelist = [filelist1,filelist2]
endif

nfiles = n_elements(filelist)

ll2_ll1_arr = dblarr(nfiles)
sl1_ll2_arr = dblarr(nfiles)
sl2_sl1_arr = dblarr(nfiles)

for i = 0, nfiles - 1 do begin
	restore,filelist[i]
	ll2_ll1_arr[i] = ll2_ll1_frac
	sl1_ll2_arr[i] = sl1_ll2_frac
	sl2_sl1_arr[i] = sl2_sl1_frac

	if keyword_set(verbose) then begin

		f = strsplit(filelist[i],'/',/ex)
		ff = strsplit(f[n_elements(f)-1],'.',/ex)
		fff = ff[0]
		print,''
		print,fff
		print,'LL2 to LL1: ',mean(ll2_ll1_arr)
		print,'SL1 to SL2: ',mean(sl1_ll2_arr)
		print,'SL2 to SL1: ',mean(sl2_sl1_arr)

	endif
endfor

plot,ll2_ll1_arr, yr=[0,4]
oplot,sl1_ll2_arr, linestyle=1
oplot,sl2_sl1_arr, linestyle=2

restore,photfile

if keyword_set(both) then begin
	restore,photfile1
	fracarr1 = lrfracarr
	restore,photfile2
	fracarr2 = lrfracarr
	
	lrfracarr = [fracarr1,fracarr2]

endif

print,''
print,'LL2 to LL1 average = ',string(mean(ll2_ll1_arr),format='(f5.2)'),' +- ', string(stddev(ll2_ll1_arr),format='(f5.2)')
print,'SL1 to SL2 average = ',string(mean(sl1_ll2_arr),format='(f5.2)'),' +- ', string(stddev(sl1_ll2_arr),format='(f5.2)')
print,'SL2 to SL1 average = ',string(mean(sl2_sl1_arr),format='(f5.2)'),' +- ', string(stddev(sl2_sl1_arr),format='(f5.2)')
print,''
print,'Photometric scaling = ',string(mean(lrfracarr[where(lrfracarr ne 1.0)]),format='(f5.2)'),' +- ', string(stddev(lrfracarr[where(lrfracarr ne 1.0)]),format='(f5.2)')
print,''

if keyword_set(stop) then stop
end
