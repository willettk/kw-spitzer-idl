;+
; NAME:
;       
;	DUSTY_ALLANALYZE
;
; PURPOSE:
;
;	
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
;       Written by K. Willett                Feb 10
;-

megatag = ohmdat('tag')
archtag = archdat('tag')
contag = condat('tag')

; Remove bad LR spectra

megatag = megatag[where(megatag ne 'mega034')]
contag = contag[where(contag ne 'control033')]

ohmtag = [megatag,transpose(archtag)]

ohmarr = fltarr(8, n_elements(ohmtag))
conarr = fltarr(8, n_elements(contag))

for i=0, n_elements(ohmtag) - 1 do begin
	
	tag, ohmtag[i], dirtag
	restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+ohmtag[i]+'_sphere_grid.sav'
	ohmarr[0:3,i] = [y_min, q_min, tauv_min, tdust]

	restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+ohmtag[i]+'_sphere_pahfitgrid.sav'
	ohmarr[4:7,i] = [y_min, q_min, tauv_min, tdust]

endfor

for i=0, n_elements(contag) - 1 do begin

	tag, contag[i], dirtag
	restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+contag[i]+'_sphere_grid.sav'
	conarr[0:3,i] = [y_min, q_min, tauv_min, tdust]

	restore, '~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+contag[i]+'_sphere_pahfitgrid.sav'
	conarr[4:7,i] = [y_min, q_min, tauv_min, tdust]

endfor

!p.multi=[0,4,2]

red  = fsc_color("Red")
blue = fsc_color("Blue")

plothist, ohmarr[0,*], /half, charsize = 2.0, xtitle='Y', xr=[0,1100], bin=10, title='Raw'
plothist, conarr[0,*], /half, xr=[0,1100], bin=10, /overplot, datacolor=blue
ver, mean(ohmarr[0,*]), linestyle=2, color=red
ver, mean(conarr[0,*]), linestyle=2, color=blue

plothist, ohmarr[1,*], /half, charsize = 2.0, xtitle='q', xr=[-1,4], bin=0.5, title='Raw'
plothist, conarr[1,*], /half, bin=0.5, /overplot, datacolor=blue
ver, mean(ohmarr[1,*]), linestyle=2, color=red
ver, mean(conarr[1,*]), linestyle=2, color=blue

plothist, ohmarr[2,*], /half, charsize = 2.0, xtitle='!7s!3!IV!N', xr=[0,500], bin=15, title='Raw'
plothist, conarr[2,*], /half, bin=15, /overplot, datacolor=blue
ver, mean(ohmarr[2,*]), linestyle=2, color=red
ver, mean(conarr[2,*]), linestyle=2, color=blue

plothist, ohmarr[3,*], /half, charsize = 2.0, xtitle='T!Idust!N [K]', xr=[0,150], bin=15, title='Raw'
plothist, conarr[3,*], /half, bin=10, /overplot, datacolor=blue
ver, mean(ohmarr[3,*]), linestyle=2, color=red
ver, mean(conarr[3,*]), linestyle=2, color=blue

plothist, ohmarr[4,*], /half, charsize = 2.0, xtitle='Y', xr=[0,1100], bin=10, title='PAHFIT'
plothist, conarr[4,*], /half, xr=[0,1100], bin=10, /overplot, datacolor=blue
ver, mean(ohmarr[4,*]), linestyle=2, color=red
ver, mean(conarr[4,*]), linestyle=2, color=blue

plothist, ohmarr[5,*], /half, charsize = 2.0, xtitle='q', xr=[-1,4], bin=0.5, title='PAHFIT'
plothist, conarr[5,*], /half, bin=0.5, /overplot, datacolor=blue
ver, mean(ohmarr[5,*]), linestyle=2, color=red
ver, mean(conarr[5,*]), linestyle=2, color=blue

plothist, ohmarr[6,*], /half, charsize = 2.0, xtitle='!7s!3!IV!N', xr=[0,500], bin=15, title='PAHFIT'
plothist, conarr[6,*], /half, bin=15, /overplot, datacolor=blue
ver, mean(ohmarr[6,*]), linestyle=2, color=red
ver, mean(conarr[6,*]), linestyle=2, color=blue

plothist, ohmarr[7,*], /half, charsize = 2.0, xtitle='T!Idust!N [K]', xr=[0,150], bin=15, title='PAHFIT'
plothist, conarr[7,*], /half, bin=10, /overplot, datacolor=blue
ver, mean(ohmarr[7,*]), linestyle=2, color=red
ver, mean(conarr[7,*]), linestyle=2, color=blue

!p.multi=[0,1,1]

; Print results to screen

print,''
print,'-----------------------------------------'

print,''
print,'Mean OHM Y (raw spectra): ',mean(ohmarr[0,*]),' +- ',stddev(ohmarr[0,*])
print,'Mean control Y (raw spectra): ',mean(conarr[0,*]),' +- ',stddev(conarr[0,*])

print,''
print,'Mean OHM Y (PAHFIT spectra): ',mean(ohmarr[4,*]),' +- ',stddev(ohmarr[4,*])
print,'Mean control Y (PAHFIT spectra): ',mean(conarr[4,*]),' +- ',stddev(conarr[4,*])

print,''
print,'-----------------------------------------'

print,''
print,'Mean OHM q (raw spectra): ',mean(ohmarr[1,*]),' +- ',stddev(ohmarr[1,*])
print,'Mean control q (raw spectra): ',mean(conarr[1,*]),' +- ',stddev(conarr[1,*])

print,''
print,'Mean OHM q (PAHFIT spectra): ',mean(ohmarr[5,*]),' +- ',stddev(ohmarr[5,*])
print,'Mean control q (PAHFIT spectra): ',mean(conarr[5,*]),' +- ',stddev(conarr[5,*])

print,''
print,'-----------------------------------------'

print,''
print,'Mean OHM tau_V (raw spectra): ',mean(ohmarr[2,*]),' +- ',stddev(ohmarr[2,*])
print,'Mean control tau_V (raw spectra): ',mean(conarr[2,*]),' +- ',stddev(conarr[2,*])

print,''
print,'Mean OHM tau_V (PAHFIT spectra): ',mean(ohmarr[6,*]),' +- ',stddev(ohmarr[6,*])
print,'Mean control tau_V (PAHFIT spectra): ',mean(conarr[6,*]),' +- ',stddev(conarr[6,*])

print,''
print,'-----------------------------------------'

print,''
print,'Mean OHM T_dust (raw spectra): ',mean(ohmarr[3,*]),' +- ',stddev(ohmarr[3,*])
print,'Mean control T_dust (raw spectra): ',mean(conarr[3,*]),' +- ',stddev(conarr[3,*])

print,''
print,'Mean OHM T_dust (PAHFIT spectra): ',mean(ohmarr[7,*]),' +- ',stddev(ohmarr[7,*])
print,'Mean control T_dust (PAHFIT spectra): ',mean(conarr[7,*]),' +- ',stddev(conarr[7,*])

stop

end
