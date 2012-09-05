
;+
; NAME:
;       
;	DUSTY_BESTFITSAVE
;
; PURPOSE:
;
;	Re-format the save files from DUSTY_COMPARE_BATCH to include S_MIN, TAUV, TDUST
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

; Currently contain:

; errorgrid
; y_min
; q_min
; tauv_min
; tdust

; Need:
;
; s_min
; taudir
; taustring

pahfit_string=['pahfit','']
alltag = [transpose(ohmdat('tag')),transpose(archdat('tag')),transpose(condat('tag'))]

for j=0, 1 do begin
	for i=0, n_elements(alltag)-1 do begin
	
		tag, alltag[i], dirtag
		obj = alltag[i]
		
		file = '~/Astronomy/Research/Spitzer/'+dirtag+'/DUSTY/sphere_sav/'+obj+'_sphere_'+pahfit_string[j]+'grid.sav'
		
		restore, file
		
		; Find best-fit model
		
		minerror = where(errorgrid eq min(errorgrid[where(errorgrid gt 0.)]))
		if n_elements(minerror) gt 1 then print, 'Number of minima in errorgrid: ',n_elements(minerror)
		minerrorval = errorgrid_parse(minerror[0],/dusty,/alltau)
		
			Y_min      = minerrorval[0]
			q_min      = minerrorval[1]
			s_min_temp = minerrorval[2] + 1
		
			if s_min_temp gt 90 then begin
				s_min = s_min_temp - 90
				taudir = 'hightau'
				taustring = '_hightau'
				hightau = 1
			endif else begin
				s_min = s_min_temp
				taudir = 'regtau'
				taustring = ''
				hightau = 0
			endelse

		save, filename=file, $
			errorgrid, Y_min, q_min, tauv_min, tdust, s_min, taudir, taustring

	endfor	
endfor	

end
