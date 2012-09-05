function ir_lines, hr=hr, lr=lr
;+
; NAME:
;       
;	IR_LINES
;
; PURPOSE:
;
; 	Function to call IR lines for plotting
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
;	IDL> array = ir_lines(/lr)
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                
;	OH, not OH- 		Jan 09
;-

if keyword_set(hr) then begin

	lines = [8.025, 8.6, 8.99138, 9.665, 10.511, 11.3, 12.279, 12.372, 12.6, 12.814, 13.7, $
	14.0, 14.2, 14.322, 14.368, 15.0, 15.555, 16.4, 17.035, 17.1, 17.4, 17.934, 18.713, 24.318, 25.890, $
	25.988, 28.218, 33.481, 34.616, 34.815]

	line_id = ['H!I2!N S(4)', 'PAH', 'Ar III', 'H!I2!N S(3)', 'S IV', 'PAH', 'H!I2!N S(2)', 'HI 7-6', 'PAH', 'Ne II', $
	'C!I2!NH!I2!N', 'HCN', 'PAH', 'Ne V', 'Cl II', 'CO!I2!N', 'Ne III', 'PAH', 'H!I2!N S(1)', 'PAH', 'PAH', 'Fe II', 'S III', $
	'Ne V', 'O IV', 'Fe II', 'H!I2!N S(0)', 'S III', 'OH', 'Si II']

; Added OH feature; wavelength taken from Skinner et al (1997) (actually a blended set of two lambda doublets with separation of ~ 0.03 um)

endif else if keyword_set(lr) then begin

	lines = [5.511, 6.1088, 6.2, 6.909, 6.98, 7.7, 8.025, 8.6, 8.99138, 9.665, $
	10.511, 11.3, 12.279, 12.6, 12.814, $
	14.2,  14.322, 15.555, 16.4, 17.035, 17.1, 17.4, $
	18.713,  24.318, 25.890, 28.218]

	line_id = ['H!I2!N S(7)', 'H!I2!N S(6)', 'PAH', 'H!I2!N S(5)', 'Ar II', 'PAH', 'H!I2!N S(4)', 'PAH', 'Ar III', 'H!I2!N S(3)', $
	'S IV', 'PAH', 'H!I2!N S(2)', 'PAH', 'Ne II', $
	'PAH',  'Ne V', 'Ne III', 'PAH', 'H!I2!N S(1)', 'PAH', 'PAH', $
	'S III', 'Ne V', 'O IV', 'H!I2!N S(0)']
endif else begin
	message,'Must set /HR or /LR keywords'
	return,0
endelse



return,[[string(lines)],[line_id]]
end
